#!/usr/bin/env python3
"""Run the review's literature search and build one deduplicated file for screening.

Sources:
  PubMed       the protocol's query through NCBI E-utilities: hit count, records
               and the run date
  seeds        checks that the seed papers are retrieved; for any that are not,
               checks whether PubMed indexes them at all (query problem or
               database problem)
  OpenAlex     backward and forward citation chasing from the seed papers, plus a
               preprint search that approximates searching bioRxiv, ChemRxiv and
               arXiv directly
  RIS imports  exports from Scopus and Web of Science, whose queries you run by
               hand with your institution's access

Outputs in --out: pubmed.csv, citations.csv, preprints.csv, combined.csv,
combined.ris (import into Rayyan or Zotero for screening) and run_log.md (the
dates and counts the protocol asks you to record).

OpenAlex has expected an API key since February 2026: without one there is a
small daily allowance and some request styles are refused. A free key (create
an account at openalex.org) goes in --openalex-key or the OPENALEX_API_KEY
environment variable. Results are read page by page, which works either way,
and a failing source is reported without stopping the run.

Uses only the Python standard library.

Examples:
  python3 review/review_search.py --email you@example.org --out review/run
  python3 review/review_search.py --email you@example.org --out review/run --import-ris scopus.ris --import-ris wos.ris
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OPENALEX = "https://api.openalex.org"
HERE = Path(__file__).resolve().parent
FIELDS = ["sources", "year", "title", "authors", "journal", "doi", "pmid", "abstract"]
PREPRINT_TERMS = "psychoactive mass spectrometry"
OPENALEX_KEY = ""  # set from --openalex-key or OPENALEX_API_KEY in main()


def oa(params: dict, email: str) -> dict:
    """OpenAlex query parameters with the contact email and, when given, the API key."""
    out = {**params, "mailto": email}
    if OPENALEX_KEY:
        out["api_key"] = OPENALEX_KEY
    return out


def http_get(url: str, params: dict | None = None, retries: int = 3) -> str | None:
    """GET with retries for server errors; None for 404; immediate error for other client errors."""
    if params:
        url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "nps-openworld-review/0.1"})
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except urllib.error.HTTPError as err:
            if err.code == 404:
                return None
            if 400 <= err.code < 500 and err.code != 429:
                raise
            if attempt == retries - 1:
                raise
        except urllib.error.URLError:
            if attempt == retries - 1:
                raise
        time.sleep(2 * (attempt + 1))
    return None


def norm(text: str | None) -> str:
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


# ---------------------------------------------------------------- PubMed

def pubmed_search(term: str, email: str, retmax: int = 10000) -> tuple[int, list[str], str]:
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax,
              "tool": "nps-openworld-review", "email": email}
    result = json.loads(http_get(f"{EUTILS}/esearch.fcgi", params))["esearchresult"]
    time.sleep(0.4)  # NCBI allows three requests a second without an API key
    return int(result["count"]), result.get("idlist", []), result.get("querytranslation", "")


def parse_pubmed_xml(xml: str) -> list[dict]:
    records = []
    for article in ET.fromstring(xml).iter("PubmedArticle"):
        citation = article.find("MedlineCitation")
        art = citation.find("Article") if citation is not None else None
        if art is None:
            continue
        title_el = art.find("ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        abstract = " ".join("".join(x.itertext()).strip() for x in art.iter("AbstractText"))
        year = art.findtext("Journal/JournalIssue/PubDate/Year") or \
            (art.findtext("Journal/JournalIssue/PubDate/MedlineDate") or "")[:4]
        authors = []
        for author in art.findall("AuthorList/Author"):
            last, initials = author.findtext("LastName"), author.findtext("Initials")
            if last:
                authors.append(f"{last}, {initials}" if initials else last)
            elif author.findtext("CollectiveName"):
                authors.append(author.findtext("CollectiveName"))
        doi = next((e.text.strip() for e in art.findall("ELocationID")
                    if e.get("EIdType") == "doi" and e.text), "")
        if not doi:  # the article's own IDs only, not those in its reference list
            doi = next((e.text.strip() for e in article.findall("PubmedData/ArticleIdList/ArticleId")
                        if e.get("IdType") == "doi" and e.text), "")
        records.append({"source": "pubmed", "pmid": citation.findtext("PMID", ""), "doi": doi,
                        "title": title, "authors": "; ".join(authors), "year": year,
                        "journal": art.findtext("Journal/Title", ""), "abstract": abstract})
    return records


def pubmed_fetch(pmids: list[str], email: str) -> list[dict]:
    records = []
    for i in range(0, len(pmids), 200):
        xml = http_get(f"{EUTILS}/efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids[i:i + 200]),
                                                 "retmode": "xml", "tool": "nps-openworld-review",
                                                 "email": email})
        records += parse_pubmed_xml(xml)
        time.sleep(0.4)
    return records


def split_blocks(query: str) -> list[str]:
    """Top-level AND blocks of a PubMed query (AND inside parentheses is left alone)."""
    blocks, depth, start, i = [], 0, 0, 0
    while i < len(query):
        ch = query[i]
        depth += ch == "("
        depth -= ch == ")"
        if depth == 0 and query[i:i + 5].upper() == " AND ":
            blocks.append(query[start:i].strip())
            start, i = i + 5, i + 5
            continue
        i += 1
    blocks.append(query[start:].strip())
    return [b for b in blocks if b]


def seed_pmid(seed: dict, email: str) -> str:
    if seed.get("pmid"):
        return seed["pmid"]
    term = f"{seed['doi']}[doi]" if seed.get("doi") else f"\"{seed['title_fragment']}\"[ti]"
    ids = pubmed_search(term, email)[1]
    return ids[0] if ids else ""


def failing_blocks(seed: dict, query: str, email: str) -> list[str]:
    """The query blocks that do not match a seed PubMed indexes (empty if it cannot be found)."""
    pmid = seed_pmid(seed, email)
    if not pmid:
        return []
    return [block for block in split_blocks(query)
            if pubmed_search(f"({block}) AND {pmid}[pmid]", email)[0] == 0]


def try_terms(query: str, base_count: int, seed_rows: list[dict], terms: list[str], email: str) -> tuple[int, list[dict]]:
    """For each candidate term added to the block that loses a seed: PubMed records and seeds recovered."""
    blocks = split_blocks(query)
    targets = []
    for s in seed_rows:
        for block in filter(None, (s.get("missed_by") or "").split(" | ")):
            if block in blocks:
                targets.append((s["label"], seed_pmid(s, email), blocks.index(block)))
    if not targets:
        return -1, []
    index = targets[0][2]

    def widened(extra: str) -> str:
        block = blocks[index]
        return block[:-1] + f" OR {extra})" if block.endswith(")") else f"({block} OR {extra})"

    rows = []
    for label, extra in [(t, t) for t in terms] + [(f"all {len(terms)} together", " OR ".join(terms))]:
        block = widened(extra)
        count = pubmed_search(" AND ".join(blocks[:index] + [block] + blocks[index + 1:]), email, retmax=0)[0]
        recovered = [name for name, pmid, i in targets
                     if i == index and pmid and pubmed_search(f"({block}) AND {pmid}[pmid]", email, retmax=0)[0]]
        rows.append({"term": label, "records": count, "extra": count - base_count, "recovers": ", ".join(recovered)})
    return index, rows


def check_seeds(seeds: list[dict], records: list[dict], email: str, query: str = "") -> list[dict]:
    dois = {r["doi"].lower() for r in records if r.get("doi")}
    pmids = {r["pmid"] for r in records if r.get("pmid")}
    titles = [norm(r["title"]) for r in records]
    rows = []
    for seed in seeds:
        fragment = norm(seed.get("title_fragment"))
        found = bool((seed.get("doi") and seed["doi"].lower() in dois)
                     or (seed.get("pmid") and seed["pmid"] in pmids)
                     or (fragment and any(fragment in t for t in titles)))
        indexed = "yes"
        if not found:
            if seed.get("doi"):
                term = f"{seed['doi']}[doi]"
            elif seed.get("pmid"):
                term = f"{seed['pmid']}[pmid]"
            else:
                term = f"\"{seed['title_fragment']}\"[ti]"
            indexed = "yes" if pubmed_search(term, email)[0] else "no"
        missed_by = ""
        if not found and indexed == "yes" and query:
            missed_by = " | ".join(failing_blocks(seed, query, email))
        rows.append({**seed, "retrieved": "yes" if found else "no", "in_pubmed": indexed, "missed_by": missed_by})
    return rows


# ---------------------------------------------------------------- OpenAlex

def abstract_from_index(inverted: dict | None) -> str:
    if not inverted:
        return ""
    position = {i: word for word, places in inverted.items() for i in places}
    return " ".join(position[i] for i in sorted(position))


def openalex_record(work: dict, source: str) -> dict:
    location = work.get("primary_location") or {}
    return {
        "source": source,
        "pmid": (work.get("ids", {}).get("pmid") or "").rstrip("/").rsplit("/", 1)[-1],
        "doi": (work.get("doi") or "").replace("https://doi.org/", ""),
        "title": work.get("display_name") or work.get("title") or "",
        "authors": "; ".join(a["author"]["display_name"] for a in work.get("authorships", []) if a.get("author")),
        "year": str(work.get("publication_year") or ""),
        "journal": ((location.get("source") or {}).get("display_name")) or "",
        "abstract": abstract_from_index(work.get("abstract_inverted_index")),
    }


def openalex_list(params: dict, email: str, cap: int = 2000, per_page: int = 200) -> list[dict]:
    """All results of a works query, read page by page (OpenAlex refuses keyless cursor paging)."""
    works: list[dict] = []
    page = 1
    while len(works) < cap and page * per_page <= 10000:
        data = json.loads(http_get(f"{OPENALEX}/works", oa({**params, "per-page": per_page, "page": page}, email)))
        results = data.get("results", [])
        works += results
        if len(results) < per_page:
            break
        page += 1
        time.sleep(0.2)
    return works[:cap]


def openalex_seed(seed: dict, email: str) -> dict | None:
    """Find a seed in OpenAlex by DOI, falling back to a title search that must contain the fragment."""
    if seed.get("doi"):
        text = http_get(f"{OPENALEX}/works/doi:{seed['doi']}", oa({}, email))
        if text:
            return json.loads(text)
    fragment = seed.get("title_fragment")
    if fragment:
        text = http_get(f"{OPENALEX}/works", oa({"search": fragment, "per-page": 5}, email))
        for work in json.loads(text).get("results", []) if text else []:
            if norm(fragment) in norm(work.get("display_name")):
                return work
    return None


def openalex_by_ids(ids: list[str], email: str) -> list[dict]:
    works = []
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        try:
            works += openalex_list({"filter": "openalex:" + "|".join(chunk)}, email)
        except urllib.error.HTTPError:  # fall back to one request per work
            for work_id in chunk:
                text = http_get(f"{OPENALEX}/works/{work_id}", oa({}, email))
                if text:
                    works.append(json.loads(text))
                time.sleep(0.1)
    return works


def chase_citations(seeds: list[dict], email: str) -> tuple[list[dict], list[dict], list[str], list[dict]]:
    found, missing = [], []
    for seed in seeds:
        work = openalex_seed(seed, email)
        (found.append(work) if work else missing.append(seed["label"]))
    ref_ids = sorted({r.rsplit("/", 1)[-1] for w in found for r in w.get("referenced_works", [])})
    backward = [openalex_record(w, "backward citation") for w in openalex_by_ids(ref_ids, email)]
    forward = []
    for work in found:
        forward += [openalex_record(w, "forward citation")
                    for w in openalex_list({"filter": f"cites:{work['id'].rsplit('/', 1)[-1]}"}, email)]
    return backward, forward, missing, [openalex_record(w, "seed") for w in found]


def preprint_search(email: str, from_year: int) -> list[dict]:
    works = openalex_list({"search": PREPRINT_TERMS,
                           "filter": f"type:preprint,from_publication_date:{from_year}-01-01"}, email, cap=500)
    return [openalex_record(w, "preprint search") for w in works]


# ---------------------------------------------------------------- RIS

def read_ris(path: Path, source: str) -> list[dict]:
    records, current = [], {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        match = re.match(r"^([A-Z][A-Z0-9])  -\s?(.*)$", line)
        if not match:
            continue
        tag, value = match.group(1), match.group(2).strip()
        if tag == "ER":
            if current:
                first = lambda *tags: next((current[t][0] for t in tags if current.get(t)), "")
                records.append({"source": source, "pmid": first("AN") if first("AN").isdigit() else "",
                                "doi": first("DO").replace("https://doi.org/", ""),
                                "title": first("TI", "T1"), "authors": "; ".join(current.get("AU", []) + current.get("A1", [])),
                                "year": first("PY", "Y1")[:4], "journal": first("T2", "JO", "JF", "J2"),
                                "abstract": first("AB", "N2")})
            current = {}
            continue
        current.setdefault(tag, []).append(value)
    return records


def write_ris(path: Path, records: list[dict]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for r in records:
            fh.write("TY  - JOUR\n")
            fh.write(f"TI  - {r['title']}\n")
            for author in filter(None, (a.strip() for a in r["authors"].split(";"))):
                fh.write(f"AU  - {author}\n")
            for tag, key in (("PY", "year"), ("JO", "journal"), ("DO", "doi"), ("AN", "pmid"), ("AB", "abstract")):
                if r.get(key):
                    fh.write(f"{tag}  - {r[key]}\n")
            if r.get("doi"):
                fh.write(f"UR  - https://doi.org/{r['doi']}\n")
            fh.write(f"N1  - found by: {r['sources']}\n")
            fh.write("ER  - \n\n")


# ---------------------------------------------------------------- merge and output

def combine(groups: list[list[dict]], from_year: int) -> tuple[list[dict], int]:
    """Deduplicate by DOI, then by normalised title; drop records dated before from_year."""
    out: list[dict] = []
    by_doi: dict[str, int] = {}
    by_title: dict[str, int] = {}
    for records in groups:
        for r in records:
            doi, title = (r.get("doi") or "").lower(), norm(r.get("title"))[:150]
            index = by_doi.get(doi) if doi else None
            if index is None and title:
                index = by_title.get(title)
            if index is None:
                out.append({**r, "sources": r["source"]})
                index = len(out) - 1
            else:
                merged = out[index]
                merged["sources"] = "; ".join(sorted(set(merged["sources"].split("; ")) | {r["source"]}))
                for field in ("pmid", "doi", "abstract", "journal", "authors", "year"):
                    if not merged.get(field) and r.get(field):
                        merged[field] = r[field]
            if doi:
                by_doi[doi] = index
            if title:
                by_title[title] = index
    kept = [r for r in out if not (r.get("year", "")[:4].isdigit() and int(r["year"][:4]) < from_year)]
    return kept, len(out) - len(kept)


def write_csv(path: Path, records: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for r in records:
            writer.writerow({**r, "sources": r.get("sources", r.get("source", ""))})


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True, help="your email, sent to NCBI and OpenAlex as their usage policies ask")
    ap.add_argument("--query-file", type=Path, default=HERE / "pubmed_query.txt")
    ap.add_argument("--seeds", type=Path, default=HERE / "seeds.csv")
    ap.add_argument("--import-ris", type=Path, action="append", default=[],
                    help="a Scopus or Web of Science RIS export; repeatable")
    ap.add_argument("--no-openalex", action="store_true", help="skip citation chasing and the preprint search")
    ap.add_argument("--openalex-key", default=os.environ.get("OPENALEX_API_KEY", ""),
                    help="free OpenAlex API key (default: the OPENALEX_API_KEY environment variable)")
    ap.add_argument("--from-year", type=int, default=2010)
    ap.add_argument("--try-terms", default="",
                    help="candidate terms for the block that loses a seed, separated by ';'. Reports the extra "
                         "PubMed records and the seeds each recovers, then stops (OpenAlex is not called)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.import_ris = [path.expanduser() for path in args.import_ris]
    missing = [str(path) for path in args.import_ris if not path.is_file()]
    if missing:
        sys.exit("RIS file not found: " + ", ".join(missing)
                 + "\nCheck the export names with: ls -lt ~/Downloads | head")
    global OPENALEX_KEY
    OPENALEX_KEY = args.openalex_key.strip()
    args.out.mkdir(parents=True, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc)

    query = " ".join(args.query_file.read_text(encoding="utf-8").split())
    count, pmids, translation = pubmed_search(query, args.email)
    print(f"PubMed: {count} hits")
    if count > len(pmids):
        print(f"  warning: only the first {len(pmids)} PMIDs were returned; narrow the query or page it")
    pubmed = pubmed_fetch(pmids, args.email)
    write_csv(args.out / "pubmed.csv", pubmed)

    with open(args.seeds, newline="", encoding="utf-8") as fh:
        seeds = [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(fh)]
    seed_rows = check_seeds(seeds, pubmed, args.email, query)
    for s in seed_rows:
        note = "" if s["retrieved"] == "yes" else ("  <- query misses it" if s["in_pubmed"] == "yes"
                                                   else "  <- not in PubMed; must come from Scopus, WoS or citations")
        print(f"  seed {s['label']}: retrieved {s['retrieved']}{note}")
        if s.get("missed_by"):
            for block in s["missed_by"].split(" | "):
                print(f"      fails this block: {block[:110]}{'...' if len(block) > 110 else ''}")

    if args.try_terms:
        terms = [t.strip() for t in args.try_terms.split(";") if t.strip()]
        index, trials = try_terms(query, count, seed_rows, terms, args.email)
        if index < 0:
            print("No seed is lost by a query block, so there is nothing to refine.")
            return
        print(f"\nWidening block {index + 1} (now {count} records):")
        print(f"  {'term':58s} {'records':>8s} {'extra':>6s}  recovers")
        for t in trials:
            print(f"  {t['term'][:58]:58s} {t['records']:8d} {t['extra']:+6d}  {t['recovers'] or '-'}")
        md = [f"# Query refinement {started:%Y-%m-%d %H:%M} UTC", "",
              f"Block {index + 1} loses a seed; the query returns {count} records.", "",
              "| Term added to the block | Records | Extra | Seeds recovered |", "| --- | --- | --- | --- |"]
        md += [f"| `{t['term']}` | {t['records']} | {t['extra']:+d} | {t['recovers'] or '-'} |" for t in trials]
        (args.out / "query_refinement.md").write_text("\n".join(md) + "\n", encoding="utf-8")
        print(f"Wrote {args.out / 'query_refinement.md'}")
        return

    backward = forward = preprints = seed_records = []
    missing_seeds: list[str] = []
    problems: list[str] = []
    if not args.no_openalex:
        if not OPENALEX_KEY:
            print("OpenAlex: no API key given; running on the small keyless allowance")
        try:
            backward, forward, missing_seeds, seed_records = chase_citations(seeds, args.email)
            print(f"OpenAlex: {len(backward)} backward and {len(forward)} forward citations of the seeds")
            if missing_seeds:
                print(f"  seeds not found in OpenAlex: {', '.join(missing_seeds)}")
            write_csv(args.out / "citations.csv", backward + forward)
        except (urllib.error.HTTPError, urllib.error.URLError) as err:
            problems.append(f"citation chasing failed ({err}); rerun with an OpenAlex key")
            print(f"  warning: {problems[-1]}")
        try:
            preprints = preprint_search(args.email, args.from_year)
            print(f"OpenAlex preprint search: {len(preprints)} records")
            write_csv(args.out / "preprints.csv", preprints)
        except (urllib.error.HTTPError, urllib.error.URLError) as err:
            problems.append(f"preprint search failed ({err}); rerun with an OpenAlex key")
            print(f"  warning: {problems[-1]}")

    imports = []
    for path in args.import_ris:
        records = read_ris(path, f"import: {path.name}")
        print(f"Imported {len(records)} records from {path.name}")
        imports.append((path.name, records))

    combined, old = combine([seed_records, pubmed, backward, forward, preprints] + [r for _, r in imports], args.from_year)
    write_csv(args.out / "combined.csv", combined)
    write_ris(args.out / "combined.ris", combined)
    print(f"Combined: {len(combined)} unique records ({old} dated before {args.from_year} left out)")

    lines = [f"# Search run {started:%Y-%m-%d %H:%M} UTC", "",
             "## PubMed", "", f"- Query: `{query}`", f"- PubMed's translation: `{translation}`",
             f"- Hits: {count}", "", "## Seed check", "",
             "| Seed | Retrieved by the query | In PubMed | Query blocks that miss it |", "| --- | --- | --- | --- |"]
    lines += [f"| {s['label']} | {s['retrieved']} | {s['in_pubmed']} | {s.get('missed_by', '')[:200].replace('|', '/')} |"
              for s in seed_rows]
    lines += ["", "## Other sources", ""]
    if not args.no_openalex:
        lines += [f"- OpenAlex backward citations of the seeds: {len(backward)}",
                  f"- OpenAlex forward citations of the seeds: {len(forward)}",
                  f"- OpenAlex preprint search (\"{PREPRINT_TERMS}\", type preprint, from {args.from_year}): {len(preprints)}"]
        if missing_seeds:
            lines.append(f"- Seeds not found in OpenAlex: {', '.join(missing_seeds)}")
    lines += [f"- Imported from {name}: {len(records)}" for name, records in imports]
    lines += [f"- Problem: {p}" for p in problems]
    lines += ["", "## Result", "", f"- Unique records after deduplication: {len(combined)}",
              f"- Left out as dated before {args.from_year}: {old}",
              "- Screening file: combined.ris (import into Rayyan or Zotero)"]
    (args.out / "run_log.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.out}/: pubmed.csv, citations.csv, preprints.csv, combined.csv, combined.ris, run_log.md")


if __name__ == "__main__":
    main()
