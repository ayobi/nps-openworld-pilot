#!/usr/bin/env python3
"""Count independent replicate EI spectra across forensic spectral libraries.

The NPS open-world benchmark needs, for every "known compound" query, a second
spectrum of the same compound that did not come from the same source as the
library entry. This script groups library entries by compound, works out who
contributed each spectrum, and labels every pair of entries as:

  independent         different known contributors, spectra not copies
  likely_copy         different stated contributors but the same spectrum
  same_contributor    both spectra from one contributor
  unknown_provenance  at least one contributor could not be identified
  low_similarity      different contributors but dissimilar spectra (cosine below
                      --min-cos): check by hand, since this is how isomers merged
                      by mistake show up, though a degraded spectrum can too

Input: MSP text or JCAMP-DX text, such as SWGDRUG's .HPJ download. The format
is detected from the file contents. Libraries that only ship in NIST MS Search
or Agilent formats must be exported to MSP first (NIST Lib2NIST, Windows only).

Examples:
  python3 count_replicates.py --lib swgdrug="jcamp/SWGDRUG 314.HPJ" --out replicates
  python3 count_replicates.py --lib swgdrug="jcamp/SWGDRUG 314.HPJ" --lib cayman=CSL.msp --out replicates
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Contributor tokens from the SWGDRUG MS library FAQ (https://www.swgdrug.org/ms.htm).
# Order matters: specific DEA laboratories are checked before the generic DEA token.
CONTRIBUTORS = [
    ("Cayman", r"cayman|\bitem\s*\d{4,6}\s*batch"),
    ("NPS Discovery (CFSRE)", r"nps\s*discovery|cfsre|center for forensic science research"),
    ("KGI (Krstenansky lab)", r"2016-r2-cx-0059|keck\s*graduate|krstenansky"),
    ("ForensicDB (RTI)", r"forensicdb|\brti\b"),
    ("ENFSI-DWG", r"enfsi"),
    ("AAFS", r"\baafs\b"),
    ("QHFSS", r"qhfss"),
    ("Broward Sheriff Office", r"broward"),
    ("US Army (J. Nawyn)", r"nawyn"),
    ("DEA SFL1", r"\bsfl\s*1\b"),
    ("DEA SWL", r"\bswl\b"),
    ("DEA SCL", r"\bscl\b"),
    ("DEA WRL", r"\bwrl\b"),
    ("DEA", r"\bdea\b"),
    ("Agilent (T. Gluodenis)", r"gluodenis"),
    ("HSA Singapore", r"\bhsa\b"),
    ("North Carolina OCME", r"north\s*carolina\s*office\s*of\s*the\s*chief"),
    ("UK Forensic Science Service", r"forensic\s*science\s*service"),
    ("EFSI Estonia", r"\befsi\b"),
    ("BKA Germany", r"bundeskriminalamt|\bbka\b"),
    ("J. Gerlits", r"gerlits"),
]

NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
CAS = re.compile(r"\b(\d{1,10})-(\d{2})-(\d)\b")
SALT = re.compile(
    r"\(\s*(?:mono|di|tri)?\s*(?:hydrochloride|hydrobromide|hydroiodide|hcl|hbr|citrate|"
    r"fumarate|oxalate|tartrate|maleate|succinate|sulfate|phosphate|salt)\s*\)"
    r"|\b(?:hydrochloride|hydrobromide|hcl|hbr)\b",
    re.IGNORECASE,
)
GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta"}
PRIMES = re.compile(r"['′’`]")
STATUSES = ["independent", "low_similarity", "likely_copy", "same_contributor", "unknown_provenance"]
SOURCE_FIELDS = ("comment", "comments", "source", "contributor")


def read_lines(path: Path) -> list[str]:
    """Decode as UTF-8, falling back to Windows-1252 for older library exports."""
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = data.decode("cp1252", errors="replace")
    return text.splitlines()


def is_jcamp(lines: list[str]) -> bool:
    for line in lines:
        if line.strip():
            return line.lstrip().startswith("##")
    return False


def parse_msp(lines: list[str]) -> list[dict]:
    """MSP records as dicts of lower-cased fields plus 'peaks' [(mz, intensity)]."""
    records: list[dict] = []
    rec: dict = {}
    peaks: list[tuple[float, float]] = []
    in_peaks = False

    def flush() -> None:
        nonlocal rec, peaks, in_peaks
        if rec.get("name"):
            rec["peaks"] = peaks
            records.append(rec)
        rec, peaks, in_peaks = {}, [], False

    for raw in lines:
        line = raw.strip()
        if not line:
            flush()
            continue
        starts_record = line.lower().startswith("name:")
        if in_peaks and not starts_record:
            values = [float(x) for x in NUMBER.findall(re.sub(r'"[^"]*"', "", line))]
            peaks.extend(zip(values[0::2], values[1::2]))
            continue
        if starts_record and rec:
            flush()
        key, sep, value = line.partition(":")
        if not sep:
            continue
        key = key.strip().lower()
        if key == "num peaks":
            in_peaks = True
            continue
        value = value.strip()
        rec[key] = f"{rec[key]} | {value}" if key in rec else value
    flush()
    return records


def parse_jcamp(lines: list[str]) -> list[dict]:
    """JCAMP-DX mass spectra, one block per compound, mapped onto the MSP-style dicts.

    Text values may wrap onto following lines (SWGDRUG breaks them at a fixed
    width, sometimes mid-word), so continuation lines are joined back on.
    """
    records: list[dict] = []
    fields: dict[str, str] = {}
    peaks: list[tuple[float, float]] = []
    in_data = False
    current = ""  # label whose value may continue on the next line

    def flush() -> None:
        nonlocal fields, peaks, in_data, current
        clean = {k: " ".join(v.split()) for k, v in fields.items()}
        name = clean.get("CASNAME") or clean.get("TITLE", "")
        if peaks and name:
            records.append({
                "name": name,
                "formula": clean.get("MOLFORM", ""),
                "cas#": clean.get("CASREGISTRYNO", ""),
                "inchikey": clean.get("$INCHIKEY", ""),
                "source": " | ".join(clean[k] for k in ("SAMPLEDESCRIPTION", "ORIGIN", "OWNER")
                                     if clean.get(k)),
                "peaks": peaks,
            })
        fields, peaks, in_data, current = {}, [], False, ""

    for raw in lines:
        body = raw.rstrip("\r\n").split("$$", 1)[0]  # $$ starts a JCAMP comment
        if not body.strip():
            continue
        head = body.lstrip()
        if head.startswith("##"):
            label, _, value = head[2:].partition("=")
            label = re.sub(r"[\s\-_/]", "", label.upper())
            if label == "END":
                flush()
                continue
            if label == "TITLE" and (fields or peaks):
                flush()
            # Only (XY..XY) peak lists are EI spectra; (X++(Y..Y)) is equally spaced data.
            in_data = label in ("XYDATA", "PEAKTABLE") and "XY..XY" in value.replace(" ", "").upper()
            fields[label] = value.lstrip()
            current = "" if in_data else label
            continue
        if in_data:
            values = [float(x) for x in NUMBER.findall(body)]
            peaks.extend(zip(values[0::2], values[1::2]))
        elif current:
            fields[current] += body  # wrapped value: join without adding a space
    flush()
    return records


def norm_name(name: str) -> str:
    text = name.lower()
    for letter, word in GREEK.items():
        text = text.replace(letter, word)
    text = PRIMES.sub("prime", text)  # 3- and 3'- are different positions, so the prime must survive
    return re.sub(r"[^a-z0-9]", "", SALT.sub("", text))


def identifiers(rec: dict) -> list[tuple[str, str]]:
    """InChIKey skeleton, CAS number (leading zeros dropped) and name+formula, where present."""
    ids = []
    match = re.match(r"\s*([A-Z]{14})-[A-Z]{10}-[A-Z]", rec.get("inchikey", ""))
    if match:
        ids.append(("inchikey", match.group(1)))
    for field in ("cas#", "casno", "cas"):
        match = CAS.search(rec.get(field, ""))
        if match:
            head = int(match.group(1))
            if head:  # 0-00-0 and all-zero placeholders mean "no CAS"
                ids.append(("cas", f"{head}-{match.group(2)}-{match.group(3)}"))
            break
    name = norm_name(rec.get("name", ""))
    if name:
        formula = re.sub(r"\s", "", rec.get("formula", "")).upper()
        ids.append(("name", f"{name}|{formula}"))
    return ids


def source_text(rec: dict) -> str:
    return " ".join(rec.get(k, "") for k in SOURCE_FIELDS).strip()


def contributor(rec: dict, lib: str, fixed: dict[str, str]) -> str:
    if lib in fixed:
        return fixed[lib]
    text = source_text(rec).lower()
    for label, pattern in CONTRIBUTORS:
        if re.search(pattern, text):
            return label
    return "unknown"


def normalised(peaks: list[tuple[float, float]]) -> dict[int, float]:
    """Unit-mass spectrum scaled to a base peak of 1."""
    spec: dict[int, float] = {}
    for mz, intensity in peaks:
        if intensity > 0:
            mass = int(round(mz))
            spec[mass] = max(spec.get(mass, 0.0), intensity)
    top = max(spec.values(), default=0.0)
    return {m: v / top for m, v in spec.items()} if top else {}


def cosine(a: dict[int, float], b: dict[int, float]) -> float:
    """Square-root-intensity cosine, reported for information only."""
    num = sum(math.sqrt(a[m] * b[m]) for m in a.keys() & b.keys())
    den = math.sqrt(sum(a.values()) * sum(b.values()))
    return num / den if den else 0.0


def likely_copy(a: dict[int, float], b: dict[int, float], tol: float) -> bool:
    """Same peaks and every relative intensity within tol: a re-deposited spectrum."""
    return bool(a) and a.keys() == b.keys() and all(abs(a[m] - b[m]) <= tol for m in a)


def pair_status(ea: dict, eb: dict, tol: float, cos: float, min_cos: float) -> str:
    if "unknown" in (ea["contributor"], eb["contributor"]):
        return "unknown_provenance"
    if ea["contributor"] == eb["contributor"]:
        return "same_contributor"
    if likely_copy(ea["spectrum"], eb["spectrum"], tol):
        return "likely_copy"
    if cos < min_cos:
        return "low_similarity"
    return "independent"


def matched_by(ea: dict, eb: dict) -> str:
    shared = {kind for kind, _ in set(ea["ids"]) & set(eb["ids"])}
    for kind in ("inchikey", "cas", "name"):
        if kind in shared:
            return kind
    return "transitive"


def group_entries(entries: list[dict]) -> list[list[int]]:
    """Union entries that share any identifier."""
    root = list(range(len(entries)))

    def find(i: int) -> int:
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    first_seen: dict[tuple[str, str], int] = {}
    for i, entry in enumerate(entries):
        for ident in entry["ids"]:
            j = first_seen.setdefault(ident, i)
            ri, rj = find(i), find(j)
            if ri != rj:
                root[ri] = rj
    groups: dict[int, list[int]] = defaultdict(list)
    for i in range(len(entries)):
        groups[find(i)].append(i)
    return sorted(groups.values(), key=lambda g: entries[g[0]]["name"].lower())


def split_pair(text: str, flag: str) -> tuple[str, str]:
    label, sep, value = text.partition("=")
    if not sep or not label.strip() or not value.strip():
        sys.exit(f"{flag} expects LABEL=VALUE, got {text!r}")
    return label.strip(), value.strip()


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lib", action="append", required=True, metavar="LABEL=PATH",
                    help="an MSP or JCAMP-DX library; repeat for each library")
    ap.add_argument("--fixed-contributor", action="append", default=[], metavar="LABEL=NAME",
                    help="treat every entry of LABEL as contributed by NAME "
                         "(labels containing 'cayman' default to Cayman)")
    ap.add_argument("--copy-tol", type=float, default=0.005,
                    help="max per-peak difference, as a fraction of the base peak, for two spectra "
                         "with identical peak lists to count as copies (default 0.005)")
    ap.add_argument("--min-cos", type=float, default=0.8,
                    help="pairs below this square-root cosine are not counted as independent "
                         "replicates (default 0.8)")
    ap.add_argument("--out", type=Path, default=Path("replicates_out"))
    args = ap.parse_args(argv)

    libs = {label: Path(path).expanduser() for label, path in
            (split_pair(s, "--lib") for s in args.lib)}
    missing = [str(p) for p in libs.values() if not p.is_file()]
    if missing:
        sys.exit("File not found: " + ", ".join(missing)
                 + '\nPass the real library files, for example --lib swgdrug="jcamp/SWGDRUG 314.HPJ"')
    fixed = {label: "Cayman" for label in libs if "cayman" in label.lower()}
    fixed.update(dict(split_pair(s, "--fixed-contributor") for s in args.fixed_contributor))

    entries = []
    unknown_sources: Counter = Counter()
    for label, path in libs.items():
        lines = read_lines(path)
        fmt = "JCAMP-DX" if is_jcamp(lines) else "MSP"
        records = parse_jcamp(lines) if fmt == "JCAMP-DX" else parse_msp(lines)
        print(f"{label}: read {len(records)} spectra as {fmt} from {path.name}")
        for rec in records:
            who = contributor(rec, label, fixed)
            if who == "unknown":
                unknown_sources[source_text(rec) or "(no source field)"] += 1
            entries.append({
                "lib": label,
                "name": rec.get("name", ""),
                "ids": identifiers(rec),
                "contributor": who,
                "spectrum": normalised(rec["peaks"]),
            })
    if not entries:
        sys.exit("No spectra found; check that the files are MSP or JCAMP-DX text.")

    groups = group_entries(entries)
    args.out.mkdir(parents=True, exist_ok=True)
    statuses: Counter = Counter()
    evidence: Counter = Counter()
    usable = multi = 0
    with open(args.out / "pairs.csv", "w", newline="", encoding="utf-8") as pf, \
            open(args.out / "compounds.csv", "w", newline="", encoding="utf-8") as cf:
        pairs = csv.writer(pf)
        compounds = csv.writer(cf)
        pairs.writerow(["compound_id", "matched_by", "name_a", "library_a", "contributor_a",
                        "name_b", "library_b", "contributor_b", "cosine", "status"])
        compounds.writerow(["compound_id", "names", "entries", "libraries", "contributors",
                            "independent_pairs"])
        for cid, members in enumerate(groups, start=1):
            multi += len(members) > 1
            independent = 0
            for i, j in itertools.combinations(members, 2):
                ea, eb = entries[i], entries[j]
                cos = cosine(ea["spectrum"], eb["spectrum"])
                status = pair_status(ea, eb, args.copy_tol, cos, args.min_cos)
                how = matched_by(ea, eb)
                statuses[status] += 1
                evidence[how] += 1
                independent += status == "independent"
                pairs.writerow([cid, how, ea["name"], ea["lib"], ea["contributor"], eb["name"],
                                eb["lib"], eb["contributor"],
                                f"{cos:.4f}", status])
            usable += independent > 0
            members_e = [entries[i] for i in members]
            compounds.writerow([cid, "; ".join(sorted({e["name"] for e in members_e})), len(members),
                                "; ".join(sorted({e["lib"] for e in members_e})),
                                "; ".join(sorted({e["contributor"] for e in members_e})), independent])

    per_contrib = Counter(e["contributor"] for e in entries)
    print("Entries per contributor: " + ", ".join(f"{k} {v}" for k, v in per_contrib.most_common()))
    print(f"Compounds: {len(groups)} ({multi} with more than one entry)")
    print(f"Compounds with at least one independent pair: {usable}"
          "   <- known-compound queries available")
    print("Pairs by status:       " + ", ".join(f"{s} {statuses[s]}" for s in STATUSES))
    print("Pairs matched by:      " + ", ".join(f"{k} {evidence[k]}"
                                                 for k in ("inchikey", "cas", "name", "transitive")))
    if statuses["low_similarity"]:
        print(f"Low-similarity pairs (cosine < {args.min_cos}) are not counted as independent; "
              "check them in pairs.csv.")
    if evidence["name"] or evidence["transitive"]:
        print("Name and transitive matches are the weakest evidence; check them in pairs.csv.")
    if unknown_sources:
        print("Unrecognised source text (most common first):")
        for text, n in unknown_sources.most_common(15):
            print(f"  {n:5d}  {text[:80]}")
    print(f"Wrote {args.out / 'pairs.csv'} and {args.out / 'compounds.csv'}")


if __name__ == "__main__":
    main()
