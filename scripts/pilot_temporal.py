#!/usr/bin/env python3
"""Pilot temporal test: how often does library search confidently misname a new compound?

Searches SWGDRUG entries against the library as it stood before they were added.
Each search is scored twice: with the Stein and Scott (1994) composite match
factor (0 to 999; mass^3, intensity^0.6 weighting; close to, but not identical
with, NIST MS Search) and with spectral entropy similarity (Li et al. 2021, 0 to
1). Entropy similarity was designed for MS/MS; here it is an independently built
second score, to check that findings do not hinge on one formula.

A query whose compound is not in the reference (no shared CAS number or name) is
"novel": any name the search returns is wrong. A "known" query has its compound
in the reference, so the question is whether search finds it or calls it new.

Modes:
  default       entries with ID >= --cutoff-id are queries, lower IDs the reference
  --changelog   every release in the SWGDRUG changelog is evaluated in turn
  --replicates  also searches each spectrum that has a counterpart from a
                different contributor against the library without it (known side)
  --plot        saves the match-factor distributions as a PNG (needs matplotlib)

Extra analyses, printed and saved in --out:
  - 95% Wilson intervals on the headline rates
  - score comparison: known-compound top-1 accuracy, AUROC of the top score for
    known versus new, and a conformal novelty flag calibrated on known spectra
  - isomer groups: how often a new compound's isomer is already in the library,
    and the effect of answering with the isomer group when a different compound
    with the top hit's formula scores within DELTA of it (isomer_rule.csv)
  - provisional compound classes guessed from names (classes.csv)

Examples:
  python3 pilot_temporal.py --lib "jcamp/SWGDRUG 314.HPJ" --compound-list compound_list.xlsx --out pilot
  python3 pilot_temporal.py --lib "jcamp/SWGDRUG 314.HPJ" --changelog Changelog_01222025.xls --replicates --plot --out series
"""
from __future__ import annotations

import argparse
import bisect
import csv
import datetime as dt
import heapq
import math
import re
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
CAS = re.compile(r"\b(\d{1,10})-(\d{2})-(\d)\b")
ENTRY_ID = re.compile(r"entry\s+(\d+)", re.IGNORECASE)
SALT = re.compile(
    r"\(\s*(?:mono|di|tri)?\s*(?:hydrochloride|hydrobromide|hydroiodide|hcl|hbr|citrate|"
    r"fumarate|oxalate|tartrate|maleate|succinate|sulfate|phosphate|salt)\s*\)"
    r"|\b(?:hydrochloride|hydrobromide|hcl|hbr)\b",
    re.IGNORECASE,
)
GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta"}
PRIMES = re.compile(r"['′’`]")
ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")
THRESHOLDS = (600, 650, 700, 750, 800, 850, 900, 950)
ALPHAS = (0.10, 0.20)
DELTAS = (None, 25, 50)
LN2, LN4 = math.log(2), math.log(4)

# Same contributor patterns as count_replicates.py (SWGDRUG FAQ plus sources seen in v3.14).
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

# Provisional classes guessed from names only, first match wins. Leonardo's labels replace these.
CLASSES = [
    ("fentanyl", r"fentan[iy]l"),
    ("nitazene", r"nitazen|etazen"),
    ("synthetic cannabinoid", r"inaca|pica\b|jwh|naphthoyl|ur-144|xlr|cumyl|mdmb|\badb\b|indole-3-carbox|indazole-3-carbox"),
    ("phytocannabinoid", r"\bthc|thc\b|cannabi|\bhhc|\bcbd|\bcbn"),
    ("cathinone", r"cathinone|pyrrolidino\w*phenone|pvp\b|php\b|ylone\b|mephedrone|mdpv|pentedrone|ephedrone|\bbk-|"
                  r"propiophenone|valerophenone|hexanophenone|butiophenone|butanophenone"),
    ("tryptamine", r"tryptamine|dmt\b|dipt\b|mipt\b|nipt\b|psilocin|psilocybin|\bdet\b|\bamt\b"),
    ("arylcyclohexylamine", r"pcp\b|pce\b|ketamine|cyclidine|methoxetamine|fluorexetamine|ephenidine|diphenidine|\bmxe\b|2-oxo|pcpy"),
    ("benzodiazepine", r"azepam|azolam|zolam\b|benzodiazep"),
    ("phenethylamine", r"amphetamine|phenethylamine|nbome|\b2c-|mescaline|\bmda\b|\bmdma\b|benzofuran|\bapb\b|apdb"),
    ("steroid", r"androst|testosterone|nandrolone|estradiol|estrone|estriol|boldenone|trenbolone|stanozolol"),
]


def read_lines(path: Path) -> list[str]:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = data.decode("cp1252", errors="replace")
    return text.splitlines()


def parse_jcamp(lines: list[str]) -> list[dict]:
    """One dict per JCAMP block; wrapped text values are joined back on."""
    records: list[dict] = []
    fields: dict[str, str] = {}
    peaks: list[tuple[float, float]] = []
    in_data, current = False, ""

    def flush() -> None:
        nonlocal fields, peaks, in_data, current
        clean = {k: " ".join(v.split()) for k, v in fields.items()}
        name = clean.get("CASNAME") or clean.get("TITLE", "")
        if peaks and name:
            match = ENTRY_ID.search(clean.get("TITLE", ""))
            records.append({
                "id": int(match.group(1)) if match else None,
                "name": name,
                "formula": clean.get("MOLFORM", ""),
                "cas": clean.get("CASREGISTRYNO", ""),
                "source": clean.get("SAMPLEDESCRIPTION", ""),
                "peaks": peaks,
            })
        fields, peaks, in_data, current = {}, [], False, ""

    for raw in lines:
        body = raw.rstrip("\r\n").split("$$", 1)[0]
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
            in_data = label in ("XYDATA", "PEAKTABLE") and "XY..XY" in value.replace(" ", "").upper()
            fields[label] = value.lstrip()
            current = "" if in_data else label
            continue
        if in_data:
            values = [float(x) for x in NUMBER.findall(body)]
            peaks.extend(zip(values[0::2], values[1::2]))
        elif current:
            fields[current] += body
    flush()
    return records


def norm_name(name: str) -> str:
    text = name.lower()
    for letter, word in GREEK.items():
        text = text.replace(letter, word)
    text = PRIMES.sub("prime", text)  # 3- and 3'- are different positions
    return re.sub(r"[^a-z0-9]", "", SALT.sub("", text))


def identifiers(rec: dict) -> set[tuple[str, str]]:
    ids = set()
    match = CAS.search(rec["cas"])
    if match and int(match.group(1)):
        ids.add(("cas", f"{int(match.group(1))}-{match.group(2)}-{match.group(3)}"))
    name = norm_name(rec["name"])
    if name:
        ids.add(("name", f"{name}|{re.sub(r'[^A-Za-z0-9]', '', rec['formula']).upper()}"))
    return ids


def contributor(source: str) -> str:
    text = source.lower()
    for label, pattern in CONTRIBUTORS:
        if re.search(pattern, text):
            return label
    return "unknown"


def class_guess(name: str) -> str:
    text = name.lower()
    for label, pattern in CLASSES:
        if re.search(pattern, text):
            return label
    return "other"


def formula_key(formula: str) -> str:
    """Element counts in a fixed order, so equal keys mean the same molecular formula."""
    counts: dict[str, int] = defaultdict(int)
    for element, n in ELEMENT.findall(re.sub(r"\s", "", formula)):
        counts[element] += int(n) if n else 1
    return "".join(f"{el}{counts[el]}" for el in sorted(counts))


def unit_spectrum(peaks: list[tuple[float, float]]) -> dict[int, float]:
    """Unit-mass spectrum with the base peak scaled to 1."""
    spec: dict[int, float] = {}
    for mz, intensity in peaks:
        if intensity > 0:
            mass = int(round(mz))
            spec[mass] = max(spec.get(mass, 0.0), intensity)
    top = max(spec.values(), default=0.0)
    return {m: v / top for m, v in spec.items()} if top else {}


def sqrt_cosine(a: dict[int, float], b: dict[int, float]) -> float:
    num = sum((a[m] * b[m]) ** 0.5 for m in a.keys() & b.keys())
    den = (sum(a.values()) * sum(b.values())) ** 0.5
    return num / den if den else 0.0


def prepare_mf(unit: dict[int, float]) -> dict:
    weights = {m: i ** 0.6 * m ** 3 for m, i in unit.items()}
    return {"w": weights, "n2": sum(v * v for v in weights.values()), "n": len(weights)}


def match_factor(q: dict, r: dict) -> float:
    """Stein and Scott composite: weighted dot product plus ratio of adjacent peak pairs."""
    wq, wr = q["w"], r["w"]
    common = sorted(wq.keys() & wr.keys())
    if not common or not q["n2"] or not r["n2"]:
        return 0.0
    dot = sum(wq[m] * wr[m] for m in common)
    fd = dot * dot / (q["n2"] * r["n2"])
    fr = 0.0
    if len(common) > 1:
        total = 0.0
        for a, b in zip(common, common[1:]):
            ratio = (wr[b] / wr[a]) * (wq[a] / wq[b])
            total += ratio if ratio <= 1 else 1 / ratio
        fr = total / (len(common) - 1)
    return 999 * (q["n"] * fd + len(common) * fr) / (q["n"] + len(common))


def prepare_entropy(unit: dict[int, float]) -> dict[int, float]:
    """Intensities as probabilities, with Li et al.'s weighting for low-entropy spectra."""
    total = sum(unit.values())
    if not total:
        return {}
    p = {m: v / total for m, v in unit.items()}
    entropy = -sum(v * math.log(v) for v in p.values() if v > 0)
    if entropy < 3:
        w = 0.25 + 0.25 * entropy
        p = {m: v ** w for m, v in p.items()}
        s = sum(p.values())
        p = {m: v / s for m, v in p.items()}
    return p


def entropy_similarity(a: dict[int, float], b: dict[int, float]) -> float:
    """1 - (2 H(merged) - H(a) - H(b)) / ln 4, computed over shared peaks only."""
    common = a.keys() & b.keys()
    if not common:
        return 0.0
    shared = gain = 0.0
    for m in common:
        x, y = a[m], b[m]
        s = x + y
        gain += x * math.log(x) + y * math.log(y) - s * math.log(s / 2)
        shared += s
    return max(0.0, 1 - (LN2 * (2 - shared) + gain) / LN4)


def build_index(pool: list[dict]) -> dict[tuple[str, str], set[int]]:
    index: dict[tuple[str, str], set[int]] = defaultdict(set)
    for i, rec in enumerate(pool):
        for ident in rec["ids"]:
            index[ident].add(i)
    return index


def evaluate(q: dict, pool: list[dict], targets: set[int], skip: set[int] = frozenset()) -> dict:
    """Search q against pool (minus skip) with both scores; targets hold q's compound."""
    mf_scores, ent_scores = [], []
    for i, r in enumerate(pool):
        if i in skip:
            continue
        mf_scores.append((match_factor(q["mf"], r["mf"]), i))
        ent_scores.append((entropy_similarity(q["ent"], r["ent"]), i))
    top = heapq.nlargest(2, mf_scores)
    mf1, i1 = top[0]
    e1, j1 = max(ent_scores)
    hit, ehit = pool[i1], pool[j1]
    # Best-scoring different compound with the top hit's formula: the isomer the answer could be confused with.
    rivals = [(mf, i) for mf, i in mf_scores if i != i1 and hit["fkey"]
              and pool[i]["fkey"] == hit["fkey"] and not (pool[i]["ids"] & hit["ids"])]
    rival = max(rivals) if rivals else None
    same = lambda a, b: (a["fkey"] == b["fkey"]) if a["fkey"] and b["fkey"] else ""
    row = {
        "release": "", "date": "", "entry_id": q["id"], "name": q["name"], "class_guess": q["class"],
        "contributor": q["contributor"], "status": "known" if targets else "novel",
        "top_hit_id": hit["id"], "top_hit_name": hit["name"], "top_mf": round(mf1),
        "formula": q["formula"], "top_hit_formula": hit["formula"], "same_formula": same(q, hit),
        "isomer_in_reference": any(pool[i]["fkey"] == q["fkey"] for _, i in mf_scores) if q["fkey"] else "",
        "isomer_rival": pool[rival[1]]["name"] if rival else "",
        "isomer_rival_mf": round(rival[0]) if rival else None,
        "hit2_name": pool[top[1][1]]["name"] if len(top) > 1 else "",
        "hit2_mf": round(top[1][0]) if len(top) > 1 else "",
        "entropy_top_name": ehit["name"], "entropy_top_score": round(e1, 4), "entropy_same_formula": same(q, ehit),
        "top_hit_correct": "", "entropy_top_correct": "", "true_match_mf": "", "true_match_rank": "",
        "is_nps": "",
    }
    if targets:
        best = max(mf for mf, i in mf_scores if i in targets)
        row["top_hit_correct"] = i1 in targets
        row["entropy_top_correct"] = j1 in targets
        row["true_match_mf"] = round(best)
        row["true_match_rank"] = 1 + sum(mf > best for mf, _ in mf_scores)
    return row


def run_queries(queries: list[dict], ref: list[dict], label: str, date) -> list[dict]:
    index = build_index(ref)
    rows = []
    for q in queries:
        targets = set().union(*(index.get(ident, set()) for ident in q["ids"]))
        row = evaluate(q, ref, targets)
        row["release"], row["date"] = label, date or ""
        rows.append(row)
    return rows


def replicate_rows(records: list[dict], min_cos: float = 0.8) -> tuple[list[dict], int]:
    """Leave-one-out search for spectra with a counterpart from a different contributor.

    Same-compound spectra from the same contributor leave the pool too (possible copies),
    and counterparts whose spectra disagree (cosine < min_cos) are not counted.
    """
    index = build_index(records)
    rows, dropped = [], 0
    for a, q in enumerate(records):
        if q["contributor"] == "unknown":
            continue
        same = {j for ident in q["ids"] for j in index[ident] if j != a}
        targets = {j for j in same if records[j]["contributor"] not in ("unknown", q["contributor"])}
        dissimilar = {j for j in targets if sqrt_cosine(q["unit"], records[j]["unit"]) < min_cos}
        dropped += bool(targets) and targets == dissimilar
        targets -= dissimilar
        if targets:
            row = evaluate(q, records, targets, skip={a} | (same - targets))
            row["release"], row["status"] = "replicate pair", "known"
            rows.append(row)
    return rows, dropped


def release_ranges(path: Path) -> list[dict]:
    """Entry IDs added in each changelog sheet; each ID belongs to the first release listing it."""
    try:
        import pandas as pd
    except ImportError:
        sys.exit("--changelog needs pandas (and xlrd for .xls files).")
    releases = []
    for name, df in pd.read_excel(path, sheet_name=None).items():
        id_col = next((c for c in df.columns if str(c).strip().upper().startswith("ID")), None)
        if id_col is None:
            continue
        text_col = df.columns[-1]
        dates = pd.to_datetime(df.iloc[:, 0], errors="coerce").dropna()
        ids: set[int] = set()
        for _, row in df.iterrows():
            if not re.search(r"add\w*\s+(new\s+)?entr", str(row[text_col]).lower()):
                continue
            for token in re.split(r"[,;&]|\band\b", str(row[id_col])):
                token = re.sub(r"\.0$", "", token.strip())
                span = re.fullmatch(r"(\d+)\s*-\s*(\d+)", token)
                if span:
                    ids.update(range(int(span.group(1)), int(span.group(2)) + 1))
                elif token.isdigit():
                    ids.add(int(token))
        if ids:
            releases.append({"release": name.strip(),
                             "date": dates.min().date() if len(dates) else None, "ids": ids})
    releases.sort(key=lambda r: (r["date"] or dt.date.max, min(r["ids"])))
    seen: set[int] = set()
    for rel in releases:
        rel["ids"] -= seen
        seen |= rel["ids"]
    return [r for r in releases if r["ids"]]


def red_names(path: Path) -> set[str] | None:
    """Names the SWGDRUG compound list marks in red (new in the latest version)."""
    try:
        from openpyxl import load_workbook
    except ImportError:
        print("openpyxl is not installed, so the compound-list check is skipped.")
        return None
    names = set()
    for cell in load_workbook(path).active["A"][1:]:
        rgb = getattr(cell.font.color, "rgb", None) if cell.font and cell.font.color else None
        if not (cell.value and isinstance(rgb, str) and len(rgb) >= 6):
            continue
        try:
            r, g, b = (int(rgb[-6:][i:i + 2], 16) for i in (0, 2, 4))
        except ValueError:
            continue
        if r > 150 and g < 100 and b < 100:
            names.add(" ".join(str(cell.value).split()).casefold())
    return names


def pct(n: int, d: int) -> str:
    return f"{100 * n / d:.0f}%" if d else "n/a"


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if not n:
        return 0.0, 0.0
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def auroc(pos: list[float], neg: list[float]) -> float:
    """Probability that a known compound's top score beats a new compound's (ties count half)."""
    if not pos or not neg:
        return float("nan")
    neg = sorted(neg)
    wins = sum(bisect.bisect_left(neg, p) + 0.5 * (bisect.bisect_right(neg, p) - bisect.bisect_left(neg, p))
               for p in pos)
    return wins / (len(pos) * len(neg))


def conformal_p(calibration: list[float], score: float) -> float:
    """Conformal p-value that a query is a known compound; low = unusually poor best match."""
    return (1 + sum(c <= score for c in calibration)) / (len(calibration) + 1)


def known_outcomes(rows: list[dict], t: int) -> tuple[int, int]:
    """(called new, given a wrong name) among known queries at threshold t."""
    flagged = sum(r["top_mf"] < t for r in rows)
    wrong = sum(r["top_mf"] >= t and r["top_hit_correct"] is not True for r in rows)
    return flagged, wrong


def outcome(r: dict, t: int, delta) -> str:
    """Name the top hit, answer with its isomer group, or call the query new."""
    if r["top_mf"] < t:
        decision = "new"
    elif delta is not None and r["isomer_rival_mf"] is not None and r["isomer_rival_mf"] >= r["top_mf"] - delta:
        decision = "group"
    else:
        decision = "name"
    if decision == "group":
        return "isomer group" if r["same_formula"] is True else "wrong"
    if r["status"] == "novel":
        return "new" if decision == "new" else "wrong"
    if decision == "new":
        return "called new"
    return "correct name" if r["top_hit_correct"] is True else "wrong"


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot(novel: list[dict], known: list[dict], path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed, so the plot is skipped.")
        return
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    bins = list(range(0, 1001, 25))
    if novel:
        ax.hist([r["top_mf"] for r in novel], bins=bins, density=True, alpha=0.55,
                label=f"New compound: best (wrong) hit, n={len(novel)}")
    if known:
        ax.hist([r["true_match_mf"] for r in known], bins=bins, density=True, alpha=0.55,
                label=f"Known compound: true match, n={len(known)}")
    for t in (800, 900):
        ax.axvline(t, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel("Match factor (Stein and Scott composite, 0 to 999)")
    ax.set_ylabel("Density")
    ax.set_title("Best wrong hit for new compounds vs true match for known compounds")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    print(f"Saved {path}")


def score_comparison(novel: list[dict], known: list[dict], out: Path) -> None:
    k = sum(r["top_mf"] >= 800 for r in novel)
    print("\nScore comparison (scale-free, so the two scores can be set side by side):")
    print(f"{'score':14s} {'known top-1':>11s} {'AUROC':>6s} "
          + " ".join(f"{'flags new @' + format(a, '.0%'):>14s}" for a in ALPHAS)
          + f" {'isomer share, top ' + str(k):>22s}")
    rows = []
    for label, key, correct, same in (("match factor", "top_mf", "top_hit_correct", "same_formula"),
                                      ("entropy", "entropy_top_score", "entropy_top_correct", "entropy_same_formula")):
        cal = [r[key] for r in known]
        right = sum(r[correct] is True for r in known)
        auc = auroc(cal, [r[key] for r in novel])
        flags = [sum(conformal_p(cal, r[key]) <= a for r in novel) for a in ALPHAS]
        loo = [sum(conformal_p(cal[:i] + cal[i + 1:], s) <= a for i, s in enumerate(cal)) for a in ALPHAS]
        confident = sorted(novel, key=lambda r: -r[key])[:k]
        iso = sum(r[same] is True for r in confident)
        print(f"{label:14s} {pct(right, len(known)):>11s} {auc:6.3f} "
              + " ".join(f"{pct(f, len(novel)):>14s}" for f in flags) + f" {pct(iso, k):>22s}")
        row = {"score": label, "known_n": len(known), "known_top1": right, "auroc": round(auc, 4),
               "novel_n": len(novel), "isomer_share_top_k": round(100 * iso / k, 1) if k else "", "k": k}
        for a, f, l in zip(ALPHAS, flags, loo):
            row[f"flags_new_alpha_{a:.2f}"] = f
            row[f"known_flagged_loo_alpha_{a:.2f}"] = l
        rows.append(row)
    print("  'flags new @ a' = conformal flag calibrated on the known spectra at false-novelty rate a;")
    print("  with so few known spectra (each pair counted in both directions) this is a dry run.")
    for r in rows:
        loo = ", ".join(f"{r[f'known_flagged_loo_alpha_{a:.2f}']}/{r['known_n']} at {a:.0%}" for a in ALPHAS)
        print(f"  leave-one-out check, known spectra flagged as new ({r['score']}): {loo}")
    write_csv(out / "score_comparison.csv", rows)


def isomer_analysis(novel: list[dict], known: list[dict], out: Path) -> None:
    have = sum(r["isomer_in_reference"] is True for r in novel)
    print(f"\nIsomer groups: {have} of {len(novel)} new compounds ({pct(have, len(novel))}) "
          "already had an isomer (same formula) in the library.")
    print("Rule: answer with the isomer group when a different compound with the top hit's formula scores within delta.")
    print(f"{'threshold':>9s} {'delta':>5s} | {'new: wrong':>10s} {'new: group':>10s} | "
          f"{'known: right':>12s} {'known: group':>12s} {'known: wrong':>12s} {'known: new':>10s}")
    rows = []
    for t in (800, 900):
        for delta in DELTAS:
            n = Counter(outcome(r, t, delta) for r in novel)
            kn = Counter(outcome(r, t, delta) for r in known)
            rows.append({"threshold": t, "delta": "off" if delta is None else delta,
                         "novel_n": len(novel), "novel_wrong": n["wrong"], "novel_group": n["isomer group"],
                         "known_n": len(known), "known_correct_name": kn["correct name"],
                         "known_group": kn["isomer group"], "known_wrong": kn["wrong"],
                         "known_called_new": kn["called new"]})
            print(f"{t:9d} {('off' if delta is None else delta):>5} | {pct(n['wrong'], len(novel)):>10s} "
                  f"{pct(n['isomer group'], len(novel)):>10s} | {pct(kn['correct name'], len(known)):>12s} "
                  f"{pct(kn['isomer group'], len(known)):>12s} {pct(kn['wrong'], len(known)):>12s} "
                  f"{pct(kn['called new'], len(known)):>10s}")
    write_csv(out / "isomer_rule.csv", rows)


def class_table(novel: list[dict], out: Path) -> None:
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in novel:
        groups[r["class_guess"]].append(r)
    print("\nProvisional classes (guessed from names; Leonardo's labels replace these):")
    print(f"{'class':22s} {'new':>5s} {'wrong >=800':>12s} {'of which isomers':>17s}")
    rows = []
    for cls, rs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        wrong = [r for r in rs if r["top_mf"] >= 800]
        iso = sum(r["same_formula"] is True for r in wrong)
        rows.append({"class_guess": cls, "new_compounds": len(rs), "wrong_800": len(wrong),
                     "wrong_800_pct": round(100 * len(wrong) / len(rs), 1), "isomers": iso})
        print(f"{cls:22s} {len(rs):5d} {len(wrong):5d} ({pct(len(wrong), len(rs)):>4s}) {iso:9d} ({pct(iso, len(wrong)):>4s})")
    write_csv(out / "classes.csv", rows)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lib", required=True, type=Path, help="SWGDRUG JCAMP-DX file (.HPJ)")
    ap.add_argument("--cutoff-id", type=int, default=5814,
                    help="single-cutoff mode: first entry ID treated as a later arrival (default 5814)")
    ap.add_argument("--changelog", type=Path, help="SWGDRUG changelog (.xls): evaluate every release")
    ap.add_argument("--replicates", action="store_true", help="also run the known-compound check")
    ap.add_argument("--plot", action="store_true", help="save mf_distributions.png")
    ap.add_argument("--compound-list", type=Path,
                    help="optional SWGDRUG compound list .xlsx, to check IDs against its red names")
    ap.add_argument("--out", type=Path, default=Path("pilot_out"))
    args = ap.parse_args(argv)

    path = args.lib.expanduser()
    if not path.is_file():
        sys.exit(f"File not found: {path}")
    args.out.mkdir(parents=True, exist_ok=True)
    records = sorted((r for r in parse_jcamp(read_lines(path)) if r["id"] is not None), key=lambda r: r["id"])
    if not records:
        sys.exit("No entries with 'Library Entry <ID>' titles were found.")
    for rec in records:
        rec["unit"] = unit_spectrum(rec["peaks"])
        rec["mf"], rec["ent"] = prepare_mf(rec["unit"]), prepare_entropy(rec["unit"])
        rec["ids"], rec["fkey"] = identifiers(rec), formula_key(rec["formula"])
        rec["contributor"], rec["class"] = contributor(rec["source"]), class_guess(rec["name"])
    print(f"Library: {len(records)} entries, IDs {records[0]['id']} to {records[-1]['id']}")

    if args.compound_list:
        red = red_names(args.compound_list.expanduser())
        if red is not None:
            key = lambda r: " ".join(r["name"].split()).casefold()
            late = sum(key(r) in red for r in records if r["id"] >= args.cutoff_id)
            early = sum(key(r) in red for r in records if r["id"] < args.cutoff_id)
            print(f"Compound list marks {len(red)} names as new; {late} have IDs >= {args.cutoff_id} "
                  f"and {early} lower (the ID split is right if the second number is near 0).")

    started = time.time()
    rows: list[dict] = []
    release_summary: list[dict] = []
    if args.changelog:
        releases = release_ranges(args.changelog.expanduser())
        print(f"Changelog: {len(releases)} releases list added entry IDs\n")
        print(f"{'release':12s} {'date':10s} {'queries':>7s} {'known':>5s} "
              f"{'wrong >=800':>12s} {'wrong >=900':>12s}")
        for rel in releases:
            first = min(rel["ids"])
            ref = [r for r in records if r["id"] < first]
            queries = [r for r in records if r["id"] in rel["ids"]]
            if not ref or not queries:
                continue
            got = run_queries(queries, ref, rel["release"], rel["date"])
            rows += got
            novel = [r for r in got if r["status"] == "novel"]
            w8 = sum(r["top_mf"] >= 800 for r in novel)
            w9 = sum(r["top_mf"] >= 900 for r in novel)
            release_summary.append({"release": rel["release"], "date": rel["date"], "queries": len(got),
                                    "known": len(got) - len(novel), "novel": len(novel),
                                    "wrong_800": w8, "wrong_900": w9})
            print(f"{rel['release']:12s} {str(rel['date']):10s} {len(got):7d} {len(got) - len(novel):5d} "
                  f"{w8:5d} ({pct(w8, len(novel)):>4s}) {w9:5d} ({pct(w9, len(novel)):>4s})"
                  f"   [{time.time() - started:.0f} s]")
    else:
        ref = [r for r in records if r["id"] < args.cutoff_id]
        queries = [r for r in records if r["id"] >= args.cutoff_id]
        print(f"Reference: {len(ref)} entries with ID < {args.cutoff_id}; queries: {len(queries)}")
        if not ref or not queries:
            sys.exit("One side of the split is empty; check --cutoff-id.")
        rows = run_queries(queries, ref, f"ID >= {args.cutoff_id}", None)

    novel = [r for r in rows if r["status"] == "novel"]
    known = [r for r in rows if r["status"] == "known"]
    print(f"\nNovel queries: {len(novel)}")
    for t in (700, 800, 900):
        hits = [r for r in novel if r["top_mf"] >= t]
        lo, hi = wilson(len(hits), len(novel))
        iso = sum(r["same_formula"] is True for r in hits)
        print(f"  wrong name with match factor >= {t}: {len(hits)} ({pct(len(hits), len(novel))}, "
              f"95% CI {100 * lo:.0f} to {100 * hi:.0f}%); top hit is an isomer (same formula): "
              f"{iso} ({pct(iso, len(hits))})")

    replicates: list[dict] = []
    if args.replicates:
        replicates, dropped = replicate_rows(records)
        print(f"\nKnown-compound check: {len(replicates)} spectra with a counterpart from another contributor"
              f" ({dropped} left out because the counterpart's spectrum disagrees, cosine < 0.8)")
    if known:
        print(f"Compounds re-added in a later release (known, but possibly re-deposited copies): {len(known)}")
    # The known side uses independent replicates only; re-added entries can be copies
    # from the same contributor, which would make known compounds look easier than they are.
    known_all = replicates if args.replicates else known
    if known_all:
        right = sum(r["top_hit_correct"] is True for r in known_all)
        lo, hi = wilson(right, len(known_all))
        mfs = [r["true_match_mf"] for r in known_all]
        source = "independent replicates" if args.replicates else "re-added entries (may include copies)"
        print(f"Known side ({source}): {len(known_all)} queries; true match ranked first: {right} "
              f"({pct(right, len(known_all))}, 95% CI {100 * lo:.0f} to {100 * hi:.0f}%); "
              f"median true-match score {statistics.median(mfs):.0f}")
        missed = [r for r in known_all if r["top_hit_correct"] is not True]
        iso = sum(r["same_formula"] is True for r in missed)
        print(f"  outranked by another compound: {len(missed)}; of those, by an isomer (same formula): {iso}")

        print("\nThreshold trade-off (a lab accepts the top hit only if its match factor reaches the threshold):")
        print(f"{'threshold':>9s} {'new: wrong name':>16s} {'known: called new':>18s} {'known: wrong name':>18s}")
        tradeoff = []
        for t in THRESHOLDS:
            wrong_new = sum(r["top_mf"] >= t for r in novel)
            flagged, wrong_known = known_outcomes(known_all, t)
            tradeoff.append({"threshold": t, "novel_n": len(novel), "novel_wrong_name": wrong_new,
                             "known_n": len(known_all), "known_called_new": flagged,
                             "known_wrong_name": wrong_known})
            print(f"{t:9d} {pct(wrong_new, len(novel)):>16s} {pct(flagged, len(known_all)):>18s} "
                  f"{pct(wrong_known, len(known_all)):>18s}")
        write_csv(args.out / "tradeoff.csv", tradeoff)
        score_comparison(novel, known_all, args.out)
    isomer_analysis(novel, known_all, args.out)
    class_table(novel, args.out)

    print("\nHighest-scoring novel queries (each top hit is a wrong name):")
    for r in sorted(novel, key=lambda r: -r["top_mf"])[:20]:
        print(f"  {r['top_mf']:4d}  {r['name'][:42]:42s} -> {r['top_hit_name'][:42]}  [{r['release']}]")

    write_csv(args.out / "queries.csv", sorted(rows, key=lambda r: (r["status"], -r["top_mf"])))
    write_csv(args.out / "replicates.csv", replicates)
    write_csv(args.out / "releases.csv", release_summary)
    if args.plot:
        plot(novel, known_all, args.out / "mf_distributions.png")
    print(f"\nWrote results to {args.out}/ in {time.time() - started:.0f} s (is_nps column in queries.csv is for Leonardo).")
    print("A 'novel' query whose compound sits in the reference under another name and CAS is miscounted;"
          " check the top of the list before quoting numbers.")


if __name__ == "__main__":
    main()
