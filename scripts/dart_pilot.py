#!/usr/bin/env python3
"""DART-MS pilot: does the isomer problem carry over from GC-EI-MS to DART-MS?

Reads releases of the NIST DART-MS Forensics Database (.SDF). Each compound has
in-source CID spectra at three energies and an InChIKey, so isomers and
stereoisomers are told apart by structure rather than by name. Energies are read
from the end of each name: "+30 V", "+60 V", "+90 V" in the 2021-2022 releases,
"LOW", "MED", "HIGH" from 2023 on; both are treated as levels low, mid, high.

Modes:
  one release    leave-one-compound-out: every compound is searched against the
                 rest of the library, so its best hit is necessarily another
                 compound (a simulated new compound)
  2+ releases    temporal, oldest first: a compound is new in a release when its
                 connectivity (first InChIKey block) was absent from the release
                 before. New compounds are searched against the compounds that
                 were already there, using the NEWER release's spectra for both,
                 so a change of spectral processing between releases cannot leak
                 into the comparison. Compounds whose InChIKey changed but whose
                 connectivity did not are reported as re-keyed, not new.

Spectra: peaks below --min-rel of the base peak are dropped (the 2023+ releases
carry raw intensities with many small peaks). Similarity: square-root-intensity
cosine per energy, peaks matched within --tol Da, averaged over the energies
both compounds have (0 to 1).

Energy weighting (2+ releases): four ways of combining the per-energy scores are
compared: mean, high-weighted (low:mid:high = 1:2:3), min (must match at every
energy) and high only. Because each changes the scale, they are compared at an
equal false-novelty budget: the threshold is set so that 5% or 10% of KNOWN
compounds would be called new, and the share of new compounds given a wrong name
at that threshold is reported. Known compounds here are compounds present in two
consecutive releases whose spectra changed in between (mostly reprocessing at the
2023 format change): the newer spectrum is searched against the older library.
Compounds carried over with identical spectra are left out. This is a proxy for
replicate measurements, not a substitute.

Examples:
  python3 dart_pilot.py --release dragonfly=He_Pos_Dragonfly.SDF --out dart
  python3 dart_pilot.py --release a=OLD.SDF --release b=NEWER.SDF --release c=NEWEST.SDF --out dart
"""
from __future__ import annotations

import argparse
import csv
import heapq
import math
import re
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

TAG = re.compile(r"^>\s*<([^>]+)>")
SUFFIX = re.compile(r"\s+(?:\+\s*(\d+)\s*V|(low|med|medium|mid|hi|high))\s*$", re.IGNORECASE)
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")
WORD_LEVEL = {"low": 1, "med": 2, "medium": 2, "mid": 2, "hi": 3, "high": 3}
VOLT_LEVEL = {30: 1, 60: 2, 90: 3}
LEVEL_NAME = {1: "low", 2: "mid", 3: "high"}
THRESHOLDS = (0.80, 0.90, 0.95)
LEVEL_WEIGHT = {1: 1.0, 2: 2.0, 3: 3.0}
VARIANTS = (("mean", "mean"), ("high-weighted", "weighted"), ("min", "min"), ("high only", "high"))
BUDGETS = (0.05, 0.10)


def combine(per: dict[int, float]) -> dict[str, float]:
    """Per-energy similarities combined four ways (keys are the variant slugs)."""
    vals = list(per.values())
    wsum = sum(LEVEL_WEIGHT.get(e, 1.0) for e in per)
    return {"mean": sum(vals) / len(vals),
            "weighted": sum(LEVEL_WEIGHT.get(e, 1.0) * s for e, s in per.items()) / wsum,
            "min": min(vals),
            "high": per.get(3, 0.0)}


def read_text(path: Path) -> str:
    data = path.read_bytes()
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("cp1252", errors="replace")


def formula_key(formula: str) -> str:
    """Element counts in a fixed order, so C10H14ClN and C10Cl1H14N1 give the same key."""
    counts: dict[str, int] = defaultdict(int)
    for element, n in ELEMENT.findall(re.sub(r"\s", "", formula)):
        counts[element] += int(n) if n else 1
    return "".join(f"{el}{counts[el]}" for el in sorted(counts))


def block(key: str) -> str:
    return key.split("-")[0]


def level_name(level: int) -> str:
    return LEVEL_NAME.get(level, f"{level}")


def prepare(peaks: list[tuple[float, float]], min_rel: float) -> dict | None:
    """Peaks at or above min_rel of the base peak, sorted by m/z, square-root weights, plus the norm."""
    peaks = [(mz, i) for mz, i in peaks if i > 0]
    if not peaks:
        return None
    top = max(i for _, i in peaks)
    kept = sorted((mz, i / top) for mz, i in peaks if i / top >= min_rel)
    weights = [r ** 0.5 for _, r in kept]
    return {"mz": [mz for mz, _ in kept], "w": weights, "norm": math.sqrt(sum(w * w for w in weights))}


def parse_sdf(path: Path, min_rel: float) -> tuple[dict[str, dict], int]:
    """Compounds keyed by full InChIKey, one prepared spectrum per energy level; also returns unlabelled count."""
    compounds: dict[str, dict] = {}
    unlabelled = 0
    for chunk in read_text(path).split("$$$$"):
        fields: dict[str, list[str]] = {}
        tag = None
        for line in chunk.splitlines():
            match = TAG.match(line)
            if match:
                tag = match.group(1).strip().upper()
                fields[tag] = []
            elif tag is not None:
                if line.strip():
                    fields[tag].append(line.strip())
                else:
                    tag = None
        name = " ".join(fields.get("NAME", [])).strip()
        key = " ".join(fields.get("INCHIKEY", [])).strip()
        if not name or not key:
            continue
        suffix = SUFFIX.search(name)
        if suffix and suffix.group(1):
            level = VOLT_LEVEL.get(int(suffix.group(1)), int(suffix.group(1)))
        elif suffix and suffix.group(2):
            level = WORD_LEVEL[suffix.group(2).lower()]
        else:
            level = 0
            unlabelled += 1
        values = [float(x) for x in NUMBER.findall(" ".join(fields.get("MASS SPECTRAL PEAKS", [])))]
        spectrum = prepare(list(zip(values[0::2], values[1::2])), min_rel)
        formula = " ".join(fields.get("FORMULA", [])).strip()
        entry = compounds.setdefault(key, {"key": key, "name": SUFFIX.sub("", name).strip(),
                                           "formula": formula, "fkey": formula_key(formula), "spectra": {}})
        if spectrum:
            entry["spectra"][level] = spectrum
    return {k: c for k, c in compounds.items() if c["spectra"]}, unlabelled


def cosine(a: dict, b: dict, tol: float) -> float:
    """One-to-one peak matching by a two-pointer walk over m/z-sorted peaks."""
    am, aw, bm, bw = a["mz"], a["w"], b["mz"], b["w"]
    i = j = 0
    num = 0.0
    while i < len(am) and j < len(bm):
        diff = am[i] - bm[j]
        if abs(diff) <= tol:
            num += aw[i] * bw[j]
            i += 1
            j += 1
        elif diff < 0:
            i += 1
        else:
            j += 1
    return num / (a["norm"] * b["norm"]) if a["norm"] and b["norm"] else 0.0


def relation(q: dict, r: dict) -> str:
    """How the hit relates to the query: stereoisomer, isomer (same formula) or other compound."""
    if block(q["key"]) == block(r["key"]):
        return "stereoisomer"
    if q["fkey"] and q["fkey"] == r["fkey"]:
        return "isomer"
    return "other"


def search(q: dict, pool: list[dict], tol: float) -> dict | None:
    """Best hit by mean similarity over shared energies; None when nothing shares an energy."""
    levels = sorted(q["spectra"])
    best_per_level: dict[int, tuple[float, int]] = {}
    scored, combos = [], []
    for i, r in enumerate(pool):
        shared = [e for e in levels if e in r["spectra"]]
        if not shared:
            continue
        per = {e: cosine(q["spectra"][e], r["spectra"][e], tol) for e in shared}
        combo = combine(per)
        scored.append((combo["mean"], i))
        combos.append((combo, i))
        for e, s in per.items():
            if s > best_per_level.get(e, (-1.0, -1))[0]:
                best_per_level[e] = (s, i)
    if not scored:
        return None
    top = heapq.nlargest(2, scored)
    s1, i1 = top[0]
    hit = pool[i1]
    row = {
        "query": q["name"], "query_inchikey": q["key"], "formula": q["formula"],
        "top_hit": hit["name"], "top_hit_inchikey": hit["key"], "top_hit_formula": hit["formula"],
        "similarity": round(s1, 4), "relation": relation(q, hit),
        "isomer_in_reference": any(r["fkey"] == q["fkey"] for r in pool) if q["fkey"] else "",
        "hit2": pool[top[1][1]]["name"] if len(top) > 1 else "",
        "hit2_similarity": round(top[1][0], 4) if len(top) > 1 else "",
    }
    for _, slug in VARIANTS:
        s, i = max((c[slug], i) for c, i in combos)
        row[f"v_{slug}_score"] = round(s, 4)
        row[f"v_{slug}_top_hit"] = pool[i]["name"]
        row[f"v_{slug}_relation"] = relation(q, pool[i])
    for e in levels:
        name = level_name(e)
        row[f"sim_{name}_to_top_hit"] = round(cosine(q["spectra"][e], hit["spectra"][e], tol), 4) if e in hit["spectra"] else ""
        if e in best_per_level:
            s, i = best_per_level[e]
            row[f"top_hit_{name}_only"] = pool[i]["name"]
            row[f"relation_{name}_only"] = relation(q, pool[i])
    return row


def known_proxy(old: dict, new: dict, tol: float) -> tuple[list[dict], int]:
    """Compounds in both releases whose spectra changed: newer spectrum searched against the older library."""
    pool = list(old.values())
    by_block: dict[str, set[int]] = defaultdict(set)
    for i, c in enumerate(pool):
        by_block[block(c["key"])].add(i)
    rows, unchanged = [], 0
    for key, q in new.items():
        targets = by_block.get(block(key))
        if not targets:
            continue
        first = pool[min(targets)]
        if all(e in first["spectra"] and cosine(q["spectra"][e], first["spectra"][e], tol) > 0.9999 for e in q["spectra"]):
            unchanged += 1
            continue
        combos = []
        for i, r in enumerate(pool):
            per = {e: cosine(q["spectra"][e], r["spectra"][e], tol) for e in q["spectra"] if e in r["spectra"]}
            if per:
                combos.append((combine(per), i))
        row = {"query": q["name"], "query_inchikey": key}
        for _, slug in VARIANTS:
            vals = [(c[slug], i) for c, i in combos]
            top_s, top_i = max(vals)
            true_s = max(s for s, i in vals if i in targets)
            row[f"v_{slug}_top_score"] = round(top_s, 4)
            row[f"v_{slug}_correct"] = top_i in targets
            row[f"v_{slug}_true_score"] = round(true_s, 4)
            row[f"v_{slug}_top_relation"] = "same compound" if top_i in targets else relation(q, pool[top_i])
        rows.append(row)
    return rows, unchanged


def compare_variants(novel: list[dict], known: list[dict], out: Path) -> None:
    """Each combination at the threshold where a fixed share of known compounds would be called new."""
    print(f"\nEnergy weighting at an equal false-novelty budget (known side: {len(known)} compounds whose "
          f"spectra changed between releases; new side: {len(novel)} new compounds)")
    header = f"{'combination':14s} {'known right':>11s} {'known: isomer won':>17s}"
    for b in BUDGETS:
        header += f" | {'threshold':>9s} {'new: wrong':>10s} {'isomers':>8s}  @{b:.0%}"
    print(header)
    rows = []
    for label, slug in VARIANTS:
        tops = sorted(r[f"v_{slug}_top_score"] for r in known)
        right = sum(r[f"v_{slug}_correct"] is True for r in known)
        iso_won = sum(r[f"v_{slug}_top_relation"] in ("isomer", "stereoisomer") for r in known)
        line = f"{label:14s} {pct(right, len(known)):>11s} {pct(iso_won, len(known)):>17s}"
        row = {"combination": label, "known_n": len(known), "known_right": right, "known_isomer_won": iso_won,
               "novel_n": len(novel)}
        for b in BUDGETS:
            t = tops[min(len(tops) - 1, int(b * len(tops)))] if tops else float("nan")
            wrong = [r for r in novel if r[f"v_{slug}_score"] >= t]
            iso = sum(r[f"v_{slug}_relation"] in ("isomer", "stereoisomer") for r in wrong)
            line += f" | {t:9.3f} {pct(len(wrong), len(novel)):>10s} {pct(iso, len(wrong)):>8s}      "
            row[f"threshold_{b:.2f}"] = round(t, 4)
            row[f"novel_wrong_{b:.2f}"] = len(wrong)
            row[f"novel_wrong_isomers_{b:.2f}"] = iso
        print(line)
        rows.append(row)
    print("  'known right' = true compound ranked first; 'isomer won' = an isomer outranked it.")
    print("  Lower 'new: wrong' at the same budget = better at keeping new compounds from getting a name.")
    with open(out / "energy_variants.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pct(n: int, d: int) -> str:
    return f"{100 * n / d:.0f}%" if d else "n/a"


def summarise(rows: list[dict], label: str) -> None:
    n = len(rows)
    print(f"\n{label}: {n} queries")
    if not n:
        return
    for t in THRESHOLDS:
        hits = [r for r in rows if r["similarity"] >= t]
        rel = Counter(r["relation"] for r in hits)
        print(f"  wrong name with similarity >= {t:.2f}: {len(hits)} ({pct(len(hits), n)}); "
              f"isomers {rel['isomer']} ({pct(rel['isomer'], len(hits))}), "
              f"stereoisomers {rel['stereoisomer']} ({pct(rel['stereoisomer'], len(hits))})")
    have = sum(r["isomer_in_reference"] is True for r in rows)
    print(f"  an isomer (same formula) was in the reference for {have} ({pct(have, n)})")
    names = [level_name(e) for e in (1, 2, 3, 0)]
    present = [nm for nm in names if any(f"relation_{nm}_only" in r for r in rows)]
    if present:
        print("  searching one energy at a time, how often is the top hit an isomer or stereoisomer?")
        for nm in present:
            rel = Counter(r[f"relation_{nm}_only"] for r in rows if f"relation_{nm}_only" in r)
            used = sum(rel.values())
            iso = rel["isomer"] + rel["stereoisomer"]
            print(f"    {nm:>4}: {iso} of {used} ({pct(iso, used)})")
        iso_rows = [r for r in rows if r["relation"] in ("isomer", "stereoisomer")]
        parts = []
        for nm in present:
            vals = [r[f"sim_{nm}_to_top_hit"] for r in iso_rows if r.get(f"sim_{nm}_to_top_hit", "") != ""]
            if vals:
                parts.append(f"{nm} {statistics.median(vals):.2f}")
        if parts:
            print(f"  when the top hit is an isomer, median similarity to it by energy: {', '.join(parts)}")


def run_pool(queries: list[dict], pool: list[dict], tol: float, mode: str, started: float) -> tuple[list[dict], int]:
    rows, skipped = [], 0
    for k, q in enumerate(queries, 1):
        row = search(q, pool, tol)
        if row is None:
            skipped += 1
            continue
        row["mode"] = mode
        rows.append(row)
        if k % 200 == 0:
            print(f"  searched {k}/{len(queries)} ({time.time() - started:.0f} s)")
    if skipped:
        print(f"  skipped {skipped} queries that share no energy level with any reference spectrum")
    return rows, skipped


def split_pair(text: str) -> tuple[str, Path]:
    label, sep, value = text.partition("=")
    if not sep or not label.strip() or not value.strip():
        sys.exit(f"--release expects LABEL=PATH, got {text!r}")
    return label.strip(), Path(value.strip()).expanduser()


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", action="append", required=True, metavar="LABEL=PATH",
                    help="a DART-MS .SDF release; repeat, oldest first")
    ap.add_argument("--tol", type=float, default=0.01, help="m/z matching tolerance in Da (default 0.01)")
    ap.add_argument("--min-rel", type=float, default=0.01,
                    help="drop peaks below this fraction of the base peak (default 0.01)")
    ap.add_argument("--out", type=Path, default=Path("dart_out"))
    args = ap.parse_args(argv)

    releases = [split_pair(s) for s in args.release]
    missing = [str(p) for _, p in releases if not p.is_file()]
    if missing:
        sys.exit("File not found: " + ", ".join(missing))
    args.out.mkdir(parents=True, exist_ok=True)

    libs = []
    for label, path in releases:
        compounds, unlabelled = parse_sdf(path, args.min_rel)
        spectra = sum(len(c["spectra"]) for c in compounds.values())
        levels = Counter(e for c in compounds.values() for e in c["spectra"])
        groups = Counter(c["fkey"] for c in compounds.values() if c["fkey"])
        in_groups = sum(n for n in groups.values() if n > 1)
        print(f"{label}: {len(compounds)} compounds, {spectra} spectra; energies "
              + ", ".join(f"{level_name(e)} x{n}" for e, n in sorted(levels.items()))
              + f"; {in_groups} compounds share their formula with another ({pct(in_groups, len(compounds))})")
        if unlabelled:
            print(f"  warning: {unlabelled} spectra had no energy at the end of their name")
        libs.append((label, compounds))

    started = time.time()
    all_rows: list[dict] = []
    if len(libs) == 1:
        label, compounds = libs[0]
        pool = list(compounds.values())
        print(f"\nLeave-one-compound-out on {label}: each compound is searched against the other {len(pool) - 1}.")
        for k, q in enumerate(pool, 1):
            row = search(q, pool[:k - 1] + pool[k:], args.tol)
            if row is None:
                continue
            row["mode"] = f"leave-one-out {label}"
            all_rows.append(row)
            if k % 200 == 0:
                print(f"  searched {k}/{len(pool)} ({time.time() - started:.0f} s)")
        summarise(all_rows, f"Leave-one-compound-out ({label})")
    else:
        known_rows: list[dict] = []
        for (old_label, old), (new_label, new) in zip(libs, libs[1:]):
            proxy, unchanged = known_proxy(old, new, args.tol)
            print(f"\n{new_label} vs {old_label}, known side: {len(proxy)} compounds with changed spectra, "
                  f"{unchanged} carried over unchanged (left out)")
            known_rows += proxy
            old_blocks = {block(k) for k in old}
            new_blocks = {block(k) for k in new}
            ref = [c for k, c in new.items() if block(k) in old_blocks]
            added = [c for k, c in new.items() if block(k) not in old_blocks]
            rekeyed = sum(k not in old and block(k) in old_blocks for k in new)
            gone = sum(block(k) not in new_blocks for k in old)
            print(f"\n{new_label} vs {old_label}: {len(added)} new by structure, {rekeyed} re-keyed, "
                  f"{gone} removed; reference = {len(ref)} compounds, spectra from {new_label}")
            rows, _ = run_pool(added, ref, args.tol, f"{new_label} vs {old_label}", started)
            summarise(rows, f"New in {new_label}")
            all_rows += rows
        if len(libs) > 2:
            summarise(all_rows, "All releases pooled")
        if known_rows and all_rows:
            compare_variants(all_rows, known_rows, args.out)
            with open(args.out / "known_proxy.csv", "w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(known_rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(known_rows)

    print("\nHighest-scoring wrong names:")
    for r in sorted(all_rows, key=lambda r: -r["similarity"])[:20]:
        levels = ", ".join(f"{k[len('sim_'):-len('_to_top_hit')]} {r[k]}" for k in r if k.startswith("sim_") and r[k] != "")
        print(f"  {r['similarity']:.3f}  {r['query'][:38]:38s} -> {r['top_hit'][:38]:38s} [{r['relation']}; {levels}]")

    if all_rows:
        fields = list(dict.fromkeys(k for r in all_rows for k in r))
        with open(args.out / "queries.csv", "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted(all_rows, key=lambda r: -r["similarity"]))
    print(f"\nWrote {args.out / 'queries.csv'} in {time.time() - started:.0f} s.")


if __name__ == "__main__":
    main()
