"""Class-level fallback in the SWGDRUG temporal test.

For every compound that was new at its release (absent from the library it was searched against),
asks whether answering with the top hit's isomer group (same formula) or drug class, instead of its
name, would have been right, by score band. Reads the per-query output of pilot_temporal.py (SWGDRUG,
match factor) or dart_pilot.py (DART-MS, similarity; detected from the columns) and reuses the
provisional, name-based classes of pilot_temporal.py (to be replaced by expert labels).

Usage: python scripts/pilot_fallback.py [results/series/queries.csv] [out_dir]
       python scripts/pilot_fallback.py results/dart_temporal/queries.csv results/dart_temporal
"""
import csv, re, sys
from pathlib import Path
import importlib.util
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
QUERIES = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "series" / "queries.csv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("pilot_temporal", ROOT / "scripts" / "pilot_temporal.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)

def drug_class(name: str):
    n = (name or "").lower().replace(".alpha.", "alpha")
    for label, pattern in pt.CLASSES:
        if re.search(pattern, n):
            return label
    return None

reader = csv.DictReader(open(QUERIES, encoding="utf-8"))
DART = "similarity" in (reader.fieldnames or [])
if DART:  # dart_pilot.py temporal output: every row is a compound new at its release
    rows = list(reader)
    for r in rows:
        r["mf"] = float(r["similarity"] or 0)
        r["isomer"] = r["relation"] in ("isomer", "stereoisomer")
        r["q_class"], r["hit_class"] = drug_class(r["query"]), drug_class(r["top_hit"])
    BANDS = [("0.95 or more", 0.95, 2), ("0.90 to 0.95", 0.90, 0.95), ("0.80 to 0.90", 0.80, 0.90), ("0.70 to 0.80", 0.70, 0.80), ("below 0.70", 0, 0.70)]
    CONFIDENT, SCORE, LIBRARY, PREFIX = 0.80, "similarity", "DART-MS", "dart_"
else:  # pilot_temporal.py output
    rows = [r for r in reader if r["status"] == "novel"]
    for r in rows:
        r["mf"] = float(r["top_mf"] or 0)
        r["isomer"] = r["same_formula"] == "True"
        r["q_class"], r["hit_class"] = drug_class(r["name"]), drug_class(r["top_hit_name"])
    BANDS = [("900 or more", 900, 10_000), ("800 to 899", 800, 900), ("700 to 799", 700, 800), ("600 to 699", 600, 700), ("below 600", 0, 600)]
    CONFIDENT, SCORE, LIBRARY, PREFIX = 800, "match factor", "SWGDRUG", ""
table = []
for label, lo, hi in BANDS:
    band = [r for r in rows if lo <= r["mf"] < hi]
    assess = [r for r in band if r["q_class"]]
    table.append({"score_band": label, "novel": len(band), "isomer_group_right": sum(r["isomer"] for r in band),
                  "class_assessable": len(assess), "isomer_group_right_assessable": sum(r["isomer"] for r in assess),
                  "class_right_only": sum((not r["isomer"]) and r["q_class"] == r["hit_class"] for r in assess),
                  "neither": sum((not r["isomer"]) and r["q_class"] != r["hit_class"] for r in assess)})
with open(OUT / f"{PREFIX}fallback_by_band.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(table[0])); w.writeheader(); w.writerows(table)

conf = [r for r in rows if r["mf"] >= CONFIDENT]; ca = [r for r in conf if r["q_class"]]
print(f"{LIBRARY}: novel compounds: {len(rows)}; confident wrong names ({SCORE} {CONFIDENT}+): {len(conf)} ({len(conf) / len(rows):.0%})")
print(f"  top hit an isomer: {sum(r['isomer'] for r in conf)} ({sum(r['isomer'] for r in conf) / len(conf):.0%})")
print(f"  class right, where a class is known: {sum(r['q_class'] == r['hit_class'] for r in ca)} of {len(ca)}")
print(f"  isomer group or class right, where a class is known: {sum(r['isomer'] or r['q_class'] == r['hit_class'] for r in ca)} of {len(ca)}")

BLUE, LIGHT, GREY, INK, SOFT = "#2c6fbb", "#9dbfe3", "#c8c8bf", "#1f1f1f", "#5a5a5a"
fig, ax = plt.subplots(figsize=(6.6, 3.2))
ys = list(range(len(table)))[::-1]
for y, t in zip(ys, table):
    n = t["class_assessable"]; a, b, c = (t["isomer_group_right_assessable"] / n, t["class_right_only"] / n, t["neither"] / n)
    ax.barh(y, a, color=BLUE, height=0.6); ax.barh(y, b, left=a, color=LIGHT, height=0.6); ax.barh(y, c, left=a + b, color=GREY, height=0.6)
    ax.text(1.02, y, f"n = {n}", va="center", color=SOFT, fontsize=8)
ax.set_yticks(ys); ax.set_yticklabels([t["score_band"] for t in table], color=INK); ax.set_xlim(0, 1.12)
ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0%", "50%", "100%"], color=SOFT)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_ylabel(f"top {SCORE}", color=SOFT)
fig.text(0.012, 0.975, "For a new compound, a high score points to its class, never its name", fontsize=10.5, weight="bold", color=INK, va="top")
fig.text(0.012, 0.905, f"Compounds new at their {LIBRARY} release with a name-based class; the exact name is wrong for every one", fontsize=8, color=SOFT, va="top")
for x, col, lab in ((0.25, BLUE, "isomer group right"), (0.47, LIGHT, "class right, not an isomer"), (0.75, GREY, "neither")):
    fig.patches.append(plt.Rectangle((x, 0.80), 0.014, 0.035, color=col, transform=fig.transFigure))
    fig.text(x + 0.02, 0.817, lab, fontsize=8, color=SOFT, va="center")
plt.subplots_adjust(left=0.2, right=0.93, top=0.74, bottom=0.12)
fig.savefig(OUT / f"{PREFIX}fallback_by_band.png", dpi=300); fig.savefig(OUT / f"{PREFIX}fallback_by_band.pdf")
