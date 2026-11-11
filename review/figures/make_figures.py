"""Draw the figures of the NPS validation review from review/extraction/abstract_coding.csv.

Run from anywhere: python review/figures/make_figures.py
Writes PDF and PNG files next to this script, in review/figures/. The flow-diagram counts for identification and
title-and-abstract screening come from the search run logs and screening decisions in review/.
"""
import csv, math
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "review" / "extraction" / "abstract_coding.csv"
OUT = ROOT / "review" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

FAM = "Liberation Sans" if any(f.name == "Liberation Sans" for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
plt.rcParams.update({"font.family": FAM, "font.size": 8.5, "pdf.fonttype": 42})
BLUE, GREY, INK, SOFT, LINE = "#2c6fbb", "#c8c8bf", "#1f1f1f", "#5a5a5a", "#9a9a92"
CLASS = {"a": "#e8833a", "b": "#3a9e6f", "n": "#8e5cc4"}

rows = list(csv.DictReader(open(DATA, encoding="utf-8")))
elig = [r for r in rows if r["eligibility"] == "eligible"]
tested = lambda r: r["tested_absent_compounds"] == "yes"

def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf"); fig.savefig(OUT / f"{name}.png", dpi=300); plt.close(fig)

def wilson(k, n, z=1.96):
    p = k / n; c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h

# Figure 1: the five evaluation designs
def fig_designs():
    fig = plt.figure(figsize=(6.8, 4.6)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    designs = [
        ("Random split of spectra", "other spectra of the same compounds", list("AABCD"), "aaabb", [("A", "a", 0), ("C", "b", 0)]),
        ("External or real samples", "new samples of known compounds", list("ABCDE"), "aabba", [("B", "a", 0), ("D", "b", 0)]),
        ("Compound-disjoint", "new compounds of known classes", list("ABCD"), "aabb", [("E", "a", 1), ("F", "b", 1)]),
        ("Class-disjoint", "a whole compound class held out", list("ABCD"), "aabb", [("G", "n", 1), ("H", "n", 1)]),
        ("Temporal", "compounds first reported after a cut-off date", list("ABCD"), "aabb", [("E", "a", 1), ("G", "n", 1)]),
    ]
    ax.text(50, 95, "Reference library or training data", ha="center", color=SOFT, fontsize=8)
    ax.text(77.5, 95, "Test set", ha="center", color=SOFT, fontsize=8)
    ax.text(92, 97.5, "Unseen\ncompounds\nin the test set?", ha="center", va="top", color=SOFT, fontsize=7.2, linespacing=1.15)
    for k, (name, desc, ref, cls, test) in enumerate(designs):
        y = 82 - k * 16
        ax.text(2, y + 2.2, name, weight="bold", color=INK, fontsize=8.6, va="center")
        ax.text(2, y - 3.0, desc, color=SOFT, fontsize=7.6, va="center")
        xs = [38 + 5.4 * i for i in range(len(ref))]
        ax.scatter(xs, [y] * len(ref), s=230, c=[CLASS[c] for c in cls], edgecolors="white", linewidths=1.0, zorder=2)
        for x, letter in zip(xs, ref):
            ax.text(x, y, letter, ha="center", va="center", color="white", fontsize=7.6, weight="bold", zorder=3)
        ax.annotate("", xy=(71.5, y), xytext=(66.5, y), arrowprops=dict(arrowstyle="->", color=LINE, lw=0.9))
        for i, (letter, c, new) in enumerate(test):
            x = 75 + 5.4 * i
            ax.scatter([x], [y], s=230, c=[CLASS[c]], edgecolors=INK if new else "white", linewidths=2.0 if new else 1.0, zorder=2)
            ax.text(x, y, letter, ha="center", va="center", color="white", fontsize=7.6, weight="bold", zorder=3)
        unseen = any(t[2] for t in test)
        ax.text(92, y, "yes" if unseen else "no", ha="center", va="center", fontsize=8.6,
                color=BLUE if unseen else SOFT, weight="bold" if unseen else "normal")
        if name == "Temporal":
            ax.text(49, y - 5.6, "first reported before the cut-off", ha="center", color=SOFT, fontsize=6.8)
            ax.text(77.7, y - 5.6, "after", ha="center", color=SOFT, fontsize=6.8)
        if k < len(designs) - 1:
            ax.plot([2, 98], [y - 8, y - 8], color="#e4e4de", lw=0.6)
    ly = 5
    for i, (lab, c) in enumerate([("class 1", "a"), ("class 2", "b"), ("a class absent from the reference", "n")]):
        x = [4, 16, 28][i]
        ax.scatter([x], [ly], s=60, c=[CLASS[c]], edgecolors="white", zorder=2); ax.text(x + 1.8, ly, lab, va="center", color=SOFT, fontsize=7.4)
    ax.scatter([66], [ly], s=60, c=["white"], edgecolors=INK, linewidths=1.8); ax.text(67.8, ly, "compound absent from the reference", va="center", color=SOFT, fontsize=7.4)
    save(fig, "fig1_designs")

# Figure 2: flow of records
def fig_flow():
    fig = plt.figure(figsize=(6.8, 6.6)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    def box(x, y, w, h, title, body=""):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc="#f3f3ef", ec=LINE, lw=0.8))
        ax.text(x + 1.8, y + h - 2.0, title, weight="bold", va="top", color=INK, fontsize=8.4)
        if body: ax.text(x + 1.8, y + h - 5.6, body, va="top", color=SOFT, fontsize=7.6, linespacing=1.45)
    def arrow(x1, y1, x2, y2): ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", color=LINE, lw=0.9))
    n_ft = sum(r["coded_from"] == "full text" for r in rows); n_ab = len(rows) - n_ft
    e_ft = sum(r["coded_from"] == "full text" and r["eligibility"] == "eligible" for r in rows); e_ab = len(elig) - e_ft
    cnt = Counter((r["coded_from"], r["eligibility"]) for r in rows)
    x_ft_excl = n_ft - e_ft
    ident = {"PubMed": 370, "Scopus": 472, "References of the seed studies": 140, "Studies citing the seed studies": 301, "Preprint search": 216, "Seed studies": 5}
    dup, early, screened = 483, 11, 1010
    excl = [("Not about NPS", 299), ("Case report or survey", 181), ("Review or opinion", 166), ("No identification step", 95),
            ("No mass-spectrometry arm", 51), ("Toxicity or QSAR only", 29), ("Quantification only", 15)]
    total = sum(ident.values()); assert total - dup - early == screened and screened - sum(n for _, n in excl) == len(rows)
    for lab, y in [("Identification", 88), ("Screening", 64.5), ("Eligibility", 39), ("Included", 10.5)]:
        ax.text(2.2, y, lab, rotation=90, ha="center", va="center", color=SOFT, fontsize=8)
    box(6, 79, 46, 18, f"Records identified (n = {total:,})",
        f"PubMed {ident['PubMed']} · Scopus {ident['Scopus']}\nReferences of the seed studies {ident['References of the seed studies']}\n"
        f"Studies citing the seed studies {ident['Studies citing the seed studies']}\nPreprint search {ident['Preprint search']} · Seed studies {ident['Seed studies']}")
    box(57, 82, 41, 12, f"Removed before screening (n = {dup + early})", f"Duplicates {dup}\nPublished before 2010 {early}")
    arrow(52, 88, 57, 88)
    box(6, 60, 46, 9, f"Screened on title and abstract", f"n = {screened:,}")
    box(57, 47, 41, 27, f"Excluded (n = {sum(n for _, n in excl)})", "\n".join(f"{a} {n}" for a, n in excl))
    arrow(52, 64.5, 57, 64.5)
    box(6, 33, 46, 12, f"Assessed for eligibility (n = {len(rows)})", f"Full text {n_ft} · abstract only, paywalled {n_ab}")
    box(57, 22, 41, 21, f"Not eligible (n = {len(rows) - len(elig)})",
        f"From full text ({x_ft_excl}): no performance figure\non a test set or case report\n"
        f"From abstracts ({n_ab - e_ab}): not eligible {cnt[('abstract', 'not eligible')]},\nno figure in the abstract {cnt[('abstract', 'unclear')]}, no abstract {cnt[('abstract', 'no abstract')]}")
    arrow(52, 39, 57, 39)
    box(6, 5, 46, 11, f"Eligible studies (n = {len(elig)})", f"{e_ft} coded from full text · {e_ab} from abstracts")
    arrow(29, 79, 29, 69); arrow(29, 60, 29, 45); arrow(29, 33, 29, 16)
    save(fig, "fig2_flow")

# Figure 3: split type against unseen-compound testing
def fig_split():
    groups = [("External or real samples only", ("external", "external only")), ("Not stated", ("not stated", "not reported")),
              ("Compound-disjoint", ("compound-disjoint",)), ("Random split", ("random", "random spectra")), ("Temporal", ("temporal",)), ("Class-disjoint", ("class-disjoint",))]
    cats = [(lab, sum(tested(r) for r in elig if r["split_type"] in keys), sum(not tested(r) for r in elig if r["split_type"] in keys)) for lab, keys in groups]
    t_all = sum(c[1] for c in cats); n_all = sum(c[1] + c[2] for c in cats); assert n_all == len(elig)
    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    ys = list(range(len(cats)))[::-1]
    for y, (lab, t, u) in zip(ys, cats):
        ax.barh(y, t, color=BLUE, height=0.58); ax.barh(y, u, left=t, color=GREY, height=0.58)
        ax.text(t + u + 0.8, y, "none" if t + u == 0 else f"{t} of {t + u}", va="center", color=SOFT, fontsize=8.5)
    ax.set_yticks(ys); ax.set_yticklabels([c[0] for c in cats], color=INK); ax.set_xlim(0, 62); ax.set_xticks([])
    for s in ("top", "right", "bottom"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#999999"); ax.tick_params(axis="y", length=0, pad=6)
    fig.text(0.012, 0.965, "Few evaluations test compounds absent from the reference", fontsize=11, weight="bold", color=INK, va="top")
    fig.text(0.012, 0.885, f"{t_all} of {n_all} eligible studies ({round(100 * t_all / n_all)}%) tested any compound absent from the reference library or training data", fontsize=8.5, color=SOFT, va="top")
    fig.patches.extend([plt.Rectangle((0.012, 0.79), 0.014, 0.035, color=BLUE, transform=fig.transFigure), plt.Rectangle((0.16, 0.79), 0.014, 0.035, color=GREY, transform=fig.transFigure)])
    fig.text(0.032, 0.807, "tested", fontsize=8.5, color=SOFT, va="center"); fig.text(0.18, 0.807, "not tested or unclear", fontsize=8.5, color=SOFT, va="center")
    plt.subplots_adjust(left=0.33, right=0.97, top=0.74, bottom=0.04)
    save(fig, "fig3_split")

# Figure 4: eligible studies by publication year
def fig_years():
    years = list(range(2010, 2027))
    t = [sum(tested(r) for r in elig if int(r["year"]) == y) for y in years]
    o = [sum(not tested(r) for r in elig if int(r["year"]) == y) for y in years]
    pre_t, pre_n = sum(a for y, a in zip(years, t) if y < 2020), sum(a + b for y, a, b in zip(years, t, o) if y < 2020)
    post_t, post_n = sum(a for y, a in zip(years, t) if y >= 2020), sum(a + b for y, a, b in zip(years, t, o) if y >= 2020)
    lo, hi = wilson(post_t, post_n)
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    ax.bar(years, t, color=BLUE, width=0.7); ax.bar(years, o, bottom=t, color=GREY, width=0.7)
    for y, a, b in zip(years, t, o):
        if a + b: ax.text(y, a + b + 0.25, str(a + b), ha="center", color=SOFT, fontsize=7.5)
    ax.axvline(2019.5, color="#999999", lw=0.8, ls="--")
    ax.text(2014.5, 13.4, f"before 2020: {pre_t} of {pre_n} tested", ha="center", color=SOFT, fontsize=8)
    ax.text(2023, 13.4, f"2020 onwards: {post_t} of {post_n} ({round(100 * post_t / post_n)}%)", ha="center", color=SOFT, fontsize=8)
    ax.set_xticks(years); ax.set_xticklabels([str(y) if y % 2 == 0 else "" for y in years], color=INK)
    ax.set_ylim(0, 14.5); ax.set_yticks([0, 5, 10]); ax.tick_params(colors=SOFT, length=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color("#999999")
    ax.set_ylabel("eligible studies", color=SOFT)
    fig.text(0.012, 0.975, "Tests on unseen compounds became more common after 2020, but remain a minority", fontsize=10.5, weight="bold", color=INK, va="top")
    fig.text(0.012, 0.905, "Eligible studies by year of publication (2026 up to November); blue: tested any compound absent from the reference", fontsize=8.2, color=SOFT, va="top")
    plt.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.1)
    save(fig, "fig4_years")
    return pre_t, pre_n, post_t, post_n, lo, hi

# Figure S1: what each search source contributed
def fig_sources():
    order = [("import: scopus.ris", "Scopus"), ("pubmed", "PubMed"), ("forward citation", "Studies citing the seed studies"),
             ("backward citation", "References of the seed studies"), ("seed", "Seed studies"), ("preprint search", "Preprint search")]
    src = lambda r: {s.strip() for s in r["sources"].split(";") if s.strip()}
    vals = [(lab, sum(key in src(r) for r in elig), sum(key in src(r) for r in rows)) for key, lab in order]
    cit_only = sum(src(r) <= {"forward citation", "backward citation"} for r in elig)
    fig, ax = plt.subplots(figsize=(6.6, 2.8))
    ys = list(range(len(vals)))[::-1]
    for y, (lab, e, c) in zip(ys, vals):
        ax.barh(y, e, color=BLUE, height=0.55)
        ax.text(e + 0.8, y, f"{e} of {c} candidates", va="center", color=SOFT, fontsize=8)
    ax.set_yticks(ys); ax.set_yticklabels([v[0] for v in vals], color=INK); ax.set_xlim(0, 112); ax.set_xticks([])
    for s in ("top", "right", "bottom"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#999999"); ax.tick_params(axis="y", length=0, pad=6)
    fig.text(0.012, 0.96, "Eligible studies by search source", fontsize=10.5, weight="bold", color=INK, va="top")
    fig.text(0.012, 0.875, f"A study can come from several sources; {cit_only} eligible studies were found only by citation chasing", fontsize=8.2, color=SOFT, va="top")
    plt.subplots_adjust(left=0.33, right=0.97, top=0.76, bottom=0.04)
    save(fig, "figS1_sources")
    return cit_only

fig_designs(); fig_flow(); fig_split(); stats = fig_years(); cit = fig_sources()
print("before 2020 %d of %d; 2020 onwards %d of %d, 95%% CI %.0f-%.0f%%" % (stats[0], stats[1], stats[2], stats[3], 100 * stats[4], 100 * stats[5]))
print("eligible found only by citation chasing:", cit)
