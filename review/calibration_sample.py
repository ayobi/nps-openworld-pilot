#!/usr/bin/env python3
"""Draw the protocol's calibration sample from a search run, reproducibly.

Reads combined.csv from a run folder (written by review_search.py), draws a
random sample (default 10%) with a fixed seed, and writes, in the same folder:

  calibration.ris            both reviewers screen these, blind, then compare
  remaining_part1.ris        the rest, split in two for divided screening
  remaining_part2.ris
  calibration_sample.csv     which records were drawn, for the record

The same run folder, fraction and seed always give the same sample, so the
seed belongs in the protocol.

Example:
  python3 review/calibration_sample.py --run review/run2 --seed 20261109
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from review_search import FIELDS, write_ris  # noqa: E402


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=Path, required=True, help="run folder holding combined.csv")
    ap.add_argument("--fraction", type=float, default=0.10, help="share drawn for calibration (default 0.10)")
    ap.add_argument("--seed", type=int, required=True, help="random seed; record it in the protocol")
    args = ap.parse_args(argv)

    source = args.run / "combined.csv"
    if not source.is_file():
        sys.exit(f"{source} not found; run review_search.py first")
    with open(source, newline="", encoding="utf-8") as fh:
        records = list(csv.DictReader(fh))
    if not records:
        sys.exit(f"{source} has no records")

    k = math.ceil(args.fraction * len(records))
    rng = random.Random(args.seed)
    drawn = sorted(rng.sample(range(len(records)), k))
    drawn_set = set(drawn)
    calibration = [records[i] for i in drawn]
    rest = [r for i, r in enumerate(records) if i not in drawn_set]
    rng.shuffle(rest)
    half = math.ceil(len(rest) / 2)

    write_ris(args.run / "calibration.ris", calibration)
    write_ris(args.run / "remaining_part1.ris", rest[:half])
    write_ris(args.run / "remaining_part2.ris", rest[half:])
    with open(args.run / "calibration_sample.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["row"] + FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for i in drawn:
            writer.writerow({"row": i + 1, **records[i]})

    print(f"{len(records)} records: {k} drawn for calibration (seed {args.seed}), "
          f"{half} and {len(rest) - half} in the two remaining parts")
    print(f"Wrote calibration.ris, remaining_part1.ris, remaining_part2.ris and calibration_sample.csv in {args.run}/")


if __name__ == "__main__":
    main()
