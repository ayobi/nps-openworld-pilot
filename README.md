# Open-world identification of new psychoactive substances: pilot analyses

Pilot analyses for a planned review and benchmark on how mass-spectral library search behaves when a new psychoactive substance (NPS) is not in the library. Every number in the project doc's "Pilot results" section is reproduced by `run_all.sh`.

## What the pilot found (run of 8 November 2026)

- **New compounds get confident wrong names, mostly isomers.** Each of 13 SWGDRUG releases from 2017 to 2024 was searched against the library as it stood before it: of 1,176 new compounds, 19% (95% CI 17 to 22%) got a wrong name with a match factor of 800 or more and 5% at 900 or more. 67% and 84% of those wrong top hits share the new compound's molecular formula.
- **Known compounds lose to isomers too.** For 35 spectra with a counterpart from another lab, the true compound ranked first 74% of the time (95% CI 58 to 86%), and all 9 misses lost to a compound with the same formula.
- **Public data cannot measure the other error.** Only 17 of 3,800 SWGDRUG compounds have a second spectrum from a different contributor.
- **No match-factor threshold keeps both errors low,** and answering with the isomer group fixes known compounds but few new ones.
- **The pattern holds for DART-MS.** In the NIST DART-MS Forensics Database, 504 compounds were new by structure across four release comparisons (2021 to 2024); 10% got a wrong name at a similarity of 0.80 or more, and 48 of those 49 were isomers. Reweighting the collision energies does not remove the isomer problem.

## Layout

| Path | What it is |
| --- | --- |
| `get_data.sh` | Downloads the public libraries into `data/` and checks the DART-MS checksums |
| `run_all.sh` | Runs every analysis; writes `results/` and a log per step |
| `scripts/count_replicates.py` | Finds compounds with spectra from more than one contributor and flags copies |
| `scripts/pilot_temporal.py` | SWGDRUG temporal test, known side, threshold trade-off, entropy comparison, conformal flag, isomer-group rule, provisional classes |
| `scripts/dart_pilot.py` | DART-MS leave-one-compound-out and temporal tests, energy weighting |
| `scripts/pilot_fallback.py` | Class-level fallback: for compounds new at their release, how often the top hit's isomer group or drug class is right, by match factor (run after `pilot_temporal.py`) |
| `review/review_search.py` | Literature search for the review: PubMed, seed check, citation chasing, RIS merge |
| `review/pubmed_query.txt`, `review/seeds.csv` | The protocol's PubMed query and seed papers |
| `review/screening/` | Title-and-abstract screening decisions for the 1,010 records, each with its reason, and the 174 records taken forward |
| `review/extraction/` | Full-text extraction workbook (codebook, double coding, progress) and the abstract-level coding of all 174 candidates (`abstract_coding.csv`) |
| `review/figures/` | `make_figures.py` and the review's figures, drawn from `abstract_coding.csv` |
| `data/` | Libraries (not committed; see `data/README.md` for sources and licences) |

## Setup

```bash
python3 -m pip install -r requirements.txt
./get_data.sh
./run_all.sh
```

Python 3.9 or later. The SWGDRUG run takes about four minutes and the DART-MS leave-one-out about five.

## Licences of the data

SWGDRUG's library is copyrighted: the scripts download it from SWGDRUG and nothing derived from it beyond summary results belongs in this repository. The NIST DART-MS Forensics Database is under the NIST open licence. NIST's server has been returning errors for the 2025 and 2026 DART-MS releases; `run_all.sh` uses whatever downloaded.

## Caveats

- The match factor is the Stein and Scott (1994) composite, close to but not identical with NIST MS Search; confirm headline numbers in MS Search.
- Searches are spectrum-only, without retention data.
- SWGDRUG compounds are matched by CAS number and name (keeping primes, so 3- and 3'- stay distinct); DART-MS compounds by InChIKey connectivity.
- The known side is small (35 spectra in SWGDRUG; 41 reprocessed spectra as a proxy in DART-MS), so the calibrated flag and the energy-weighting comparison are directional.
- Compound classes are guessed from names until expert labels replace them.

## Review search

`review/review_search.py --email you@example.org --out review/run` runs the protocol's PubMed query, checks the seed papers are found, chases their citations through OpenAlex, runs a preprint search, merges any Scopus or Web of Science RIS exports given with `--import-ris`, and writes one deduplicated `combined.ris` for screening plus a `run_log.md` with the dates and counts to record.

`python review/figures/make_figures.py` redraws the review's figures from `review/extraction/abstract_coding.csv` (PDF and PNG, written to `review/figures/`).

## Citation

See `CITATION.cff`.
