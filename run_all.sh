#!/usr/bin/env bash
# Reproduce every pilot number in the project doc. Logs go to results/*.log.
set -euo pipefail
cd "$(dirname "$0")"
S=data/swgdrug
LIB="$S/jcamp/SWGDRUG 314.HPJ"
[ -f "$LIB" ] || { echo "Run ./get_data.sh first ($LIB is missing)"; exit 1; }
mkdir -p results

echo "1/5 Replicate spectra in SWGDRUG"
python3 scripts/count_replicates.py --lib swgdrug="$LIB" --out results/replicates | tee results/replicates.log

echo "2/5 SWGDRUG v3.14 additions against the earlier library"
python3 scripts/pilot_temporal.py --lib "$LIB" --compound-list "$S/compound_list.xlsx" --out results/pilot_v314 | tee results/pilot_v314.log

echo "3/5 SWGDRUG release by release, known side, trade-off, score comparison, isomer rule, classes"
python3 scripts/pilot_temporal.py --lib "$LIB" --changelog "$S/Changelog_01222025.xls" --replicates --plot --out results/series | tee results/series.log

releases=()
for d in $(ls -d data/dartms/[0-9]*_He*Pos* 2>/dev/null | grep -v '\.zip$' | sort); do
  s=$(find "$d" -iname '*.sdf' -not -path '*__MACOSX*' | head -1)
  [ -n "$s" ] && releases+=(--release "$(basename "$d" | cut -d_ -f1)=$s")
done
# releases holds pairs ("--release" "label=path"), so [0] and [1] are the first release
if [ ${#releases[@]} -ge 2 ]; then
  echo "4/5 DART-MS leave-one-compound-out on the first release"
  python3 scripts/dart_pilot.py "${releases[0]}" "${releases[1]}" --out results/dart_loo | tee results/dart_loo.log
fi
if [ ${#releases[@]} -ge 4 ]; then
  echo "5/5 DART-MS release by release and energy weighting"
  python3 scripts/dart_pilot.py "${releases[@]}" --out results/dart_temporal | tee results/dart_temporal.log
fi
echo "Done: results/"
