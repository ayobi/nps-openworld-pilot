#!/usr/bin/env bash
# Download the public spectral libraries into data/. Safe to rerun: existing files are kept.
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p data/swgdrug data/dartms

echo "SWGDRUG MS library v3.14"
cd data/swgdrug
for f in JCAMP_01222025.zip Changelog_01222025.xls; do
  [ -f "$f" ] || curl -fLO "https://www.swgdrug.org/MS/$f" || echo "  could not download $f"
done
[ -f compound_list.xlsx ] || curl -fL -o compound_list.xlsx "https://www.swgdrug.org/MS/Compound%20List%2001222025.xlsx" || echo "  could not download the compound list"
[ -d jcamp ] || unzip -oq JCAMP_01222025.zip -d jcamp

echo "NIST DART-MS Forensics Database"
cd ../dartms
for r in 210915_He_Pos_Library_Dragonfly 220120_He_Pos_Library_Earthworm 230104_He_Pos_Library_Grasshopper \
         230919_HePos_Hornet_v1 240319_HePos_Inchworm_v1 250101_HePos_Joro 251229_HePos_Katydid 260414_HePos_Locust; do
  [ -f "$r.zip" ] || curl -fLO "https://data.nist.gov/od/ds/mds2-2313/$r.zip" || echo "  could not download $r"
done
[ -f 260914_HePos_Monarch.zip ] || curl -fLO "https://data.nist.gov/od/ds/mds2-2313/260914_HePos_Monarch.zip" || curl -fLO "https://data.nist.gov/od/ds/mds2-4327/260914_HePos_Monarch.zip" || echo "  could not download 260914_HePos_Monarch"

sha() { if command -v shasum > /dev/null; then shasum -a 256 "$1" | cut -d' ' -f1; else sha256sum "$1" | cut -d' ' -f1; fi; }
while read -r hash file; do
  [ -f "$file" ] || continue
  if [ "$(sha "$file")" = "$hash" ]; then echo "  checksum ok   $file"; else echo "  CHECKSUM MISMATCH  $file"; fi
done < ../dartms.sha256

for z in [0-9]*.zip; do
  [ -f "$z" ] || continue
  [ -d "${z%.zip}" ] || unzip -oq "$z" -d "${z%.zip}"
done
echo "Done. Releases unpacked: $(ls -d [0-9]*_He*Pos* 2>/dev/null | grep -vc '\.zip$')"
