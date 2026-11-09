# Data

Nothing in this folder except this file and `dartms.sha256` is committed. Run `./get_data.sh` from the repository root to fill it.

| Folder | Source | Licence | Notes |
| --- | --- | --- | --- |
| `swgdrug/` | [SWGDRUG Mass Spectral Library](https://www.swgdrug.org/ms.htm) v3.14 (22 Jan 2025): JCAMP-DX file, changelog, compound list | Copyright SWGDRUG, all rights reserved | Never commit or redistribute; download from SWGDRUG |
| `dartms/` | [NIST DART-MS Forensics Database](https://data.nist.gov/od/id/mds2-2313), releases from Sep 2021 on | [NIST open licence](https://www.nist.gov/open/license) | Zips are checked against the SHA-256 sums NIST publishes (`dartms.sha256`) |

The Cayman Spectral Library needs a free registration and ships only in NIST MS Search and Agilent formats, so it is not part of the automated download.
