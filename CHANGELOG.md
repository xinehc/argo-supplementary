# Changelog
## 2025-05-01
### Changed
- Fix diamond evalue to 1e-25 (from 1e-15).


## 2024-10-22
### Changed
- Make plasmid identification cutoffs stringent by setting `--conservative` and ensuring *zero* USCG found in the first round of `geNomad` filtering.
- Account for *circular* topology of *complete genome* and *chromosome* level GTDB assemblies when extracting ARG-containing sequences.