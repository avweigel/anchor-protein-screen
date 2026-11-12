# Anchor protein screen

A scored list of every reviewed human protein against the Alpha Cell anchor protein criteria, and a page to explore it.

Anchor proteins are a handful of proteins, 5 to 10, tagged in the U2-OS line so the same molecule can be seen by light microscopy and found in the EM volume of the same cell. Nothing in this list is excluded and nothing is a must-have. The score is a guide for the working group, not the decision.

## Explore

The page is served from this repository with GitHub Pages: open `index.html` on the Pages URL. It loads `anchor_scores.tsv` and lets you filter by compartment family, by each criterion and by total score, sort any column, open a protein to see where every score comes from, and download the filtered rows.

## Files

- `index.html`: the explorer, one file, no build step.
- `anchor_scores.tsv`: the table, one row per UniProt entry (20,661 rows).
- `score.py`: builds the table from `data/`. Run `python3 score.py` in this folder.
- `fetch_screens.py`: downloads the loss-of-function screens and writes small per-gene extracts to `data/screens/` (needs network, pandas, openpyxl, pyarrow, h5py; big raw files go to `data/screens/raw/`, not committed). Only needed to refresh the screens; `score.py` reads the extracts.
- `hpa_em_lookup.csv`: how each HPA location term maps to an EM-visible structure (yes, partial, no).
- `data/`: the inputs, see `data/README.md` for sources and versions.

## Scores

Six criteria, each 0 to 2, summed with equal weight to a total out of 12.

| Criterion | 2 | 1 | 0 |
|---|---|---|---|
| Structure identity | clear EM structure | broad class (vesicles, cytoskeleton, chromosome) | diffuse or no location |
| Sources agree | HPA and UniProt in the same family | only one source has a location | disagree |
| HPA antibody | Enhanced or Supported | Approved | Uncertain or no IF data |
| Expressed in U2-OS | 10 nTPM or more | 1 nTPM or more | below 1 |
| Copy number known | in Beck 2011 | not measured | |
| Tagged line exists | OpenCell or Allen line | none | |

A seventh column, **Functional (proposed)**, was suggested at the Space WG on 11 Nov 2026: does losing the protein change the cell in published loss-of-function screens? It is shown in the table but is not in `score_total`; `score_total_with_function` (out of 14) includes it, and the page has a switch for it.

| Functional (proposed) | clear imaging phenotype in any one screen (JUMP p < 0.001, PERISCOPE 2+ conditions, Funk 6+ features or mitotic arrest, MitoCheck validated) | any other phenotype | tested, no phenotype (blank if in no screen) |

Screens: JUMP Cell Painting (CRISPR, U2-OS), PERISCOPE (CRISPR, HeLa and A549), Funk 2022 (CRISPR, HeLa), MitoCheck (RNAi, HeLa), Perturb-seq (CRISPRi, K562 and RPE1), DepMap 26Q3 (U2-OS fitness), BioGRID ORCS 2.0.18 (U2-OS and HeLa CRISPR screens). Per-screen results are in the `fx_*` columns.

Details and known limits are in the page under "How the score is made".
