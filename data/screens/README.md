# Loss-of-function screens

Per-gene extracts used by `score.py` for the proposed Functional column. Made by `fetch_screens.py` on 11 Nov 2026; large raw downloads are in `raw/` (not committed).

- `extract_jump_crispr_u2os.tsv`: JUMP Cell Painting, CRISPR knockouts in U2OS, corrected p of phenotypic activity (Zenodo 14642166). Hit: p < 0.05.
- `extract_periscope.tsv`: PERISCOPE (Ramezani et al. 2025), HeLa DMEM, HeLa HPLM, A549; 1 = in the published whole-cell or compartment hit list. Tested genes: `periscope_HeLa_DMEM_genes_tested.txt`.
- `extract_funk2022_hela.tsv`: Funk et al. 2022, Table S2A (5,049 fitness genes, HeLa). Hit: 3+ of 32 features at FDR < 0.05, or mitotic index above non-targeting mean + 3 SD.
- `extract_mitocheck_hela.tsv`: MitoCheck (Neumann et al. 2010), IDR idr0013 screen A annotation. Hit: phenotype reproducible across siRNA replicates, or validated mitotic hit.
- `extract_replogle_k562_gwps.tsv`, `extract_replogle_rpe1.tsv`: Replogle et al. 2022 Perturb-seq bulk obs. Hit: energy-distance p < 0.001 and 50+ Anderson-Darling DE genes.
- `extract_orcs_screens.tsv`, `extract_orcs_genes.tsv`: BioGRID ORCS 2.0.18, screen metadata (from each screen page) and per-gene hit counts in loss-of-function screens (knockout or CRISPRi) in U2OS, HeLa, and all lines.
- `extract_hgnc.tsv`: HGNC symbols, previous symbols, aliases, Entrez, Ensembl, UniProt, to match old gene names.
- Inputs copied from GitHub: `periscope__*` hit lists and READMEs, `funk2022_*` tables, `jump_crispr_metadata.csv.gz`.
