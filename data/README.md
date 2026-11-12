# Anchor proteome scoring, data folder

Downloaded 5 Nov 2026 for scoring every human protein with subcellular data against the Alpha Cell anchor criteria.

- subcellular_location.tsv: HPA v25.1 subcellular location, 13,603 genes (proteinatlas.org/about/download)
- proteinatlas.tsv: HPA v25.1 full per-gene table (protein class, antibodies, UniProt, Ensembl)
- uniprot_human_reviewed.tsv: UniProt reviewed human proteome, 20,429 entries, with signal peptide, transmembrane, topological domain, lipidation, subcellular location (rest.uniprot.org stream)
- opencell_lines.json: all OpenCell tagged lines (HEK293T), with terminus, copy number, localization categories (opencell.sf.czbiohub.org/api/lines)
- depmap_26Q3_gene_effect_U2OS_and_summary.tsv: Chronos gene effect for U2OS (ACH-000364) and the fraction of lines below -0.5, extracted from DepMap Public 26Q3 CRISPRGeneEffect.csv
- Model.csv: DepMap model metadata
- allen_cell_collection_lines.tsv: Allen Cell Collection hiPSC lines, from cell-catalog.allencell.org

- hpa_rna_U2OS.tsv: U2OS rows (TPM, pTPM, nTPM) extracted from HPA rna_celline.tsv (the full 1 GB file is in Downloads as rna_celline.tsv.zip)
- beck2011_U2OS_copies_per_cell.tsv: Table S1 of Beck et al. 2011, Mol Syst Biol (doi 10.1038/msb.2011.82), 7,309 proteins, copies per cell in U2OS
- beck2011_U2OS_mitotic_copy_ratio.tsv: Table S2 of the same paper, copy-number ratio control vs mitotically synchronized cells

All inputs are now in place (5 Nov 2026, 22:30).

- screens/: loss-of-function screens for the proposed Functional column (JUMP, PERISCOPE, Funk 2022, MitoCheck, Perturb-seq, BioGRID ORCS), added 11 Nov 2026. See screens/README.md; rebuilt by fetch_screens.py.
- complexportal_9606.tsv: EBI Complex Portal, curated human complexes (members, stoichiometry, expanded participant list), downloaded 12 Nov 2026 from ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv
- uniprot_subunit.tsv: UniProt reviewed human, "Subunit structure" text (rest.uniprot.org stream, fields accession,cc_subunit), 12 Nov 2026
