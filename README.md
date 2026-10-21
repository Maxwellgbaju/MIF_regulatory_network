# MIF_regulatory_network

Analysis code for the manuscript (in preparation):

**Cross-pipeline integrative analysis nominates HIF1A and TFAP2A as candidate regulators of MIF and links MIF to T-cell exhaustion across human cancers**

The study nominates candidate transcriptional and post-transcriptional regulators of macrophage migration inhibitory factor (MIF) from promoter conservation, motif prediction, curated interactions, ChIP-seq occupancy and two microRNA screens; evaluates them in tumour transcriptomes from 32 TCGA cancer types; and tests every expression-based result across independently processed versions of the same TCGA data (original RSEM quantification, PanCanAtlas batch-corrected, Toil re-quantification), after tumour-composition adjustment, by target-based activity inference, against expression-matched genes and with cross-database machine learning. Direct evidence comes from ChIP-Atlas and public perturbation signatures.

A companion spatial transcriptomics study of MIF and immune-cell distribution is available separately (repository `Cancer_atlas_study`).

## Repository structure

```
original_analysis/                      Analyses performed in September 2026
  1_promoter_conservation_and_motifs/   Promoter definition, liftOver, phyloP, JASPAR scans, TRRUST/ENCODE tiers, microRNA scans (R)
  2_mirna_lncrna_network/               microRNA-lncRNA network (multiMiR, ENCORI) - Supplementary Fig. S1
  3_six_cohort_cbioportal_and_ml/       Six-cohort cBioPortal analyses: correlation, purity, immune, CNA/mutation, survival,
                                        four-gene LASSO, cross-database machine learning (cBioPortal arm) - Supp. Figs. S3, S5-S7
  4_hypothesis_free_motif_screen/       Twelve-factor JASPAR screen - Supplementary Fig. S9
  5_cytargetlinker_network/             Independent microRNA screen and lncRNA partners
  6_xena_validation_and_ml/             Tumour vs normal, lncRNA tests, pan-cancer survival screen, Xena machine-learning arm
                                        - Supplementary Figs. S4, S8
revision_2026/                          Analyses performed in October 2026
  analyses/                             mifrev_utils.py (shared helpers) and rev1-rev8 analyses
  tables/                               Tables 1-4 and the immune-marker check
  figures/                              Main Figures 1-8 (run_figures.sh)
  supplementary_figures/                Collection of Supplementary Figures S1-S9; corrected S2 network
results/                                Result tables (CSV) produced by the revision scripts
```

## Revision analyses, in run order

| Script | Purpose |
|---|---|
| `rev1_coexpression.py` | Candidates vs MIF, 32 cancer types (Toil) |
| `rev1b_rela_crosscheck.py` | Correlations across three processing pipelines |
| `rev1c_rela_survival.py` | Survival across three pipelines |
| `rev1d_rela_adjudicate.py` | Patient-level agreement and NF-kB target check |
| `rev2_perturbation.py` | Perturbation signatures (GEO, CREEDS, LINCS) |
| `rev3_tf_activity.py` | Transcription-factor activity (CollecTRI) vs MIF |
| `rev4_mirna_cerna.py` | microRNA-lncRNA expression test (TCGA microRNA-seq) |
| `rev5_chipatlas.py` | ChIP-Atlas binding at the MIF promoter |
| `rev6_hypoxia_mif.py` | Hypoxia programme vs random gene sets |
| `rev7_hallmark_specificity.py` | All 50 Hallmark programmes |
| `rev8_myc_specificity.py` | Expression-matched controls and MYC |
| `table1`–`table3`, `immune_bridge_check.py` | Manuscript Tables 1-4 |
| `run_figures.sh` | Main Figures 1-8 |
| `collect_supp_figures.sh`, `figS2_cytargetlinker_network.py` | Supplementary Figures S1-S9 |

Run revision scripts from the folder containing `mifrev_utils.py`, with a `results/` folder alongside. Analyses download public data on first use and cache it in `cache/`.

## Data sources

All data are public: UCSC Xena (Toil TCGA-TARGET-GTEx, TCGA and PanCanAtlas hubs, including TCGA survival data), cBioPortal (TCGA PanCancer Atlas 2018), JASPAR2024, TRRUST v2, ENCODE, ENCORI/starBase, multiMiR, ChIP-Atlas, CollecTRI/OmniPath, MSigDB Hallmark gene sets and Enrichr perturbation libraries (via GSEApy). Raw data are not redistributed here.

## Software

Python 3 with pandas, NumPy, SciPy, statsmodels, matplotlib, xenaPython, gseapy, lifelines and requests (versions in `requirements.txt`); R/Bioconductor (multiMiR) for the microRNA scans.

## Notes

- The original analyses are provided as run in September 2026; some scripts were superseded by later versions in the same folder (for example `*_v2`, `*_final`). The figures used in the manuscript are produced by the scripts listed in `revision_2026/supplementary_figures/collect_supp_figures.sh`.
- Machine-learning values for Figure 5 are written by `make_ml_values.py` from the reported results of the original machine-learning scripts in `original_analysis/3_...` and `6_...`.

## Contact

Questions about the code can be raised through the Issues tab of this repository.
