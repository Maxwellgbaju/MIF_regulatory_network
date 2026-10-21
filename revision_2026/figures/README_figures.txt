MIF manuscript figure scripts. Copy all files into /shared/rc/malariarna/MIF_revision_2026 (same folder as mifrev_utils.py).

Script                          Input files (in results/)                                   Output (figures_final/)
fig1_study_design.py            none (schematic)                                            Figure1_study_design
fig2_mif_tumor_normal.py        Xena Toil data (fetched); writes fig2_MIF_tumor_vs_normal.csv Figure2_MIF_tumor_vs_normal
fig3_coexpression.py            rev1_MIF_coexpression_all.csv                               Figure3_coexpression
fig4_pipeline_dependence.py     rev1b_rela_crosscheck.csv, rev1c_survival_by_pipeline.csv,
                                rev1d_agreement.csv, rev1d_nfkb_targets.csv                 Figure4_pipeline_dependence
fig5_machine_learning.py        ml_stability_selection.csv, ml_cindex.csv, ml_rsf_importance.csv  Figure5_machine_learning
fig6_composition_specificity.py rev3_tf_activity_vs_MIF.csv, rev6_hypoxia_vs_MIF.csv,
                                rev7_hallmark_specificity_summary.csv, rev8_program_specificity.csv  Figure6_composition_specificity
fig7_binding.py                 rev5_chipatlas_MIF_promoter.csv                             Figure7_binding
fig8_mirna_lncrna.py            rev4_mirna_cerna_summary.csv                                Figure8_mirna_lncrna

Notes
- rev1d_rela_adjudicate.py replaces the earlier rev1d; it now saves its results as CSV (needed by Figure 4).
- make_ml_values.py writes the ML results reported in manuscript Tables 12-16 (September analysis). If the original
  ML output files are located, replace the three ml_*.csv files with values from those files.
- fig2 recomputes the tumour-vs-normal comparison; if the number of cancer types differs from the 23 in the
  previous manuscript, report what this script produces and update the text.
- figstyle.py sets fonts (Arial/Helvetica, falling back to DejaVu Sans), 600-dpi PNG and EPS output, and lowercase panel labels.
