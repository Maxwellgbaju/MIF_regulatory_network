#!/bin/bash
# Run from /shared/rc/malariarna/MIF_revision_2026 with: conda activate myenv; bash run_figures.sh
# Outputs: figures_final/*.png (600 dpi) and figures_final/*.eps
set -e
python3 -u rev1d_rela_adjudicate.py     # re-runs rev1d so its results are saved as CSV (uses cached data)
python3 -u make_ml_values.py            # writes the ML values from manuscript Tables 12-16 to CSV
python3 -u fig1_study_design.py
python3 -u fig2_mif_tumor_normal.py     # recomputes MIF tumour vs normal; also writes Supplementary Table S1
python3 -u fig3_coexpression.py
python3 -u fig4_pipeline_dependence.py
python3 -u fig5_machine_learning.py
python3 -u fig6_composition_specificity.py
python3 -u fig7_binding.py
python3 -u fig8_mirna_lncrna.py
echo "All figures written to figures_final/"
