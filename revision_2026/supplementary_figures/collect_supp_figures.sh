#!/bin/bash
# Collects supplementary figures from the original September scripts into supp_figures/ (PNG + EPS).
# Reuses existing outputs; reruns the original script in its own folder only if PNG or EPS is missing.
BASE=/shared/rc/malariarna
OUT=$BASE/MIF_revision_2026/supp_figures
mkdir -p "$OUT"
get () {   # folder  script  output_basename  new_name
  local dir="$BASE/$1" script="$2" stem="$3" name="$4"
  if [ ! -f "$dir/$stem.png" ] || [ ! -f "$dir/$stem.eps" ]; then
    echo ">> $name: rerunning $1/$script"
    ( cd "$dir" && python3 -u "$script" ) 2>&1 | tail -3
  fi
  for ext in png eps; do
    if [ -f "$dir/$stem.$ext" ]; then cp "$dir/$stem.$ext" "$OUT/$name.$ext"; echo "   $name.$ext"
    else echo "   MISSING: $dir/$stem.$ext"; fi
  done
}
get MIF_lncRNA_analysis               build_network_final.py       MIF_regulatory_network_v2            S1_candidate_A_network
get MIF_promoter_comparative_analysis network_v3_final.py          MIF_network_v3_final                 S1_candidate_B_network
get MIF_student_network_v2            build_student_network_v2.py  student_MIF_network_v2               S2_cytargetlinker_network
get MIF_Phase3_pancancer              phase3_4_figure_final.py     Phase3_4_figure_final                S3_mutation_cna
get MIF_student_network_v2/gepia_test lncrna_combined_figure.py    lncRNA_combined_figure               S4_norad_al356488
get MIF_Phase3_pancancer              make_final_figure.py         Phase3_MIF_correlation_figure_final  S5_initial_correlation
get MIF_Phase3_pancancer              immune_figures.py            Immune_partial_heatmaps              S6a_immune_heatmaps
get MIF_Phase3_pancancer              pdcd1_slopeplot_v2.py        PDCD1_slopeplot_v2                   S6b_pdcd1_shift
get MIF_Phase3_pancancer              phase4_2_final_figure.py     phase4_2_final_figure                S7_initial_lasso
get MIF_student_network_v2/gepia_test survival_forest_plot.py      survival_forest_plot                 S8a_survival_screen
get MIF_student_network_v2/gepia_test mif_luad_km.py               MIF_LUAD_KM_xena                     S8b_mif_luad_km
get MIF_PROJECT_PY_REBECCA            generate_mif_figure.py       mif_tf_jaspar_figure                 S9_motif_screen
echo; echo "Collected in $OUT:"; ls "$OUT"
