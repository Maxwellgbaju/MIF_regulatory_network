"""Writes the cross-database machine-learning results reported in the manuscript (Tables 12-16, September analysis)
to CSV so Figure 5 can be drawn from files. If the original ML output files are available, replace these CSVs with them."""
import os, pandas as pd
os.makedirs("results", exist_ok=True)
stab = pd.DataFrame([
    # gene, cBioPortal stability %, cBioPortal coefficient, Xena stability %, Xena coefficient   (Tables 12 and 13)
    ("TP53", 100.0, 0.147, 66.0, 0.021), ("IRF1", 100.0, 0.107, 83.0, 0.047), ("RELA", 100.0, 0.135, 7.0, 0.0),
    ("KLF4", 99.0, -0.199, 100.0, -0.116), ("NR1H2", 98.0, 0.113, 8.0, 0.0), ("TFAP2A", 97.0, 0.044, 100.0, 0.222),
    ("NORAD", 89.0, -0.059, 100.0, -0.194), ("EGR1", 86.0, -0.057, 14.0, 0.0), ("ESR1", 85.0, 0.054, 53.0, 0.0),
    ("NFIL3", 36.0, 0.0, 2.0, 0.0), ("MIF", 16.0, 0.0, 17.0, 0.0), ("RXRA", 16.0, 0.0, 91.0, -0.068),
    ("SP1", 14.0, 0.0, 0.0, 0.0), ("IRF2", 13.0, 0.0, 16.0, 0.0), ("AL356488.2", None, None, 64.0, 0.013)],
    columns=["gene", "stability_cbioportal", "coef_cbioportal", "stability_xena", "coef_xena"])
stab.to_csv("results/ml_stability_selection.csv", index=False)
cidx = pd.DataFrame([  # Table 14
    ("cBioPortal (n=510)", "Multi-gene LASSO", 0.614, 0.047), ("cBioPortal (n=510)", "RELA alone", 0.590, 0.046),
    ("cBioPortal (n=510)", "MIF alone", 0.572, 0.047), ("cBioPortal (n=510)", "Random survival forest", 0.663, 0.045),
    ("Xena (n=530)", "Multi-gene LASSO", 0.643, 0.044), ("Xena (n=530)", "RELA alone", 0.540, 0.044),
    ("Xena (n=530)", "MIF alone", 0.545, 0.046)], columns=["dataset", "model", "c_index", "sd"])
cidx.to_csv("results/ml_cindex.csv", index=False)
rsf = pd.DataFrame([  # Table 16 (cBioPortal, 15 repeats)
    ("TFAP2A", .0810, .0108), ("TP53", .0266, .0055), ("RXRA", .0227, .0024), ("RELA", .0206, .0054), ("NORAD", .0142, .0031),
    ("NR1H2", .0137, .0021), ("IRF1", .0134, .0024), ("KLF4", .0131, .0038), ("EGR1", .0117, .0023), ("ESR1", .0106, .0010),
    ("IRF2", .0096, .0015), ("NFIL3", .0086, .0008), ("SP1", .0072, .0007), ("MIF", .0061, .0013)],
    columns=["gene", "importance", "sd"])
rsf.to_csv("results/ml_rsf_importance.csv", index=False)
print("wrote results/ml_stability_selection.csv, ml_cindex.csv, ml_rsf_importance.csv")
