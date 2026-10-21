library(multiMiR)
val_df <- read.csv("MIF_validated_all.csv")

mtb <- unique(val_df$mature_mirna_id[val_df$database == "mirtarbase"])
tb  <- unique(val_df$mature_mirna_id[val_df$database == "tarbase"])
mrec <- unique(val_df$mature_mirna_id[val_df$database == "mirecords"])

overlap <- intersect(mtb, tb)
cat("miRTarBase hits:", length(mtb), "\n")
cat("TarBase hits:", length(tb), "\n")
cat("Overlap (validated by both):", length(overlap), "\n\n")
print(overlap)

writeLines(overlap, "MIF_miRNA_validated_overlap.txt")
