library(multiMiR)

pred_df <- read.csv("MIF_predicted_symbol_allsites.csv")
val_df <- read.csv("MIF_validated_all.csv")

ts   <- unique(pred_df$mature_mirna_id[pred_df$database == "targetscan"])
pita <- unique(pred_df$mature_mirna_id[pred_df$database == "pita"])
mir  <- unique(pred_df$mature_mirna_id[pred_df$database == "miranda"])
mc   <- unique(pred_df$mature_mirna_id[pred_df$database == "microcosm"])
mtb  <- unique(val_df$mature_mirna_id[val_df$database == "mirtarbase"])
tb   <- unique(val_df$mature_mirna_id[val_df$database == "tarbase"])

all_mirnas <- unique(c(ts, pita, mir, mc, mtb, tb))
tbl <- data.frame(
  mirna = all_mirnas,
  targetscan = all_mirnas %in% ts,
  pita = all_mirnas %in% pita,
  miranda = all_mirnas %in% mir,
  microcosm = all_mirnas %in% mc,
  mirtarbase = all_mirnas %in% mtb,
  tarbase = all_mirnas %in% tb
)
tbl$n_predicted <- rowSums(tbl[,c("pita","miranda","microcosm")])
tbl$n_validated <- rowSums(tbl[,c("mirtarbase","tarbase")])

# Final candidate rule: in TargetScan AND (>=1 other predicted DB OR >=1 validated DB)
final <- tbl[tbl$targetscan & (tbl$n_predicted >= 1 | tbl$n_validated >= 1), ]
final <- final[order(-final$n_validated, -final$n_predicted), ]

write.csv(tbl, "MIF_miRNA_full_agreement.csv", row.names = FALSE)
write.csv(final, "MIF_miRNA_final_candidates.csv", row.names = FALSE)

cat("Total unique miRNAs across all sources:", nrow(tbl), "\n")
cat("Final candidates (TargetScan + corroboration):", nrow(final), "\n\n")
print(final)
