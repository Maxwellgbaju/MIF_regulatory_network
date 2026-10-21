library(multiMiR)

# Predicted targets from TargetScan + miRDB (remote query, live database)
predicted <- get_multimir(
  org = "hsa",
  target = "MIF",
  table = "predicted",
  predicted.cutoff.type = "p",
  predicted.cutoff = 20,      # top 20th percentile — adjust for stringency
  summary = TRUE
)
pred_df <- predicted@data
cat("Predicted-table source breakdown:\n")
print(table(pred_df$database))
write.csv(pred_df, "MIF_predicted_all.csv", row.names = FALSE)

# Validated targets from miRTarBase (+ TarBase, miRecords)
validated <- get_multimir(
  org = "hsa",
  target = "MIF",
  table = "validated",
  summary = TRUE
)
val_df <- validated@data
cat("\nValidated-table source breakdown:\n")
print(table(val_df$database))
write.csv(val_df, "MIF_validated_all.csv", row.names = FALSE)

# --- Cross-database agreement: require the miRNA in >=2 of {targetscan, mirdb, mirtarbase} ---
ts  <- unique(pred_df$mature_mirna_id[pred_df$database == "targetscan"])
mdb <- unique(pred_df$mature_mirna_id[pred_df$database == "mirdb"])
mtb <- unique(val_df$mature_mirna_id[val_df$database == "mirtarbase"])

all_mirnas <- unique(c(ts, mdb, mtb))
agreement <- data.frame(
  mirna = all_mirnas,
  in_targetscan = all_mirnas %in% ts,
  in_mirdb = all_mirnas %in% mdb,
  in_mirtarbase = all_mirnas %in% mtb
)
agreement$n_sources <- rowSums(agreement[,2:4])
agreement <- agreement[order(-agreement$n_sources), ]

write.csv(agreement, "MIF_miRNA_agreement.csv", row.names = FALSE)

cat("\nmiRNAs supported by >=2 of 3 databases:\n")
print(agreement[agreement$n_sources >= 2, ])
