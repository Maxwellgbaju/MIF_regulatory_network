pred_df <- read.csv("MIF_predicted_symbol_allsites.csv")
val_df <- read.csv("MIF_validated_all.csv")

cat("Columns available in predicted data:\n")
print(colnames(pred_df))
cat("\nColumns available in validated data:\n")
print(colnames(val_df))

# Check: does each accession map to exactly one name, or multiple?
if ("mature_mirna_acc" %in% colnames(pred_df)) {
  acc_to_names <- aggregate(mature_mirna_id ~ mature_mirna_acc, data = pred_df, FUN = function(x) paste(unique(x), collapse = " | "))
  multi_name_accs <- acc_to_names[sapply(strsplit(acc_to_names$mature_mirna_id, " \\| "), length) > 1, ]
  cat("\nAccessions mapped to MULTIPLE different name strings (predicted table):\n")
  print(multi_name_accs)
} else {
  cat("\nNo mature_mirna_acc column found in predicted table — cannot cross-check by accession.\n")
}

if ("mature_mirna_acc" %in% colnames(val_df)) {
  acc_to_names_v <- aggregate(mature_mirna_id ~ mature_mirna_acc, data = val_df, FUN = function(x) paste(unique(x), collapse = " | "))
  multi_name_accs_v <- acc_to_names_v[sapply(strsplit(acc_to_names_v$mature_mirna_id, " \\| "), length) > 1, ]
  cat("\nAccessions mapped to MULTIPLE different name strings (validated table):\n")
  print(multi_name_accs_v)
}

# Specifically check your 10 candidates for accession-level presence across all databases,
# re-deriving the agreement table by accession instead of name
targets <- c("hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p", "hsa-miR-23a-5p",
             "hsa-miR-363-5p", "hsa-miR-23b-5p", "hsa-miR-541-3p", "hsa-miR-92a-1-5p",
             "hsa-miR-608", "hsa-miR-1537-3p")

cat("\n\n=== Accession lookup for your 10 candidates ===\n")
for (t in targets) {
  rows <- pred_df[pred_df$mature_mirna_id == t, c("mature_mirna_id", "mature_mirna_acc")]
  rows <- unique(rows)
  cat(t, "->", if(nrow(rows) > 0) paste(rows$mature_mirna_acc, collapse=", ") else "NOT FOUND", "\n")
}
