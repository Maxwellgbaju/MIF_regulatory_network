library(multiMiR)

predicted_entrez <- get_multimir(
  org = "hsa",
  target = "4282",
  target.type = "entrez",
  table = "predicted",
  predicted.site = "all",
  summary = TRUE
)
pred_df <- predicted_entrez@data
cat("Predicted-table source breakdown (Entrez ID, all sites):\n")
print(table(pred_df$database))
write.csv(pred_df, "MIF_predicted_entrez.csv", row.names = FALSE)
