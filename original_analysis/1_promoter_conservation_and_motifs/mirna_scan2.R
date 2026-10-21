library(multiMiR)

predicted_full <- get_multimir(
  org = "hsa",
  target = "MIF",
  table = "predicted",
  summary = TRUE
)
pred_df <- predicted_full@data
cat("Predicted-table source breakdown (no cutoff):\n")
print(table(pred_df$database))
write.csv(pred_df, "MIF_predicted_full.csv", row.names = FALSE)
