library(multiMiR)

# Try both gene symbol and Entrez ID, with predicted.site="all" (not just conserved)
predicted_symbol <- get_multimir(
  org = "hsa",
  target = "MIF",
  table = "predicted",
  predicted.site = "all",
  summary = TRUE
)
cat("=== By gene symbol 'MIF' ===\n")
print(table(predicted_symbol@data$database))

predicted_entrez <- get_multimir(
  org = "hsa",
  target = "4282",
  table = "predicted",
  predicted.site = "all",
  summary = TRUE
)
cat("\n=== By Entrez ID '4282' ===\n")
print(table(predicted_entrez@data$database))

write.csv(predicted_symbol@data, "MIF_predicted_symbol_allsites.csv", row.names = FALSE)
write.csv(predicted_entrez@data, "MIF_predicted_entrez_allsites.csv", row.names = FALSE)
