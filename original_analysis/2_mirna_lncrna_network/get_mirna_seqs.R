library(mirbase.db)

targets <- c("hsa-miR-1228-5p", "hsa-miR-23a-5p", "hsa-miR-363-5p",
             "hsa-miR-23b-5p", "hsa-miR-1537-3p", "hsa-miR-92a-1-5p")

# mirbase.db keys are typically mature accession IDs; map by name first
all_names <- unlist(as.list(mirbaseMATURE))
head(all_names)  # inspect the format before filtering — don't assume structure
