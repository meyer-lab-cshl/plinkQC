df <- read.table(
  "/Users/meyerlab/RpackageDevelopment/plinkQC/extdata/data.bim",
  header = FALSE,
  sep = "\t",
  dec = "."
)

ids <- paste(df$V1, df$V4, df$V6, df$V5, sep = ":")
df$V2 <- ids
head(df)


write.table(
  df,
  file = "/Users/meyerlab/RpackageDevelopment/plinkQC/extdata/dataIDs.bim",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
