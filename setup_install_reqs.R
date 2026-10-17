
install.packages(c(
  "tidyverse",
  "rstatix"
))


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GEOquery",
  "Biobase",
  "AnnotationDbi",
  "hgu133plus2.db"
))