# Load packages
library(tidyverse)
library(GEOquery)
library(AnnotationDbi)
library(hgu133plus2.db)


# Create folders
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)



#  Load GSE20711


gse <- getGEO(
  "GSE20711",
  GSEMatrix = TRUE,
  destdir = "data/raw"
)

eset <- gse[[1]]




expr <- Biobase::exprs(eset)

meta <- Biobase::pData(eset) %>%
  rownames_to_column("sample")




# Find all metadata columns that contain "subtypege:"
subtype_cols <- names(meta)[
  sapply(
    meta,
    function(x)
      any(
        grepl(
          "^subtypege\\s*:",
          as.character(x),
          ignore.case = TRUE
        )
      )
  )
]


# For each sample, find the cell containing subtypege:
subtype_raw <- apply(
  meta[, subtype_cols, drop = FALSE],
  1,
  function(x) {
    
    match <- x[
      grepl(
        "^subtypege\\s*:",
        x,
        ignore.case = TRUE
      )
    ]
    
    if (length(match) > 0) {
      match[1]
    } else {
      NA
    }
  }
)


# Remove "subtypege:" and keep only subtype name
subtype_data <- tibble(
  sample = meta$sample,
  
  subtype = str_remove(
    subtype_raw,
    regex(
      "^subtypege\\s*:\\s*",
      ignore_case = TRUE
    )
  )
)


#find bcl2 probes

bcl2_probes <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = "BCL2",
  keytype = "SYMBOL",
  columns = "PROBEID"
)

bcl2_probe_ids <- intersect(
  unique(bcl2_probes$PROBEID),
  rownames(expr)
)


#bcl expression

bcl2_expr <- expr[
  bcl2_probe_ids,
  ,
  drop = FALSE
]


# For this initial exploratory analysis:
# average the BCL2 probes for each sample
bcl2_mean <- colMeans(
  bcl2_expr,
  na.rm = TRUE
)




data <- tibble(
  sample = names(bcl2_mean),
  BCL2 = bcl2_mean
) %>%
  
  left_join(
    subtype_data,
    by = "sample"
  ) %>%
  
  filter(
    !is.na(subtype),
    subtype != ""
  )


# Look at the data used for plotting
print(data)


# Save clean data
write_csv(
  data,
  "results/GSE20711_BCL2_subtype.csv"
)


#Plot

p <- ggplot(
  data,
  aes(
    x = subtype,
    y = BCL2
  )
) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.6
  ) +
  
  labs(
    title = "BCL2 expression across breast cancer subtypes",
    x = "Breast cancer subtype",
    y = "BCL2 expression (log2 RMA intensity)"
  ) +
  
  theme_classic()


print(p)


# Save plot
ggsave(
  "results/GSE20711_BCL2_subtype.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)