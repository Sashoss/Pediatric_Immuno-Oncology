#!/usr/bin/env Rscript

# process_seurat_doubletfinder.R
# Description: Processes a single Seurat object with DoubletFinder and saves the output.
# Usage: Rscript process_seurat_doubletfinder.R --input <input_file.rds> --output_dir <output_directory>

# -----------------------------
# 1. Load Required Libraries
# -----------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(optparse)
})

# -----------------------------
# 2. Parse Command-Line Arguments
# -----------------------------

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to input .rds Seurat object", metavar = "FILE"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Path to output directory", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Both input file and output directory must be specified.", call. = FALSE)
}

input_file <- opt$input
output_dir <- opt$output_dir

# -----------------------------
# 3. Define the DoubletFinder Function
# -----------------------------

callDoubletFinder <- function(seuratObj) {
  # Pre-process Seurat object
  seuratObj <- NormalizeData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj)
  seuratObj <- ScaleData(seuratObj)
  seuratObj <- RunPCA(seuratObj)
  
  # Determine number of PCs to use
  stdv <- seuratObj[["pca"]]@stdev
  sum_stdv <- sum(stdv)
  percent_stdv <- (stdv / sum_stdv) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(which((percent_stdv[1:(length(percent_stdv) - 1)] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)
  
  # Continue pre-processing
  seuratObj <- RunUMAP(seuratObj, dims = 1:min_pc)
  seuratObj <- FindNeighbors(seuratObj, dims = 1:min_pc)              
  seuratObj <- FindClusters(seuratObj, resolution = 0.1)
  
  # DoubletFinder parameter sweep
  sweep.list <- paramSweep(seuratObj, PCs = 1:min_pc, num.cores = 20) # Use 1 core to manage memory
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK selection
  bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal_pk <- as.numeric(as.character(bcmvn_max$pK))
  
  # Estimate homotypic doublet proportion
  annotations <- seuratObj@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations) 
  nExp_poi <- round(0.075 * nrow(seuratObj@meta.data)) # Adjust 7.5% as needed
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  
  # Run DoubletFinder
  seuratObj <- doubletFinder(seu = seuratObj, 
                             PCs = 1:min_pc, 
                             pN = 0.25,  # pN is typically set to 0.25
                             pK = optimal_pk,
                             nExp = nExp_poi_adj)
  
  # Safely rename columns to avoid duplication
  seuratObj@meta.data <- seuratObj@meta.data %>%
    # Rename DF.classifications columns with unique suffix
    rename_with(~ paste0("DFC_", seq_along(.)), matches("^DF.classifications")) %>%
    # Rename pANN columns with unique suffix
    rename_with(~ paste0("pANN_", seq_along(.)), matches("^pANN_"))
  
  return(seuratObj)
}

# -----------------------------
# 4. Process the Seurat Object
# -----------------------------

# Read the Seurat object
cat("Reading Seurat object from:", input_file, "\n")
seurat_obj <- readRDS(input_file)

# Apply DoubletFinder
cat("Running DoubletFinder on:", input_file, "\n")
seurat_obj <- callDoubletFinder(seurat_obj)

# -----------------------------
# 5. Save the Processed Object
# -----------------------------

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define output filename
input_basename <- tools::file_path_sans_ext(basename(input_file))
output_file <- file.path(output_dir, paste0(input_basename, "_doubletfinder.rds"))

# Save the Seurat object
cat("Saving processed Seurat object to:", output_file, "\n")
saveRDS(seurat_obj, file = output_file)

# -----------------------------
# 6. Clean Up
# -----------------------------

rm(seurat_obj)
gc()

cat("Processing complete for:", input_file, "\n")
