suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(ggplot2)
  library(R6)
})

# DoubletFinderPipeline R6 class with Plotting options
# Example usage:
# pipeline <- DoubletFinderPipeline$new(seuratObj = your_seurat_obj, min_pc = 10)
# result_seurat_obj <- pipeline$runAnalysis(output_dir = "output/plots", generate_extra_plots = TRUE)

DoubletFinderPipeline <- R6Class("DoubletFinderPipeline",
    public = list(
                seuratObj = NULL,
                min_pc = NULL,
                doublet_rate = NULL,
                pN = NULL,
                res = NULL,
                n_cores = NULL,
                optimal_pk = NULL,
                homotypic_prop = NULL,
                nExp_poi = NULL,
                nExp_poi_adj = NULL,
    
    # Constructor: initializes the pipeline with a Seurat object and parameters.
    initialize = function(seuratObj, min_pc, doublet_rate = 0.075, pN = 0.25, res = 0.05, n_cores = 1) {
      self$seuratObj <- seuratObj
      self$min_pc <- min_pc
      self$doublet_rate <- doublet_rate
      self$pN <- pN
      self$res <- res
      self$n_cores <- n_cores
    },
    
    # Preprocess the data: normalization, variable feature selection, scaling, PCA, UMAP, neighbor finding, and clustering.
    preprocessData = function() {
      message("Preprocessing data...")
      self$seuratObj <- NormalizeData(self$seuratObj)
      self$seuratObj <- FindVariableFeatures(self$seuratObj)
      self$seuratObj <- ScaleData(self$seuratObj)
      self$seuratObj <- RunPCA(self$seuratObj, npcs = self$min_pc, features = VariableFeatures(self$seuratObj))
      self$seuratObj <- RunUMAP(self$seuratObj, dims = 1:self$min_pc)
      self$seuratObj <- FindNeighbors(self$seuratObj, dims = 1:self$min_pc)
      self$seuratObj <- FindClusters(self$seuratObj, resolution = self$res)
    },
    
    # Run a parameter sweep to determine the optimal pK value.
    runParameterSweep = function() {
      message("Running parameter sweep for optimal pK...")
      sweep.list <- paramSweep(self$seuratObj, PCs = 1:self$min_pc, num.cores = self$n_cores)
      sweep.stats <- summarizeSweep(sweep.list)
      bcmvn <- find.pK(sweep.stats)
      bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric), ]
      self$optimal_pk <- as.numeric(as.character(bcmvn_max$pK))
      message("Optimal pK determined: ", self$optimal_pk)
    },
    
    # Run the DoubletFinder algorithm using the optimal parameters.
    runDoubletDetection = function() {
      message("Estimating homotypic doublet proportion...")
      annotations <- self$seuratObj@meta.data$seurat_clusters
      self$homotypic_prop <- modelHomotypic(annotations)
      
      message("Calculating expected number of doublets...")
      self$nExp_poi <- round(self$doublet_rate * nrow(self$seuratObj@meta.data))
      self$nExp_poi_adj <- max(0, round(self$nExp_poi * (1 - self$homotypic_prop)))
      
      message("Running DoubletFinder with optimal parameters...")
      self$seuratObj <- doubletFinder(
        seu = self$seuratObj,
        PCs = 1:self$min_pc,
        pN = self$pN,
        pK = self$optimal_pk,
        nExp = self$nExp_poi_adj
      )
      
      # Rename meta data columns for clarity
      self$seuratObj@meta.data <- self$seuratObj@meta.data %>%
        rename_with(~ paste0("DFC_", seq_along(.)), matches("^DF.classifications")) %>%
        rename_with(~ paste0("pANN_", seq_along(.)), matches("^pANN_"))
    },
    
    # Generate default publication-quality plots including UMAP, pK selection plot, doublet counts, and pANN histogram.
    generatePlots = function(output_dir = "./plots") {
      message("Generating default plots...")
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      # 1. UMAP plot colored by doublet classification (assumes first DFC_ column)
      doublet_class_col <- grep("^DFC_", colnames(self$seuratObj@meta.data), value = TRUE)[1]
      p1 <- DimPlot(self$seuratObj, reduction = "umap", group.by = doublet_class_col, label = TRUE) +
        ggtitle("UMAP by Doublet Classification")
      ggsave(filename = file.path(output_dir, "umap_doublet_classification.png"), plot = p1)
      
      # 2. pK selection plot using the results from a parameter sweep
      sweep.list <- paramSweep(self$seuratObj, PCs = 1:self$min_pc, num.cores = self$n_cores)
      sweep.stats <- summarizeSweep(sweep.list)
      bcmvn <- find.pK(sweep.stats)
      p2 <- ggplot(bcmvn, aes(x = as.numeric(as.character(pK)), y = BCmetric)) +
        geom_point() +
        geom_line() +
        ggtitle("pK Selection Plot") +
        xlab("pK") +
        ylab("BCmetric")
      ggsave(filename = file.path(output_dir, "pK_selection.png"), plot = p2)
      
      # 3. Bar plot of doublet counts
      doublet_counts <- table(self$seuratObj@meta.data[[doublet_class_col]])
      df_counts <- as.data.frame(doublet_counts)
      p3 <- ggplot(df_counts, aes(x = Var1, y = Freq)) +
        geom_bar(stat = "identity") +
        ggtitle("Doublet Counts") +
        xlab("Doublet Classification") +
        ylab("Count")
      ggsave(filename = file.path(output_dir, "doublet_counts.png"), plot = p3)
      
      # 4. Histogram of pANN distribution (assumes first pANN_ column)
      pANN_col <- grep("^pANN_", colnames(self$seuratObj@meta.data), value = TRUE)[1]
      if (!is.null(pANN_col) && nchar(pANN_col) > 0) {
        p4 <- ggplot(self$seuratObj@meta.data, aes_string(x = pANN_col)) +
          geom_histogram(bins = 30) +
          ggtitle("Distribution of pANN Values") +
          xlab("pANN") +
          ylab("Frequency")
        ggsave(filename = file.path(output_dir, "pANN_distribution.png"), plot = p4)
      }
      
      message("Default plots saved to: ", output_dir)
    },
    
    # Generate Violin plots for a set of features. Default features include nFeature_RNA, nCount_RNA, and percent.mt.
    generateViolinPlots = function(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), output_dir = "./plots") {
      message("Generating Violin plots...")
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      for(feature in features) {
        p <- VlnPlot(self$seuratObj, features = feature, pt.size = 0.1) +
          ggtitle(paste("Violin Plot for", feature))
        ggsave(filename = file.path(output_dir, paste0("violin_", feature, ".png")), plot = p)
      }
      message("Violin plots saved to: ", output_dir)
    },
    
    # Generate a custom DimPlot with user-defined options.
    generateCustomDimPlot = function(reduction = "umap", group_by = "seurat_clusters", label = TRUE, output_file = "custom_dimplot.png", output_dir = "./plots") {
      message("Generating custom DimPlot...")
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      p <- DimPlot(self$seuratObj, reduction = reduction, group.by = group_by, label = label) +
           ggtitle(paste("DimPlot (", reduction, ") colored by", group_by))
      ggsave(filename = file.path(output_dir, output_file), plot = p)
      message("Custom DimPlot saved to: ", file.path(output_dir, output_file))
    },
    
    # Run the full analysis: preprocessing, parameter sweep, doublet detection, and plot generation.
    # Optionally, generate additional Violin plots and custom DimPlots.
    runAnalysis = function(output_dir = "./plots", generate_extra_plots = TRUE) {
      self$preprocessData()
      self$runParameterSweep()
      self$runDoubletDetection()
      self$generatePlots(output_dir)
      if (generate_extra_plots) {
        self$generateViolinPlots(output_dir = output_dir)
        # Example of generating a custom DimPlot for seurat_clusters:
        self$generateCustomDimPlot(group_by = "seurat_clusters", output_file = "custom_dimplot_seurat_clusters.png", output_dir = output_dir)
      }
      return(self$seuratObj)
    }
  )
)


