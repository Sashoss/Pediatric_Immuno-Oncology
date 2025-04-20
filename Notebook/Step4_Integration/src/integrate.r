suppressPackageStartupMessages({
    library(Seurat)
    library(tools) 
    library(harmony)
    library(AnnotationDbi)
    library(org.Hs.eg.db)  
    library(RCurl)
    library(dplyr)
    library(AnnotationHub)
    library(parallel)
    library(future)
    library(qs)
})

source("../Utility/Matrix.utils.R") 
source("../Utility/cell_cycle_scoring.r")

ndim = 50
nfeat = 2000
nmax = 3000

seurat_list <- qread("../Step3_Preprocessing/out/seurat_list_filtered_SCPCP000001.qs")


process_seurat_sample <- function(seurat_sample, sample_name, nfeat, ndim, add_orig_ident=TRUE) {
    DefaultAssay(seurat_sample) <- "RNA"

    if (ncol(seurat_sample) <= 50) {
        message(paste("Skipping", sample_name, "- too few cells"))
        return(NULL)
    }

    if (add_orig_ident==TRUE) {
        seurat_sample@meta.data$orig.ident <- sample_name
    }
    
    seurat_sample <- NormalizeData(seurat_sample, verbose = FALSE)
    seurat_sample <- FindVariableFeatures(
        seurat_sample,
        selection.method = "vst",
        nfeatures = nfeat,
        verbose = FALSE
    )

    seurat_sample <- cell_cycle_scoring(
        seurat_sample, 
        cc_file = "../Utility/Data/CC_Homo_sapiens.csv"
    )

    seurat_sample <- ScaleData(
        seurat_sample,
        vars.to.regress = "CC.Difference", 
        features = rownames(seurat_sample),
        verbose = FALSE
    )

    seurat_sample <- RunPCA(
        seurat_sample, 
        npcs = ndim, 
        features = VariableFeatures(seurat_sample), 
        verbose = FALSE
    )

    return(seurat_sample)
}



seurat_list <- Filter(function(obj) {
    DefaultAssay(obj) <- "RNA"
    nrow(obj) > 50
}, seurat_list)

seurat_list <- Filter(function(obj) {
    DefaultAssay(obj) <- "RNA"
    ncol(obj) > 50
}, seurat_list)



n_cores <- detectCores()
sample_names <- names(seurat_list)

seurat_list_processed <- mclapply(
    sample_names,
    function(name) process_seurat_sample(seurat_list[[name]], name, nfeat, ndim),
    mc.cores = n_cores
)
names(seurat_list_processed) <- sample_names
seurat_list_processed <- Filter(Negate(is.null), seurat_list_processed)

combined <- merge(
    x = seurat_list_processed[[1]],
    y = seurat_list_processed[-1],
    add.cell.ids = names(seurat_list_processed ),
    project = "HarmonyIntegration"
)

combined$filename <- combined$orig.ident
combined$orig.ident <- as.factor(combined$orig.ident)
Idents(combined) <- combined$orig.ident

combined <- JoinLayers(combined)
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(
                        combined, 
                        selection.method = "vst", 
                        nfeatures = nfeat, 
                        verbose = FALSE
            )

combined <- cell_cycle_scoring(combined, cc_file="../Utility/Data/CC_Homo_sapiens.csv")
combined <- ScaleData(combined, 
                    vars.to.regress = "CC.Difference", 
                    features = rownames(combined))

combined <- RunPCA(
        combined, 
        npcs = ndim, 
        features = VariableFeatures(combined), 
        verbose = FALSE
)

combined <- RunHarmony(
    object = combined,
    group.by.vars = "orig.ident",
    assay.use = "RNA",
    reduction.use = "pca",
    dims.use = 1:ndim,
    verbose = TRUE
)



combined <- RunUMAP(combined, reduction = "harmony", dims = 1:ndim)

qsave(combined, file = file.path(paste0("./out/integrated_harmony_SCPCP000001_", 
                                            as.character(ndim), "_", 
                                            as.character(nfeat), "_", 
                                            as.character(nmax), ".qs")
                                )
    )