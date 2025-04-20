source("../Utility/doublets.r")
source("../Utility/seurat_tools.r")

library(parallel)
library(ggplot2)

seurat_list <- readRDS("../Step2_Generate_Seurat_Object/out/seurat_list.rds")

for (sample_name in names(seurat_list)) {
    clean_seurat_object <- clean_seurat(seurat_list[[sample_name]])
    seurat_list[[sample_name]] <- suppressMessages(
                                                suppressWarnings(
                                                    callDoubletFinder(
                                                                clean_seurat_object, 
                                                                min_pc=30, 
                                                                doublet_rate=0.075, 
                                                                pN=0.25, 
                                                                res=0.5, 
                                                                n_cores=20
                                                    )
                                            )
    )
} 


saveRDS(seurat_list, file = "./out/seurat_list_doublets.rds")

