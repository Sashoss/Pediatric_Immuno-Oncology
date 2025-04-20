library(Seurat)
library(SeuratDisk)
library(qs)

seurat_obj <- qread("../Step5_Clustering/out/SCPCP000001/annotated_harmony_SCPCP000001_50_2000_3000.qs")
seurat_obj <- SetIdent(seurat_obj, value="cell_label")

SaveLoom(seurat_obj, filename = "./out/init_loom_harmony_SCPCP000001_50_2000_3000.loom", verbose = TRUE)
