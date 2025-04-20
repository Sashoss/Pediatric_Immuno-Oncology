
filter_scpca_data <- function(input_dir, output_filepath, cancer_types=c("Glioblastoma"), dtype="filtered") {
    sample_dirs <- list.dirs(input_dir, recursive = TRUE)
    seurat_list <- list()
    metadata_file_path <- file.path(input_dir, "single_cell_metadata.tsv")

    metadata <- read.delim(metadata_file_path, sep="\t", header=TRUE)
    sample_pick = c()
    counter = 0
    for (i in 1:length(sample_dirs)) {
        sample_id <- basename(sample_dirs[i])
        rds_files <- list.files(path = sample_dirs[i], pattern = paste0("\\_", dtype, ".rds$"), full.names = TRUE)
        
        for (file in rds_files) {
            test_sample <- TRUE
            sample_id <- strsplit(sub("\\.[[:alnum:]]+$", "", basename(file)), "_")[[1]][1]
        
            sample_metadata <- metadata[metadata$scpca_library_id == sample_id,]
    
            if (any(sample_metadata$diagnosis %in% cancer_types)) {
                counter = counter + 1
                sce_obj <- readRDS(file)
                seurat_obj <- as.Seurat(sce_obj)
                tryCatch({
                    seurat_obj <- as.Seurat(sce_obj)
                }, error = function(e) {
                    seurat_obj <- as.Seurat(sce_obj, data="counts")
                })
                if(test_sample==FALSE) next
            
                seurat_obj$diagnosis <- sample_metadata$diagnosis
                seurat_obj$age <- sample_metadata$age
                seurat_obj$sex <- sample_metadata$sex

                sample_pick <- c(sample_pick, sample_id)
                seurat_list[[sample_id]] <- seurat_obj
            }
        }
    }

    return(seurat_list)
}


