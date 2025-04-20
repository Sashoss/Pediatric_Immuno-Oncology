

create_soft_links <- function(input_dir, output_dir, include_dirs = FALSE) {
    if (!dir.exists(input_dir)) stop("Input directory does not exist.")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    all_items <- list.files(input_dir, full.names = TRUE, include.dirs = include_dirs)

    for (item in all_items) {
        item_name <- basename(item)
        target_link <- file.path(output_dir, item_name)
    
        if (!file.exists(target_link)) {
            success <- file.symlink(from = item, to = target_link)
            if (!success) warning(paste("Failed to create symlink for:", item))
        } else {
            message(paste("Link already exists:", target_link))
        }
    }
    message("Symlinks created in: ", output_dir)
}