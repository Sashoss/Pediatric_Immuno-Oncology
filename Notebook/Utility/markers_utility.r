create_marker_list <- function(marker_filepath, sheet){
    marker_df <- read.xlsx(marker_filepath, sheet=sheet)
    cell_types <- colnames(marker_df)
    marker_list <- list()
    for (cellname in cell_types) {
        markers <- as.vector(marker_df[, cellname])
        markers <- markers[!is.na(markers)]
        marker_list[[cellname]] <- markers
    }

    filtered_marker_list <- list()
    for (cellname in names(marker_list)) {
        if (length(marker_list[[cellname]]) > 1) {
            filtered_marker_list[[cellname]] <- unique(marker_list[[cellname]])
        }
    }

    return(filtered_marker_list)
}
