library(rhdf5)

loom_file <- "out/aucell.loom"
default_metadata <- "{}"  

if (!file.exists(loom_file)) {
    stop("File not found: ", loom_file)
}

fid <- H5Fopen(loom_file, flags = "H5F_ACC_RDWR")
attrs_obj <- H5Gopen(fid, "/attrs")

h5writeAttribute(default_metadata, attrs_obj, "MetaData")

H5Gclose(attrs_obj)
H5Fclose(fid)

attrs <- h5readAttributes(loom_file, "/attrs")
cat("Attributes in '/attrs' after adding MetaData:\n")
print(names(attrs))

