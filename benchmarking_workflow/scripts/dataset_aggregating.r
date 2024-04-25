
# Libraries ---------------------------------------------------------------
library(anndata)
library(Seurat)

# Parameters --------------------------------------------------------------
sc_data <- as.character(commandArgs(TRUE)[1])
gene_annotations_path <- as.character(commandArgs(TRUE)[2])
metadata_path <- as.character(commandArgs(TRUE)[3])
nb_embryo <- as.numeric(commandArgs(TRUE)[4])
output_folder <- as.character(commandArgs(TRUE)[5])

dir.create(paste0(output_folder), recursive = T)
options("scipen"=100, "digits"=4)

# Load data ---------------------------------------------------------------
sparse_matrix <- readRDS(sc_data)
gene_annotations <- read.csv(gene_annotations_path, row.names = 1)
rownames(sparse_matrix) <- gene_annotations[rownames(sparse_matrix), "gene_short_name"]
sparse_matrix <- sparse_matrix[-which(duplicated(rownames(sparse_matrix))),]

metadata <- read.csv(metadata_path, row.names = 1)
metadata <- metadata[colnames(sparse_matrix), c(1:5,9:10,19,22)]
metadata$embryo_id <- paste0("embryo_", metadata$embryo_id)
table(metadata$embryo_id)

embryo_list <- names(sort(table(metadata$embryo_id)))

ds_tmp <- embryo_list[1:nb_embryo]
print(ds_tmp)
indices <- metadata$embryo_id %in% ds_tmp

sobj <- CreateSeuratObject(counts = sparse_matrix[,indices], meta.data = metadata[indices,])

saveRDS(sobj, file = paste0(output_folder, "sc_data_", nb_embryo, ".rds"))
