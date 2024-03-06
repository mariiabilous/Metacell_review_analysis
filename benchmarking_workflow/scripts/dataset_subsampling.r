
# Libraries ---------------------------------------------------------------
library(anndata)
library(Seurat)

# Parameters --------------------------------------------------------------
sc_data <- as.character(commandArgs(TRUE)[1])
nb_cells <- as.numeric(commandArgs(TRUE)[2])
output_folder <- as.character(commandArgs(TRUE)[3])

set.seed(1234)
dir.create(paste0(output_folder), recursive = T)
options("scipen"=100, "digits"=4)


# Load data ---------------------------------------------------------------
adata <- read_h5ad(sc_data, backed = "r")
adata$var_names <- adata$var$feature_name # We will use gene short name for downstream analyses


indices <- sample(1:nrow(adata), size = nb_cells, replace = F)
adata_subset <- AnnData(X = adata[indices]$raw$X,
                        var = adata[indices]$var,
                        obs = adata[indices]$obs)
#This will allow us to construct supervised metacell for each cell type in each sample later in the tutorial
adata_subset$obs$ann <- as.character(adata_subset$obs$ann_level_3)
# For cell without an annotation at the 3rd level we will use the second level of annotation
adata_subset$obs$ann[adata_subset$obs$ann_level_3 == 'None'] = as.character(adata_subset$obs$ann_level_2[adata_subset$obs$ann_level_3 == 'None'])
adata_subset$obs$ann_sample <- paste0(adata_subset$obs$ann,"_",adata_subset$obs$sample)

write_h5ad(anndata = adata_subset, filename = paste0(output_folder, "sc_data_", nb_cells, "_cells.h5ad"))


