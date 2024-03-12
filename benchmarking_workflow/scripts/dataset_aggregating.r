
# Libraries ---------------------------------------------------------------
library(anndata)
library(Seurat)

# Parameters --------------------------------------------------------------
sc_data <- as.character(commandArgs(TRUE)[1])
dataset <- as.numeric(commandArgs(TRUE)[2])
output_folder <- as.character(commandArgs(TRUE)[3])

set.seed(1234)
dir.create(paste0(output_folder), recursive = T)
options("scipen"=100, "digits"=4)
datasets_lists <- c("Lafyatis_Rojas_2019_10Xv1", "Seibold_2020_10Xv2", "Teichmann_Meyer_2019", "Jain_Misharin_2021_10Xv1",
                    "Lafyatis_Rojas_2019_10Xv2", "Seibold_2020_10Xv3", "Jain_Misharin_2021_10Xv2", "Meyer_2019",
                    "Misharin_Budinger_2018" , "Krasnow_2020", "Misharin_2021", "Nawijn_2021", "Barbry_Leroy_2020", "Banovich_Kropski_2020")

ds_tmp <- datasets_lists[1:dataset]
print(ds_tmp)

# Load data ---------------------------------------------------------------
adata <- read_h5ad(sc_data, backed = "r")
print("read done")
adata$var_names <- adata$var$feature_name # We will use gene short name for downstream analyses


#adata_subset <- adata[adata$obs$dataset %in% ds_tmp, ]

indices <- adata$obs$dataset %in% ds_tmp
adata_subset <- AnnData(X = adata[indices]$raw$X,
                        var = adata[indices]$var,
                        obs = adata[indices]$obs)
#This will allow us to construct supervised metacell for each cell type in each sample later in the tutorial
adata_subset$obs$ann <- as.character(adata_subset$obs$ann_level_3)
# For cell without an annotation at the 3rd level we will use the second level of annotation
adata_subset$obs$ann[adata_subset$obs$ann_level_3 == 'None'] = as.character(adata_subset$obs$ann_level_2[adata_subset$obs$ann_level_3 == 'None'])
adata_subset$obs$ann_sample <- paste0(adata_subset$obs$ann,"_",adata_subset$obs$sample)
adata_subset$obs <- adata_subset$obs[,c("sample","dataset","study","ann_level_3","ann","ann_sample")]
print("subsetting done")
write_h5ad(anndata = adata_subset, filename = paste0(output_folder, "sc_data_", dataset, ".h5ad"))