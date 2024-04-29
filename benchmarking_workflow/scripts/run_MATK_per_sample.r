
# Libraries ---------------------------------------------------------------
library(anndata)
library(Seurat)
library(foreach)
library(getopt)

# Parameters --------------------------------------------------------------
# tool <- as.character(commandArgs(TRUE)[1])
# dataset <- as.numeric(commandArgs(TRUE)[2])
# gamma <- as.numeric(commandArgs(TRUE)[3])
# MATK_params <- as.character(commandArgs(TRUE)[4])
# nb_cores <- as.numeric(commandArgs(TRUE)[5])
# output_path <- as.character(commandArgs(TRUE)[6])

spec = matrix(c(
  'tool',  't', 1, "character", "Metacell construction tool",
  'dataset',  'd', 1, "numeric", "Dataset id",
  'gamma',  'g', 1, "numeric", "Graining level",
  'nb_cores',  'n', 1, "numeric", "Number of cores to parallelize metacells construction",
  'output_path',  'o', 1, "character", "Output folder"
), byrow=TRUE, ncol=5);

opt = getopt(spec)

args <- commandArgs()

tool <- opt$tool
dataset <- opt$dataset
gamma <- opt$gamma
nb_cores <- opt$nb_cores
output_path <- opt$output_path

if(tool == "SuperCell"){
  tool_name <- "SuperCell" 
  MATK_params <- "-n 50 -f 2000 -k 30"
}else if(tool == "SEACells"){
  tool_name <- "SEACells" 
  MATK_params <- "-n 50 -f 2000 -k 30 -e 1e-3"
}else if(tool == "MetaCell"){
  tool_name <- "MetaCell" 
  MATK_params <- ""
}else if(tool == "SEACellsGPU"){
  tool_name <- "SEACells" 
  MATK_params <- "-n 50 -f 2000 -k 30 -u -e 1e-3"
}

# tool <- "SuperCell"
# dataset <- 3
# gamma <- 25
# MATK_params <- "-n 50 -f 2000 -k 30"
# nb_cores <- 5
# output_path <- "results/MC_anno_embryo_id/SuperCell/gamma25/"

cl <- parallel::makeCluster(nb_cores)
doParallel::registerDoParallel(cl)

samples <- 1:dataset

print(tool)
print(samples)

print(MATK_params)

foreach::foreach(sample = 1:length(samples)) %dopar% {
  system(paste0("MATK -t ", tool_name, " -i results/data_split/sc_data/sc_data_noBPcells_", sample, ".rds",
                " -o ", output_path, "/sample", sample, "/ -g ", gamma, " -s seurat ", MATK_params))
}
parallel::stopCluster(cl)

metacell.files <- list.files(paste0(output_path), recursive = T, pattern = "mc_Seurat.rds", full.names = T)

metacell.objs <- lapply(X = metacell.files, function(X){
  # adata <- read_h5ad(X)
  # countMatrix <- Matrix::t(adata$X)
  # colnames(countMatrix) <- adata$obs_names
  # rownames(countMatrix) <- adata$var_names
  # sobj <- Seurat::CreateSeuratObject(counts = countMatrix, meta.data = adata$obs)
  sobj <- readRDS(X)
  sobj <- Seurat::CreateSeuratObject(counts = GetAssayData(sobj, slot = "counts"), meta.data = sobj@meta.data) # to have seurat V5 assays
  sobj <- RenameCells(sobj, add.cell.id = unique(sobj$embryo_id)) # we give unique name to metacells
  return(sobj)
})
merged.mc <- merge(metacell.objs[[1]], metacell.objs[-1])
saveRDS(merged.mc, file = paste0(output_path, "/merged_metacells.rds")) 