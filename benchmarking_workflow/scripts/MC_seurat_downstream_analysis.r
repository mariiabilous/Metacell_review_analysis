
# Libraries ---------------------------------------------------------------
library(anndata)
library(Seurat)
library(foreach)
library(getopt)

# Parameters --------------------------------------------------------------
color.celltypes <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                     '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                     '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000", '#E4C755', '#F7F398',
                     '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                     '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                     '#968175', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
                     '#476D87', '#E95C59', '#AB3282', '#23452F', '#BD956A', '#585658', '#9FA3A8',
                     '#E0D4CA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863')

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

unintegrated.data <- merge(metacell.objs[[1]], metacell.objs[-1])

unintegrated.data <- NormalizeData(unintegrated.data)
unintegrated.data <- FindVariableFeatures(unintegrated.data)
unintegrated.data <- ScaleData(unintegrated.data)
unintegrated.data <- RunPCA(unintegrated.data)

integrated.data <- IntegrateLayers(
  object = unintegrated.data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

reduc_name <- "harmony"
integrated.data <- FindNeighbors(integrated.data, reduction = reduc_name, dims = 1:30)
integrated.data <- FindClusters(integrated.data, resolution = 2, cluster.name = "clusters")

integrated.data <- RunUMAP(integrated.data, reduction = reduc_name, dims = 1:30, reduction.name = "umap.integrated")
p1 <- DimPlot(
  integrated.data,
  reduction = "umap.integrated",
  group.by = c("embryo_id", "Main_cell_type"),
  combine = FALSE, label.size = 2, cols = color.celltypes
)

# Join layers to perform differential analysis
integrated.data <- JoinLayers(integrated.data)

Idents(integrated.data) <- "Main_cell_type"
markers <- FindAllMarkers(integrated.data)

pdf(paste0(output_path, "datasets_integrated_umap.pdf"), w=15, h=10)
patchwork::wrap_plots(p1)
dev.off()
