# Libraries ---------------------------------------------------------------
library(Seurat)
library(BPCells)


# Parameters --------------------------------------------------------------
nb_embryo <- as.numeric(commandArgs(TRUE)[1])
input_type <- as.character(commandArgs(TRUE)[2])
data_path <- as.character(commandArgs(TRUE)[3])
output_path <- as.character(commandArgs(TRUE)[4])

color.celltypes <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                     '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                     '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000", '#E4C755', '#F7F398',
                     '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                     '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                     '#968175', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
                     '#476D87', '#E95C59', '#AB3282', '#23452F', '#BD956A', '#585658', '#9FA3A8',
                     '#E0D4CA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863')
                       
# nb_embryo <- 4
# input_type <- "sc"
# data_path <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/Metacell_review/revisions/Metacell_review_analysis/benchmarking_workflow/results/by_dataset/"

if(input_type == "sc_BPCells"){
  seurat_pattern <- "/sc_data_BPcells_"
}else{
  seurat_pattern <- "/sc_data_noBPcells_"
}

# Load data ---------------------------------------------------------------

seurat.files <- sapply(1:nb_embryo, FUN = function(x){paste0(data_path, seurat_pattern, x, ".rds")}) 

seurat.objs <- lapply(X = 1:length(seurat.files), function(X){
  sobj <- readRDS(seurat.files[X])

  if(grepl("mc_data", seurat_pattern)){
    sobj <- RenameCells(sobj, add.cell.id = unique(sobj$embryo_id)) # we give unique name to metacells
  }
  sobj$dataset <- names(seurat.files)[X]
  return(sobj)
})

# Merge seurat objects
unintegrated.data <- merge(seurat.objs[[1]], seurat.objs[-1])

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

pdf(paste0(output_path, input_type, "_", nb_embryo, "datasets_integrated_umap.pdf"), w=15, h=10)
patchwork::wrap_plots(p1)
dev.off()

# saveRDS(integrated.data, file = paste0(data_path, input_type, "_", nb_embryo, "datasets_integrated_data.rds"))