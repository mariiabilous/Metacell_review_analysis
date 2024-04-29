library(ggplot2)
library(scales)
library(lubridate)

path_to_benchmarking_res <- "./"

sc_data <- paste0(path_to_benchmarking_res, "results/data_full/sc_data.rds")
gene_annotations_path <- paste0(path_to_benchmarking_res, "results/data_full/gene.csv")
metadata_path <- paste0(path_to_benchmarking_res, "results/data_full/metadata.csv")

sparse_matrix <- readRDS(sc_data)
gene_annotations <- read.csv(gene_annotations_path, row.names = 1)
rownames(sparse_matrix) <- gene_annotations[rownames(sparse_matrix), "gene_short_name"]
sparse_matrix <- sparse_matrix[-which(duplicated(rownames(sparse_matrix))),]

metadata <- read.csv(metadata_path, row.names = 1)
metadata <- metadata[colnames(sparse_matrix), c(1:5,9:10,19,22)]
metadata$embryo_id <- paste0("embryo_", metadata$embryo_id)
table(metadata$embryo_id)

nb_cells_all <- cumsum(sort(table(metadata$embryo_id)))
nb_cells_all <- nb_cells_all[c(3,7,10,12,15,20,25,30,35,40,50)]
names(nb_cells_all) <- as.character(c(3,7,10,12,15,20,25,30,35,40,50))

# Using time --------------------------------------------------------------
data_folder <- paste0(path_to_benchmarking_res, "results/benchmarks_time/")

prof_files = list.files(data_folder, pattern = "gamma")

prof_summary = data.frame(gamma = vector(),
                          nb_cells = vector(),
                          method = vector(),
                          MaxRSS = vector(),
                          cpu_time = vector(),
                          elapsed_time = vector())
for(i in prof_files){
  print(i)
  gamma <- as.numeric(gsub("gamma","",strsplit(i, "_")[[1]][1]))
  dataset <- strsplit(i, "_")[[1]][2] 
  method <- gsub(".txt","",paste(strsplit(i, "_")[[1]][-c(1:2)], collapse = "_"))
  nb_cells <- nb_cells_all[dataset]
  
  if(!file.info(paste0(data_folder, i))$size == 0){
    time_res = read.table(paste0(data_folder, i), header = F, sep = "\t")$V2
    
    MaxRSS = time_res[10]
    MaxRSS = as.numeric(unlist(strsplit(MaxRSS,": "))[[2]])/1000000
    
    e_time = unlist(strsplit(time_res[5],": "))[[2]]
    if(length(unlist(strsplit(e_time,":")))==2){
      e_time <- paste0("0:",e_time)
    }
    running_time <- hms(e_time) # format to 'hours:minutes:seconds'
    running_time <- hour(running_time)*60*60 + minute(running_time)*60 + second(running_time)
    running_time <- running_time / 60
    
    prof_summary = rbind(prof_summary,data.frame(gamma = gamma,
                                                 nb_cells = nb_cells,
                                                 method = method,
                                                 MaxRSS = MaxRSS,
                                                 cpu_time = as.numeric(unlist(strsplit(time_res[2],": "))[[2]])/60,
                                                 elapsed_time = running_time))
  }
  
}


pdf("time_profiling_final.pdf", w=5.5, h=4)


prof_summary$method <- gsub("MetaCell", "MC2",prof_summary$method)

colors <- c("MC2_multithreads" = "#FF9900",
            "MC2" = "#FF6600", 
            "SEACells_gpu" = "#FF99CC",
            "SEACells" = "#FF3399",
            "SuperCell_approx" = "#99CCFF",
            "SuperCell" = "#3399FF",
            "SuperCell_extrapolation" = "#3399FF")
shapes <- c("MC2_multithreads" = 15,
            "MC2" = 15, 
            "SEACells_gpu" = 16,
            "SEACells" = 16,
            "SuperCell_approx" = 17,
            "SuperCell" = 17,
            "SuperCell_extrapolation" = 17)

p1 <- ggplot(prof_summary, aes(x=nb_cells, y=cpu_time, colour=method, group = method, shape = method)) + 
  geom_line() + 
  geom_point(size = 2) + theme_bw() + 
  xlab("Number of cells") + ylab("Time (minutes)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colors) + scale_shape_manual(values = shapes)+
  theme(legend.title=element_blank())

print(p1 + ggtitle("Metacell construction (Time)")) 

p2 <- ggplot(prof_summary, aes(x=nb_cells, y=MaxRSS, colour=method, group = method, shape = method)) + 
  geom_line() + 
  geom_point(size = 2) + xlab("Number of cells") +  ylab("Max RSS (Gb)") + theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colors) + scale_shape_manual(values = shapes)+
  theme(legend.title=element_blank())
print(p2 + ggtitle("Metacell construction (Memory)"))


dev.off()


# With downstream analyses -----------------------------------------------------
path_to_downstream_res <- "./"
data_folder <- paste0(path_to_downstream_res, "results/MC_anno_embryo_id/benchmarks_time/")

prof_files = list.files(data_folder)

prof_summary = data.frame(nb_cells = vector(),
                          method = vector(),
                          MaxRSS = vector(),
                          cpu_time = vector(),
                          elapsed_time = vector())
for(i in prof_files){
  print(i)
  if(grepl("downstream", i)){
    dataset <- gsub("sc_|sc_BPCells_|datasets_downstream.txt|MC_gamma75_MetaCell_|MC_gamma75_SEACells_|MC_gamma75_SuperCell_", "", i)
    method <- ifelse(grepl("sc_", i), "sc_downstream", "MC+downstream")
    method <- ifelse(grepl("BPCells", i), "scBPCells_downstream", method)
    method <- ifelse(grepl("MC_", i), paste0(strsplit(i, "_")[[1]][3], "+downstream"), method)
  }else{
    dataset <- strsplit(i, "_")[[1]][2] # as.numeric
    method <- gsub(".txt","",paste(strsplit(i, "_")[[1]][-c(1:2)], collapse = "_"))
  }
  
  nb_cells <- nb_cells_all[dataset]
  
  if(!file.info(paste0(data_folder, i))$size == 0){
    time_res = read.table(paste0(data_folder, i), header = F, sep = "\t")$V2
    
    MaxRSS = time_res[10]
    MaxRSS = as.numeric(unlist(strsplit(MaxRSS,": "))[[2]])/1000000
    
    e_time = unlist(strsplit(time_res[5],": "))[[2]]
    if(length(unlist(strsplit(e_time,":")))==2){
      e_time <- paste0("0:",e_time)
    }
    running_time <- hms(e_time) # format to 'hours:minutes:seconds'
    running_time <- hour(running_time)*60*60 + minute(running_time)*60 + second(running_time)
    running_time <- running_time / 60
    
    
    prof_summary = rbind(prof_summary,data.frame(nb_cells = nb_cells,
                                                 method = method,
                                                 MaxRSS = MaxRSS,
                                                 cpu_time = as.numeric(unlist(strsplit(time_res[2],": "))[[2]])/60,
                                                 elapsed_time = running_time))
  }
  
}

pdf("resources_profiling_downstreamAnalyses_final.pdf", w=6, h=4)

prof_summary$method <- gsub("MetaCell", "MC2",prof_summary$method)
colors <- c("MC2" = "#FF9900",
            "MC2+downstream" = "#FF6600", 
            "SEACells" = "#FF99CC",
            "SEACells+downstream" = "#FF3399",
            "SuperCell" = "#99CCFF",
            "SuperCell+downstream" = "#3399FF",
            "sc_downstream" = "black",
            "scBPCells_downstream" = "black")
shapes <- c("SEACells" = 16,
            "MC2" = 15, 
            "SuperCell" = 17,
            "SEACells+downstream" = 16,
            "MC2+downstream" = 15, 
            "SuperCell+downstream" = 17,
            "sc_downstream" = 18,
            "scBPCells_downstream" = 8)

p1 <- ggplot(prof_summary, aes(x=nb_cells, y=cpu_time, colour=method, group = method, shape = method)) + 
  geom_line() + 
  geom_point(size = 2) + theme_bw() + 
  xlab("Number of cells") + ylab("Time (minutes)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colors) + scale_shape_manual(values = shapes) +
  theme(legend.title=element_blank()) + ggtitle("Downstream analysis (Time)")  

p1 

p2 <- ggplot(prof_summary, aes(x=nb_cells, y=MaxRSS, colour=method, group = method, shape = method)) + 
  geom_line() + 
  geom_point(size = 2) + xlab("Number of cells") +  ylab("Max RSS (Gb)") + theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colors) + scale_shape_manual(values = shapes) + ggtitle("Downstream analysis (Memory)")+
  theme(legend.title=element_blank())
p2 

dev.off()

p1 + p2
