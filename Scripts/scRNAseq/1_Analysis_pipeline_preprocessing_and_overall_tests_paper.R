###########################
# Prepare environment
set.seed(1234567)
options(stringsAsFactors = FALSE)

# DLLs controls
length(getLoadedDLLs()) 

###########################
# LOAD -> Libraries & Functions
###########################
source("seurat_melanoma_timeseries_libraries.R")
source("seurat_melanoma_timeseries_functions.R")
source("seurat_melanoma_timeseries_charts.R")

library("cowplot")
library("Matrix")
library("data.table")

###########################
# LOAD DATA
###########################
# Set your base data path
base_data_path <- "your_path_here"

# Dataset 1 & 1_re (20_days)
data_rna_control_20_days <- Read10X(data.dir = file.path(base_data_path, "output_control_20_days/job/outs/filtered_feature_bc_matrix/"))
data_rna_control_re_20_days <- Read10X(data.dir = file.path(base_data_path, "output_control_20_re_days/job/outs/filtered_feature_bc_matrix/"))

data_rna_ba_20_days <- Read10X(data.dir = file.path(base_data_path, "output_ba_20_days/job/outs/filtered_feature_bc_matrix/"))
data_rna_ba_re_20_days <- Read10X(data.dir = file.path(base_data_path, "output_ba_20_re_days/job/outs/filtered_feature_bc_matrix/"))

# Dataset 2 (7_days)
data_rna_control_7_days <- Read10X(data.dir = file.path(base_data_path, "output_control_7_days/job/outs/filtered_feature_bc_matrix/"))
data_rna_ba_7_days  <- Read10X(data.dir = file.path(base_data_path, "output_ba_7_days/job/outs/filtered_feature_bc_matrix/"))

# Dataset 3 (13_days)
data_rna_control_13_days <- Read10X(data.dir = file.path(base_data_path, "output_control_13_days/job/outs/filtered_feature_bc_matrix/"))
data_rna_control_re_13_days <- Read10X(data.dir = file.path(base_data_path, "output_control_13_re_days/job/outs/filtered_feature_bc_matrix/"))

data_rna_ba_13_days <- Read10X(data.dir = file.path(base_data_path, "output_ba_13_days/job/outs/filtered_feature_bc_matrix/"))
data_rna_ba_re_13_days <- Read10X(data.dir = file.path(base_data_path, "output_ba_13_re_days/job/outs/filtered_feature_bc_matrix/"))

# Dataset 4 (Basal)
data_rna_basal <- Read10X(data.dir = file.path(base_data_path, "output_basal/job/outs/filtered_feature_bc_matrix/"))


###########################
# FINAL SAMPLES SELECTION
###########################
# combined_seurat should be created/loaded before this
# e.g., combined_seurat <- CreateSeuratObject(...)

combined_seurat_use <- SubsetData(object = combined_seurat, ident.use = c("0","1","2","3","4","5","6","7","8","9","10","11","12","14","15"))

# Subset per protocol
data_rna_control_20_days <- combined_seurat_use@assays$RNA@counts[, combined_seurat_use@meta.data$protocol == "data_rna_control_20_days"]
data_rna_ba_20_days <- combined_seurat_use@assays$RNA@counts[, combined_seurat_use@meta.data$protocol == "data_rna_ba_20_days"]

data_rna_control_7_days <- combined_seurat_use@assays$RNA@counts[, combined_seurat_use@meta.data$protocol == "data_rna_control_7_days"]
data_rna_ba_7_days <- combined_seurat_use@assays$RNA@counts[, combined_seurat_use@meta.data$protocol == "data_rna_ba_7_days"]

data_rna_control_13_days <- combined_seurat_use@assays$RNA@counts[, combined_seurat_use@meta.data$protocol == "data_rna_control_13_days"]
data_rna_ba_13_days <- combined_seurat_use@assays$RNA@counts[, combined_seurat_use@meta.data$protocol == "data_rna_ba_13_days"]

rm(list = setdiff(ls(), c("data_rna_control_20_days","data_rna_ba_20_days",
                          "data_rna_control_7_days","data_rna_ba_7_days",
                          "data_rna_control_13_days","data_rna_ba_13_days",
                          "data_rna_basal", "combined_seurat")))
gc()

save.image("temporal.Rdata")


######################################################
#Start downstream analysis
######################################################
#Working directory
#setwd("/datos_2/Analysis/New_full_analysis/all_samples_together_with_out_noise_cells/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/sample_20/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/sample_20/cleaned")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/sample_13/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/sample_all/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/sample_final/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/sample_final_no_basal/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_remove_day20_low_quality/")
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/")
###########################################################

#AGGREGATE COUNTS OF RESEQ OF SAME SAMPLES (FOR 20 AND 13 DAYS)
cell_overlap_20_days_control= colnames(data_rna_control_20_days)[colnames(data_rna_control_20_days) %in% colnames(data_rna_control_re_20_days)]
data_rna_control_20_days = data_rna_control_20_days[,cell_overlap_20_days_control]
data_rna_control_re_20_days = data_rna_control_re_20_days[,cell_overlap_20_days_control]
data_rna_control_20_days = data_rna_control_20_days + data_rna_control_re_20_days

cell_overlap_20_days_ba= colnames(data_rna_ba_20_days)[colnames(data_rna_ba_20_days) %in% colnames(data_rna_ba_re_20_days)]
data_rna_ba_20_days = data_rna_ba_20_days[,cell_overlap_20_days_ba]
data_rna_ba_re_20_days = data_rna_ba_re_20_days[,cell_overlap_20_days_ba]
data_rna_ba_20_days = data_rna_ba_20_days + data_rna_ba_re_20_days

cell_overlap_13_days_control= colnames(data_rna_control_13_days)[colnames(data_rna_control_13_days) %in% colnames(data_rna_control_re_13_days)]
data_rna_control_13_days = data_rna_control_13_days[,cell_overlap_13_days_control]
data_rna_control_re_13_days = data_rna_control_re_13_days[,cell_overlap_13_days_control]
data_rna_control_13_days = data_rna_control_13_days + data_rna_control_re_13_days

cell_overlap_13_days_ba= colnames(data_rna_ba_13_days)[colnames(data_rna_ba_13_days) %in% colnames(data_rna_ba_re_13_days)]
data_rna_ba_13_days = data_rna_ba_re_13_days[,cell_overlap_13_days_ba]
data_rna_ba_re_13_days= data_rna_ba_re_13_days[,cell_overlap_13_days_ba]
data_rna_ba_13_days = data_rna_ba_13_days + data_rna_ba_re_13_days

#ALL together
filename="project"
file_list = c("data_rna_control_20_days","data_rna_ba_20_days","data_rna_control_13_days","data_rna_ba_13_days", "data_rna_control_7_days","data_rna_ba_7_days", "data_rna_basal")
file_list_mat = c("data_rna_control_20_days","data_rna_ba_20_days","data_rna_control_13_days","data_rna_ba_13_days", "data_rna_control_7_days","data_rna_ba_7_days", "data_rna_basal")
list_all = c(data_rna_control_20_days=SingleCellExperiment(assays = list(counts = data_rna_control_20_days)),
             data_rna_ba_20_days=SingleCellExperiment(assays = list(counts = data_rna_ba_20_days)),
             data_rna_control_13_days=SingleCellExperiment(assays = list(counts = data_rna_control_13_days)),
             data_rna_ba_13_days=SingleCellExperiment(assays = list(counts = data_rna_ba_13_days)),
             data_rna_control_7_days=SingleCellExperiment(assays = list(counts = data_rna_control_7_days)),
             data_rna_ba_7_days=SingleCellExperiment(assays = list(counts = data_rna_ba_7_days)),
             data_rna_basal=SingleCellExperiment(assays = list(counts = data_rna_basal))
)

#Form day_X and day_X_re
file_list = c("data_rna_control","data_rna_control_re","data_rna_BA","data_rna_BA_re")
file_list_mat = c("data_rna_control","data_rna_control_re","data_rna_BA","data_rna_BA_re")
list_all = c(data_rna_control_20_days=SingleCellExperiment(assays = list(counts = data_rna_control)),
             data_rna_control_20_days_re=SingleCellExperiment(assays = list(counts = data_rna_control_re)),
             data_rna_ba_20_days=SingleCellExperiment(assays = list(counts = data_rna_BA)),
             data_rna_ba_20_days_re=SingleCellExperiment(assays = list(counts = data_rna_BA_re))
)

#ADD rowData the rownames for genes
list_all_Seurat = c()
rm(list = file_list)
rm(list = file_list_mat)
rm(file_list_mat)

###########################################################
#PRE-filtering
###########################################################
#Replace possible NAs by 0 in count matrix
for (n in 1:length(file_list)){
  counts(list_all[[n]])[is.na(counts(list_all[[n]]))] <- 0
}

#Remove Feature(genes,rna,etc) not expressed in any cell
for (n in 1:length(file_list)){
  keep_feature = rowSums(as.matrix(counts(list_all[[n]])) > 0) > 0
  list_all[[n]] = list_all[[n]][keep_feature, ]
}
rm(keep_feature)
#######################################################
#Define with features(genes,rna,etc) are the ERCC spike-ins and mitochondrial genes
for (n in 1:length(file_list)){
  isSpike(list_all[[n]], "ERCC") = grepl("ERCC-", rownames(list_all[[n]]))
  isSpike(list_all[[n]], "MT") = grepl("MT-", rownames(list_all[[n]]))
}
gc()

#######################################################
#Calculate de quality metrics
for (n in 1:length(file_list)){
  list_all[[n]]  = calculateQCMetrics(list_all[[n]],
                                      feature_controls = list(
                                        ERCC = isSpike(list_all[[n]], "ERCC"), 
                                        MT = isSpike(list_all[[n]], "MT")
                                      )
  )
  print (sprintf("QC added to sample -> %s", file_list[n]))
}
###########################################################
#PRE-PROCESSING -> FILTERING BAD QUALITY INSTANCES AND VARIABLES
###########################################################
### total_counts / total_features  #

#for day 20: 3250 nCount and 1500 nFeature (SPECIFIC EFFECT THING)
cells_to_preserve_by_counts = colnames(list_all[[1]])[list_all[[1]][["total_counts"]] > 3250]
cells_to_preserve_by_features = colnames(list_all[[1]])[list_all[[1]][["total_features_by_counts"]] > 1500]
cells_to_preserve = cells_to_preserve_by_counts[cells_to_preserve_by_counts %in% cells_to_preserve_by_features]
list_all[[1]] = list_all[[1]][,cells_to_preserve]

#rest as usual
lista_dataframes_counts = filter_instances_variables("total_counts", "Filtering_total_counts_plots_Histograms.pdf", "Filtering_total_counts_plots_Tables.pdf")
lista_dataframes_features = filter_instances_variables("total_features_by_counts", "Filtering_total_features_plots_Histograms.pdf", "Filtering_total_features_plots_Tables.pdf")

#######################################################
###ERCC/ MT 
lista_dataframes_ERCC = filter_ERCC_MT("is_spike_ERCC")
#lista_dataframes_MT = filter_ERCC_MT("is_spike_MT")

#######################################################
###Apply all previous created filters and filter umi based on those filters
for (n in 1:length(file_list)){
  #unlist(lista_dataframes_features[n]) &
  list_all[[n]]$use = (unlist(lista_dataframes_features[n]) & unlist(lista_dataframes_counts[n]))
  print (table(list_all[[n]]$use))
}
gc()

#######################################################
#Gene filtering after cell filtering !!!
#Low quality genes + MT + ERCC
#[x > 0] 1 before
lista_filter_genes = list()
for (n in 1:length(file_list)){
  #
  filter_genes = apply(
    counts(list_all[[n]][ , colData(list_all[[n]])$use]), 
    1, 
    function(x) length(x[x > 0]) >= 2
  )
  lista_filter_genes = append(lista_filter_genes, list(filter_genes))
  #Filter rows (genes)
  #& !unlist(lista_dataframes_MT[n]) 
  rowData(list_all[[n]])$use = (unlist(lista_filter_genes[n]) & !unlist(lista_dataframes_ERCC[n]))
  print(table(rowData(list_all[[n]])$use))
}

#lista_dataframes_MT
rm(list=setdiff(ls(), c("list_all", "list_all_Seurat" ,"file_list", "filename")))
gc()

####
#PROCESSING -> COMBINE DATASETS SEURAT CCA 
####
#######################################################
#Creating Seurat objects
for (n in 1:length(file_list)){
  list_all_Seurat = append(list_all_Seurat, CreateSeuratObject(
    counts = counts(list_all[[n]][rowData(list_all[[n]])$use , colData(list_all[[n]])$use]), min.cells = 3, min.features = 200
  )) 
}
rm(list_all)


#Adding % mit and % rib information for each cell
for (n in 1:length(file_list)){
  list_all_Seurat[[n]][["percent.mt"]] <- PercentageFeatureSet(list_all_Seurat[[n]], pattern = "^MT-")
  list_all_Seurat[[n]][["percent.rb"]] <- PercentageFeatureSet(list_all_Seurat[[n]], pattern = c("^RPS", "^RPL"))
}

#######################################################
#Visualize data for filtering propouse
seurat_data_plot(paste0(filename,"Seurat_data_filtered.pdf"))

#REMOVE CELLS HAVING MORE THAN 5% MIT GENES
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = subset(x = list_all_Seurat[[n]], percent.mt < 10)
#  list_all_Seurat[[n]] = subset(x = list_all_Seurat[[n]], percent.rb < 5)
}

#######################################################
#Visualize data for filtering propouse
seurat_data_plot(paste0(filename,"Seurat_data_filtered_2.pdf"))

#######################################################
#Normalization and log2 transformation
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = NormalizeData(
    object = list_all_Seurat[[n]], 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
}

#######################################################
#Scale data
#usually regress nUMI if possible 
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = ScaleData(object = list_all_Seurat[[n]], 
                                   vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt")
  )
}

#######################################################
#Highly variable genes obtain for downstream analysis
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = FindVariableFeatures(
    object = list_all_Seurat[[n]],
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    #x.low.cutoff = 0.0125, 
    #x.high.cutoff = 3, 
    #y.cutoff = 0.5,
    nfeatures = 2000,
    mean.cutoff = c(0.1, 8),
    dispersion.cutoff = c(1, Inf),
    verbose = FALSE
  )
}

#######################################################
#Run PCA on each sample
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = RunPCA(object = list_all_Seurat[[n]], npcs = 30,verbose = FALSE)
}

#######################################################
#GENERATE SCT NORMALIZATION OUTCOMES

#ONLY if SCT assay its used for doublets detection
#for (n in 1:length(file_list)){
#  list_all_Seurat[[n]] = SCTransform(
#    object = list_all_Seurat[[n]],
#    assay = "RNA",
#    new.assay.name = "SCT",
#    vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"),
#  )
#}
#Run PCA on each sample
#for (n in 1:length(file_list)){
#  list_all_Seurat[[n]] = RunPCA(object = list_all_Seurat[[n]], npcs = 30,verbose = FALSE)
#}

#######################################################
#Set the protocol and name for easy identification
for (n in 1:length(file_list)){
  list_all_Seurat[[n]]@meta.data[, "protocol"] = file_list[n]
}

list_all_Seurat[[1]]@meta.data[, "protocol_2"] = "control"
list_all_Seurat[[2]]@meta.data[, "protocol_2"] = "ba"
list_all_Seurat[[3]]@meta.data[, "protocol_2"] = "control"
list_all_Seurat[[4]]@meta.data[, "protocol_2"] = "ba"
list_all_Seurat[[5]]@meta.data[, "protocol_2"] = "control"
list_all_Seurat[[6]]@meta.data[, "protocol_2"] = "ba"
list_all_Seurat[[7]]@meta.data[, "protocol_2"] = "basal"

#######################################################
#Find Doblets
library(DoubletFinder)
list_all_Seurat

#CHANGE TO RNA  // IF assay SCT its used for doublets and has been computed previously change to SCT
for (n in 1:length(file_list)){
  DefaultAssay(list_all_Seurat[[n]]) <- "RNA"
}

#DO IT ON EACH SAMPLE INDEPENDENTLY NEED FINDVAR AND PCA RUNNED ON EACH SAMPLE BEFORE IT
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
bcmvn =c()
for (n in 1:length(file_list)){
  sweep.res <- paramSweep_v3(list_all_Seurat[[n]], PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn[[n]] <- find.pK(sweep.stats)
  bcmvn[[n]]=data.frame(bcmvn[[n]])
}

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
#sweep.res.list_combined_seurat <- paramSweep_v3(list_all_Seurat[[1]], PCs = 1:10, sct = FALSE)
#gt.calls <- list_all_Seurat[[1]]@meta.data[rownames(sweep.res.list_combined_seurat[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
#sweep.stats_combined_seurat<- summarizeSweep(sweep.res.list_combined_seurat, GT = TRUE, GT.calls = gt.calls)
#bcmvn_combined_seurat <- find.pK(sweep.stats_combined_seurat)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
for (n in 1:length(file_list)){
  annotations = list_all_Seurat[[n]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)                  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(list_all_Seurat[[n]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  pk_value = as.numeric(as.vector(bcmvn[[n]][which.max(bcmvn[[n]]$BCmetric),]$pK))
  to_replace = colnames(list_all_Seurat[[n]]@meta.data)[length(colnames(list_all_Seurat[[n]]@meta.data))-1]
  
  list_all_Seurat[[n]] <- doubletFinder_v3(list_all_Seurat[[n]], PCs = 1:15, pN = 0.25, pK = pk_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  list_all_Seurat[[n]] <- doubletFinder_v3(list_all_Seurat[[n]], PCs = 1:15, pN = 0.25, pK = pk_value, nExp = nExp_poi.adj, reuse.pANN = to_replace, sct = FALSE)
}

pdf(file="Distribution_predicted_doublets.pdf", width=12, height=6)
for (n in 1:length(file_list)){
  p <- list_all_Seurat[[n]]@meta.data %>%
    ggplot( aes(x=nCount_RNA, fill= as.factor(as.character(unlist(list_all_Seurat[[n]]@meta.data[ncol(list_all_Seurat[[n]]@meta.data)]))) )) +
      geom_histogram(alpha=0.2, position = 'identity', binwidth = 30) +
      scale_fill_manual(values=c("#FF0000", "#404080")) +
      labs(fill="")
  plot(p)
}
dev.off()

#Select remaining cells for downstream analysis without PREDICTED doublets
for (n in 1:length(file_list)){
  to_filter = colnames(list_all_Seurat[[n]]@meta.data)[length(colnames(list_all_Seurat[[n]]@meta.data))]
  
  cells.use <- colnames(list_all_Seurat[[n]])[which(list_all_Seurat[[n]][[]][to_filter] != "Doublet")]
  list_all_Seurat[[n]] <- subset(list_all_Seurat[[n]], cells = cells.use)
}

#####################################################################################################################################################################
#####################################################################################################################################################################
#Redo analysis without doblets

#######################################################
#Normalization and log2 transformation
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = NormalizeData(
    object = list_all_Seurat[[n]], 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
}

#######################################################
#Scale data
#usually regress nUMI if possible 
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = ScaleData(object = list_all_Seurat[[n]], 
                                   vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt")
  )
}

#######################################################
#Highly variable genes obtain for downstream analysis
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = FindVariableFeatures(
    object = list_all_Seurat[[n]],
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    #x.low.cutoff = 0.0125, 
    #x.high.cutoff = 3, 
    #y.cutoff = 0.5,
    nfeatures = 2000,
    mean.cutoff = c(0.1, 8),
    dispersion.cutoff = c(1, Inf),
    verbose = FALSE
  )
}

#######################################################
#Run PCA on each sample
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = RunPCA(object = list_all_Seurat[[n]], npcs = 30,verbose = FALSE)
}

#######################################################
#GENERATE SCT NORMALIZATION OUTCOMES
for (n in 1:length(file_list)){
  list_all_Seurat[[n]] = SCTransform(
    object = list_all_Seurat[[n]],
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"),
  )
}

#######################################################
#Set the protocol and name for easy identification
for (n in 1:length(file_list)){
  list_all_Seurat[[n]]@meta.data[, "protocol"] = file_list[n]
}

list_all_Seurat[[1]]@meta.data[, "protocol_2"] = "control"
list_all_Seurat[[2]]@meta.data[, "protocol_2"] = "ba"
list_all_Seurat[[3]]@meta.data[, "protocol_2"] = "control"
list_all_Seurat[[4]]@meta.data[, "protocol_2"] = "ba"
list_all_Seurat[[5]]@meta.data[, "protocol_2"] = "control"
list_all_Seurat[[6]]@meta.data[, "protocol_2"] = "ba"
list_all_Seurat[[7]]@meta.data[, "protocol_2"] = "basal"


save.image("temporal.Rdata")

#Compute vectors, sustract and append datasets into one dataset
#UNLY FOR ONE SAMPLE UNCOMMENT LINE DOWN AND COMMENT REST UNTIL scale and pca for umap sentence
#combined_seurat = list_all_Seurat[[1]]
##############################################################################################################
##############################################################################################################
##############################################################################################################
#Run CCA with ANCHORS to combine and integrate datasets finding diff vectors and sustracting to each dataset)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 10240^2)

#NEEDED if we make use of SCT 
list_all_Seurat_features <- SelectIntegrationFeatures(object.list = list_all_Seurat, nfeatures = 3000)
list_all_Seurat <- PrepSCTIntegration(object.list = list_all_Seurat, anchor.features = list_all_Seurat_features, verbose = FALSE)

#Find anchors
#For not SCT
#anchors <- FindIntegrationAnchors(object.list = list_all_Seurat, dims = 1:20)
#for SCT
anchors <- FindIntegrationAnchors(object.list = list_all_Seurat, normalization.method = "SCT", anchor.features = list_all_Seurat_features, dims = 1:30, verbose = FALSE)


#Selecting overlaping genes between datasets
# NOW THIS FUNCTION EXIST FOR DOING THIS WITH SEURAT: SelectIntegrationFeatures
overlap_genes = Reduce(intersect, list(rownames(list_all_Seurat[[1]])), (rownames(list_all_Seurat[[2]])))
if (length(file_list) > 2){
  for (n in 3:length(file_list)){
    overlap_genes = Reduce(intersect, list(overlap_genes, (rownames(list_all_Seurat[[n]]))))
  }
}

#change overlap_genes by list_all_Seurat_features if speed its needed
combined_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = list_all_Seurat_features, dims = 1:30, verbose = TRUE)

##############################################################################################################
##############################################################################################################
##############################################################################################################

#DO UMAP AND CLUSTERING ON INTEGRATED ASSAY
#######################################################
#Scale and PCA for UMAP

#ONLY IF SCT WAS NOT USED
#combined_seurat = ScaleData(object = combined_seurat, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), verbose = TRUE)

combined_seurat = RunPCA(object = combined_seurat, npcs = 30,verbose = FALSE)

#######################################################
#PC selection
combined_seurat = JackStraw(
  object = combined_seurat,
  num.replicate = 100
)
combined_seurat = ScoreJackStraw(object = combined_seurat, dims = 1:20)
#Plot
pdf(file="PC_statisticPower.pdf", width=12, height=6)
p1 = JackStrawPlot(object = combined_seurat, dims = 1:20)
p2 = ElbowPlot(combined_seurat)
plot_grid(p1, p2)
dev.off()

pdf(file="Corrected_DATA.pdf", width=12, height=6)
p1 = DimPlot(object = combined_seurat, reduction = "pca")
p2 = VlnPlot(combined_seurat, features = "PC_1")
plot_grid(p1, p2)
dev.off()

#######################################################
#Automaticaly select the best CC amount based in MetageneBicorPlot result
#Metagene_result = split(p2$data$stdev, ceiling(seq_along(p2$data$stdev)/25))
#cc_list = c()
#cc_list = append(cc_list ,min(which(Metagene_result[[1]] < quantile(Metagene_result[[1]], probs = seq(0,1,0.1))[[4]])))
cc_number = 15
#rm(Metagene_result)

#######################################################
#Generate UMAP
combined_seurat = RunUMAP(object = combined_seurat, reduction = "pca", dims = 1:30)

#######################################################
#######################################################
#######################################################
#IKAP CLUSTERING METHOD
#sobj <- IKAP(combined_seurat, out.dir = "./IKAP")
#######################################################
#######################################################
#######################################################

combined_seurat = FindNeighbors(object = combined_seurat, dims = 1:30)
combined_seurat = FindClusters(object = combined_seurat, reduction.type = "pca",dims = 1:30, resolution = 0.6)

#######################################################
#PLOT UMAP
pdf(file="UMAP_clustered_yyy.pdf", width=12, height=6)
p1 = DimPlot(combined_seurat, reduction ="umap", label = TRUE, do.return = T, pt.size = 0.5, do.label = T)
plot_grid(p1)
dev.off()

#######################################################
#PLOT 2 UMAPS together
pdf(file="UMAP_clustered_2.pdf", width=12, height=6)
p1 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "protocol")
p2 = DimPlot(object = combined_seurat, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

#Plot MIT and RIB
pdf(file="Mito_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "percent.mt") 
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = combined_seurat@meta.data$percent.mt)) +
    scale_colour_gradient2()
dev.off()

pdf(file="Ribo_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "percent.rb") 
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = combined_seurat@meta.data$percent.rb)) +
    scale_colour_gradient2()
dev.off()

#########################################################################################
#########################################################################################
#FINDMARKERS all clusters agains rest of cells
DefaultAssay(combined_seurat) <- "RNA"

list_RNAseq_data_raw.markers = FindAllMarkers(object = combined_seurat, only.pos = TRUE, min.pct = 0.25)
print(x = head(x = list_RNAseq_data_raw.markers, n = 5))
list_RNAseq_data_raw.markers = (subset(list_RNAseq_data_raw.markers,(list_RNAseq_data_raw.markers['p_val_adj'] <= 0.01 )))

write.csv(list_RNAseq_data_raw.markers, file = "Integrated_Markers.csv")

#######################################################
#For heatmap
#######################################################
#Scale data
#usually regress nUMI if possible 
combined_seurat = ScaleData(object = combined_seurat, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"))

#Plot markers for each cluster
library(dplyr)
list_RNAseq_data_raw.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10

pdf(file="Integrated_Markers.pdf", width=16, height=20)
DoHeatmap(object = combined_seurat, features = top10$gene) + NoLegend()
dev.off()

save.image("temporal.Rdata")
#########################################################################################
#########################################################################################
#########################################################################################
#Cell type enrich of each cluster
setwd("/datos_2/datos_temporal/datos/Bcell/RAW_DATA/Bcells_scRNAseq/cellType_enrich/R")
working_dir = "/datos_2/datos_temporal/datos/Bcell/RAW_DATA/Bcells_scRNAseq/cellType_enrich/data/"

source("enrich.R")
source("make_adjacency_matrix.R")

#Create the cluster file
cluster=subset(list_RNAseq_data_raw.markers, select=c("p_val_adj", "avg_logFC", "cluster", "gene"))
#Export
write.csv(cluster, file = "cluster.csv")

setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/")
###
#CELL enrich for each cluster from TSNE marker genes (SEURAT)
pdf(file="Clusters_cell_type_enrich.pdf", width = 9, height = 9)
res_enrich <- enrich(clusters = "/datos_2/datos_temporal/datos/Bcell/RAW_DATA/Bcells_scRNAseq/cellType_enrich/R/cluster.csv", annoDB = 'xCell')
dev.off()

#######################################
#######################################
#Subseting clusters not low expressed (bad quality)
pdf(file="nCount_RNA_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "nCount_RNA") 
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = combined_seurat@meta.data$nCount_RNA)) +
    scale_colour_gradient2()
dev.off()

pdf(file="nFeature_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "nFeature") 
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = combined_seurat@meta.data$nFeature_RNA)) +
    scale_colour_gradient2()
dev.off()

pdf(file="Mito_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "percent.mt") 
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = combined_seurat@meta.data$percent.mt)) +
    scale_colour_gradient2()
dev.off()

pdf(file="Ribo_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = combined_seurat, reduction = "umap", group.by = "percent.rb") 
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = combined_seurat@meta.data$percent.rb)) +
    scale_colour_gradient2()
dev.off()

#Select cells to be removed (from mit removed experiment)
combined_seurat_to_remove_cells = SubsetData(object = combined_seurat, ident.use=c("0","10","13")) #Quitar(0,10,13)
control_cells_to_by_removed = rownames(combined_seurat_to_remove_cells@meta.data)[combined_seurat_to_remove_cells@meta.data$protocol=="data_rna_control"]
control_re_cells_to_by_removed = rownames(combined_seurat_to_remove_cells@meta.data)[combined_seurat_to_remove_cells@meta.data$protocol=="data_rna_control_re"]

BA_cells_to_by_removed =rownames(combined_seurat_to_remove_cells@meta.data)[combined_seurat_to_remove_cells@meta.data$protocol=="data_rna_BA"]
BA_re_cells_to_by_removed =rownames(combined_seurat_to_remove_cells@meta.data)[combined_seurat_to_remove_cells@meta.data$protocol=="data_rna_BA_re"]

'%ni%' <- Negate('%in%')
data_rna_control = data_rna_control[,colnames(data_rna_control)[colnames(data_rna_control) %ni% sub("\\_.*", "", control_cells_to_by_removed)]]
data_rna_control_re = data_rna_control_re[,colnames(data_rna_control_re)[colnames(data_rna_control_re) %ni% sub("\\_.*", "", control_re_cells_to_by_removed)]]

data_rna_BA = data_rna_BA[,colnames(data_rna_BA)[colnames(data_rna_BA) %ni% sub("\\_.*", "", BA_cells_to_by_removed)]]
data_rna_BA_re = data_rna_BA_re[,colnames(data_rna_BA_re)[colnames(data_rna_BA_re) %ni% sub("\\_.*", "", BA_re_cells_to_by_removed)]]

#######################################
#######################################
#plot proportions for each clsuter
library(Seurat) 
library(dplyr)
library(ggplot2)

meta.data <- combined_seurat@meta.data

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(count = n())
pdf(file="Proportion_1_2.pdf", width=12, height=6)
ggplot(counts, aes(seurat_clusters, count, fill = protocol)) +
  geom_bar(stat = 'identity')
dev.off()

counts$count = (counts$count / nrow(meta.data))*100
pdf(file="Proportion_2.pdf", width=12, height=6)
ggplot(counts, aes(seurat_clusters, count, fill = protocol)) +
  geom_bar(stat = 'identity')
dev.off()

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(percent_mt = mean(percent.mt))
pdf(file="Proportion_MT_3.pdf", width=12, height=6)
ggplot(counts, aes(seurat_clusters, percent_mt, fill = protocol)) +
  geom_bar(stat = 'identity')
dev.off()

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(nFeature_RNA = mean(nFeature_RNA))
pdf(file="Proportion_Features_4.pdf", width=12, height=6)
ggplot(counts, aes(seurat_clusters, nFeature_RNA, fill = protocol)) +
  geom_bar(stat = 'identity')
dev.off()

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(nCount_RNA = mean(nCount_RNA))
pdf(file="Proportion_Counts_5.pdf", width=12, height=6)
ggplot(counts, aes(seurat_clusters, nCount_RNA, fill = protocol)) +
  geom_bar(stat = 'identity')
dev.off()

#Cell_cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

marrow <- CellCycleScoring(combined_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pdf(file="cell_cycle_corrected.pdf", width=12, height=6)
DimPlot(marrow)
dev.off()

marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
pdf(file="cell_cycle_corrected_2.pdf", width=12, height=6)
DimPlot(marrow)
dev.off()

marrow <- CellCycleScoring(combined_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(marrow)
dev.off()
marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

marrow = RunUMAP(object = marrow, reduction = "pca", dims = 1:cc_number)
pdf(file="UMAP_cell_cycle_corrected.pdf", width=12, height=6)
p1 = DimPlot(object = marrow, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "integrated_snn_res.0.6")
plot_grid(p1)
dev.off()

#######################################
#######################################
#Diff Expression Analaysis
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/control_vs_BA/")
#For each cluster
#For each time points
#Control vs BA
library(dplyr)

#ADD vector to metadata with time information for filtering
length(grep("20",(combined_seurat@meta.data$protocol)))
length(grep("13",(combined_seurat@meta.data$protocol)))
length(grep("7",(combined_seurat@meta.data$protocol)))


vector_days = c()
vector_days = append(vector_days,(rep("20_days", length(grep("20",(combined_seurat@meta.data$protocol))))))
vector_days = append(vector_days,(rep("13_days", length(grep("13",(combined_seurat@meta.data$protocol))))))
vector_days = append(vector_days,(rep("7_days", length(grep("7",(combined_seurat@meta.data$protocol))))))
vector_days = append(vector_days,(rep("Basal_days", length(grep("basal",(combined_seurat@meta.data$protocol))))))
combined_seurat@meta.data$vector_days = vector_days

Idents(combined_seurat) <- "vector_days"
for (n in 1:length(levels(Idents(combined_seurat)))) { # for each time_point 20/13/7
  cells_time_point_n = SubsetData(object = combined_seurat, ident.use=c(levels(Idents(combined_seurat))[n]))
  for (m in 1:length(levels(cells_time_point_n@meta.data$seurat_clusters))-1) { # for each cluster 0-17
    Idents(cells_time_point_n) <- "seurat_clusters"
    cells_cluster_n = SubsetData(object = cells_time_point_n, ident.use=c(m))
    vector_condition = c()
    vector_condition = append(vector_condition,(rep("control", length(grep("control",(cells_cluster_n@meta.data$protocol))))))
    vector_condition = append(vector_condition,(rep("ba", length(grep("ba",(cells_cluster_n@meta.data$protocol))))))
    cells_cluster_n@meta.data$vector_condition = vector_condition
    if (length(cells_cluster_n@meta.data$vector_condition) > 0){
      if (table(cells_cluster_n@meta.data$vector_condition)[1] > 3 & table(cells_cluster_n@meta.data$vector_condition)[2] > 3){
        Idents(cells_cluster_n) <- "vector_condition"
        resulting.markers = FindMarkers(object = cells_cluster_n, ident.1 = "control", ident.2 = "ba", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.125)
        resulting.markers = (subset(resulting.markers,(resulting.markers['p_val_adj'] <= 0.01 )))
        if (nrow(resulting.markers) > 0){
          write.csv(resulting.markers, file = paste0("Integrated_",levels(Idents(combined_seurat))[n],"_cluster_",m,".csv"))
          #
          resulting.markers$cluster = m
          resulting.markers$gene = rownames(resulting.markers)
          resulting.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> top100
          #pdf(file= paste0("Integrated_",levels(Idents(combined_seurat))[n],"_cluster_",m,"_heatmap_top100",".pdf"), width=16, height=20)
          #plot(DoHeatmap(object = cells_cluster_n, features = top100$gene) + NoLegend())
          #dev.off()
        }
      }
    }
  }
}

#######################################
#######################################    
#GENE SET ANALYSIS DONE BY CLUSTERPROFILER
library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")
library("enrichplot")
library("ggnewscale")

#OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/GSA/GO/"
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/GSA/GO/"
for (m in 1:length(levels(combined_seurat@meta.data$seurat_clusters))-1) {
  markers_cluster_m = (subset(list_RNAseq_data_raw.markers,(list_RNAseq_data_raw.markers['cluster'] == m )))
  
  #Create the geneList
  geneList <- markers_cluster_m$avg_logFC
  markers_cluster_m$avg_logFC = markers_cluster_m$avg_log2FC
  
  ont = c("CC", "MF", "BP")
  for (n in 1:length(ont)) {
    ego <- enrichGO(gene       = names(geneList),
                 universe      = rownames(combined_seurat),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = ont[n],
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
    ego_table <- as.data.frame(ego) 
    write.csv(ego_table,paste0(OutPath,"ego_",ont[n],"_cluster_",m,".csv"))
    #plots
    if (nrow(ego_table) > 0){
      #ego <- pairwise_termsim(ego) 
      pdf(file=paste0(OutPath,"ego__barplot_",ont[n],"_cluster_",m,".pdf"), width=8, height=8)  
      plot(barplot(ego, showCategory=20))  
      dev.off()
      pdf(file=paste0(OutPath,"ego_cnetplot_",ont[n],"_cluster_",m,".pdf"), width=8, height=8)  
      plot(cnetplot(ego, categorySize="pvalue", foldChange=geneList))
      dev.off()  
    }  
    #if (nrow(ego_table) > 1){
    #  pdf(file=paste0(OutPath,"ego_emaplot_",ont[n],"_cluster_",m,".pdf"), width=8, height=8)  
    #  plot(emapplot(ego, showCategory=20))
    #  dev.off()
    #}
  }  
}


#Change names to entrez code
name_changer = function(especie_database, especie_name, x){
  ensembl <- useMart(especie_database, dataset = especie_name, host = "useast.ensembl.org")
  annotation = getBM(attributes=c("hgnc_symbol","entrezgene_id"), values=x, mart=ensembl)
  return(annotation)
}

OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/GSA/KEGG/"
for (m in 1:length(levels(combined_seurat@meta.data$seurat_clusters))-1) {
  markers_cluster_m = (subset(list_RNAseq_data_raw.markers,(list_RNAseq_data_raw.markers['cluster'] == m )))
  
  #Create the geneList
  geneList <- markers_cluster_m$avg_logFC
  names(geneList) <- as.character(markers_cluster_m$gene)

  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]

  #gseKEGG
  ekegg <- enrichKEGG(gene      = names(geneList),
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  ekegg_table <- as.data.frame(ekegg)
  write.csv(ekegg_table,paste0(OutPath,"ekegg_cluster_",m,".csv"))
  #plots
  if (nrow(ekegg_table) > 0){
    #ekegg <- pairwise_termsim(ekegg)
    pdf(file=paste0(OutPath,"ekegg_barplot_cluster_",m,".pdf"), width=8, height=8)
    plot(barplot(ekegg, showCategory=20))            
    dev.off()
    #pdf(file=paste0(OutPath,"ekegg_emaplot_cluster_",m,".pdf"), width=8, height=8)  
    #plot(emapplot(ekegg))
    #dev.off() 
  }
}


OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/GSA/REACTOME/"
for (m in 1:length(levels(combined_seurat@meta.data$seurat_clusters))-1) {
  markers_cluster_m = (subset(list_RNAseq_data_raw.markers,(list_RNAseq_data_raw.markers['cluster'] == m )))
  
  #Create the geneList
  geneList <- markers_cluster_m$avg_logFC
  names(geneList) <- as.character(markers_cluster_m$gene)

  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]

  #REACTOME   
  library(ReactomePA)    
  reactome <- enrichPathway(gene=names(geneList),pvalueCutoff=0.05, readable=T)
  reactome_table <- as.data.frame(reactome)
  write.csv(reactome_table,paste0(OutPath,"reactome_cluster",m,".csv"))
  #plots
    if (nrow(reactome_table) > 0){
      pdf(file=paste0(OutPath,"reactome_barplot_cluster_",m,".pdf"), width=8, height=8)
      plot(barplot(reactome, showCategory=20))            
      dev.off()
      #pdf(file=paste0(OutPath,"reactome_emaplot_cluster_",m,".pdf"), width=8, height=8)  
      #plot(emapplot(reactome))
      #dev.off() 
    }
}


for (m in 1:length(levels(combined_seurat@meta.data$seurat_clusters))-1) {
  markers_cluster_m = (subset(list_RNAseq_data_raw.markers,(list_RNAseq_data_raw.markers['cluster'] == m )))
  
  #Create the geneList
  geneList <- markers_cluster_m$avg_logFC
  names(geneList) <- as.character(markers_cluster_m$gene)

  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]
  
  OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/GSA/MSIGDBR/H/"
  
  #MSigDb
  library(msigdbr)
  H = msigdbr(species = "Homo sapiens", category = "H")
  ont = H[,3]
  gene = H[,8]
  H_table = data.frame(ont,gene)
  colnames(H_table) = c("ont","gene")
  H_signatures <- enricher(gene=names(geneList), TERM2GENE=H_table, pvalueCutoff = 0.05)
  H_signatures_table = data.frame(H_signatures)
  write.csv(H_signatures_table,paste0(OutPath,"H_signatures_cluster",m,".csv"))
  #plots
  if (nrow(H_signatures_table) > 0){    
    pdf(file=paste0(OutPath,"H_signatures_barplot_cluster_",m,".pdf"), width=8, height=8)            
    plot(barplot(H_signatures))
    dev.off()
    #pdf(file=paste0(OutPath,"H_signatures_emapplot_cluster_",m,".pdf"), width=8, height=8)  
    #plot(emapplot(H_signatures))
    #dev.off() 
  }
  
  OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/GSA/MSIGDBR/C6/"
  
  c6 = msigdbr(species = "Homo sapiens", category = "C6")
  ont = c6[,3]
  gene = c6[,8]
  c6_table = data.frame(ont,gene)
  colnames(c6_table) = c("ont","gene")
  oncogenic <- enricher(names(geneList), TERM2GENE=c6_table, pvalueCutoff = 0.05)
  oncogenic_table = data.frame(oncogenic)
  write.csv(oncogenic_table,paste0(OutPath,"oncogenic_cluster",m,".csv"))
  #plots
  if (nrow(oncogenic_table) > 0){    
    pdf(file=paste0(OutPath,"oncogenic_barplot_cluster_",m,".pdf"), width=8, height=8)            
    plot(barplot(oncogenic))
    dev.off()
    #pdf(file=paste0(OutPath,"oncogenic_emapplot_cluster_",m,".pdf"), width=8, height=8)  
    #plot(emapplot(oncogenic))
    #dev.off() 
  }
  
  OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/GSA/MSIGDBR/C7/"
  
  c7 = msigdbr(species = "Homo sapiens", category = "C7")
  ont = c7[,3]
  gene = c7[,8]
  c7_table = data.frame(ont,gene)
  colnames(c7_table) = c("ont","gene")
  immunologic_signatures <- enricher(entrezgene, TERM2GENE=c7_table, pvalueCutoff = 0.05)
  immunologic_signatures_table = data.frame(immunologic_signatures)
  write.csv(immunologic_signatures_table,paste0(OutPath,"immunologic_signatures_cluster",m,".csv"))
  #plots
  if (nrow(immunologic_signatures_table) > 0){ 
    pdf(file=paste0(OutPath,"immunologic_signatures_barplot_cluster_",m,".pdf"), width=8, height=8)            
    plot(barplot(immunologic_signatures))
    dev.off()
    #pdf(file=paste0(OutPath,"immunologic_signatures_emapplot_cluster_",m,".pdf"), width=8, height=8)  
    #plot(emapplot(immunologic_signatures))
    #dev.off() 
  }
}


############################################################################################
############################################################################################
############################################################################################
DefaultAssay(combined_seurat) <- "RNA"
#FeaturePlot  
#p2 = FeaturePlot(
#      object = combined_seurat,
#      features = c("SOX2"),
#      pt.size = 0.1,
#      max.cutoff = 'q95',
#      ncol = 1
#      )
#OR
#p2 = RidgePlot(
#      object = combined_seurat,
#      features = c("SOX2"),
#      ncol = 1
#      )
#
#Plot Genes from OPTION _1 AND _2
pdf(file="SC_RNA_interesting_genes_time_2.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("SOX2"),
      ncol = 1
      )
p3 = FeaturePlot(
      object = combined_seurat,
      features = c("SOX3"),
      ncol = 1
      )
p4 = FeaturePlot(
      object = combined_seurat,
      features = c("SOX4"),
      ncol = 1
      )
p5 = FeaturePlot(
      object = combined_seurat,
      features = c("SOX9"),
      ncol = 1
      )
p6 = FeaturePlot(
      object = combined_seurat,
      features = c("SOX11"),
      ncol = 1
      )
p7 = FeaturePlot(
      object = combined_seurat,
      features = c("DCX"),
      ncol = 1
      )
p8 = FeaturePlot(
      object = combined_seurat,
      features = c("NCAM1"),
      ncol = 1
      )
CombinePlots(list(p1, p2,p3,p4,p5,p6,p7,p8))
dev.off()
#Plot Genes from  OPTION _3
pdf(file="SC_RNA_interesting_genes_3_time.pdf", width=6, height=6)
features = c("SOX2","SOX3","SOX4","SOX9","SOX11","DCX","NCAM1")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()

pdf(file="SC_RNA_interesting_genes_FINAL.pdf", width=6, height=6)
features = c("SOX2","SOX3","SOX4","SOX11","DCX","SOX9")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()


#Plot Genes from  OPTION _4
pdf(file="SC_RNA_interesting_genes_4_umap.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("NES"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
p3 = FeaturePlot(
      object = combined_seurat,
      features = c("PAX6"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
CombinePlots(list(p1,p2,p3))
dev.off()
pdf(file="SC_RNA_interesting_genes_4_clusters_balls.pdf", width=6, height=6)
features = c("NES","PAX6")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()

#Plot Genes from OPTION _5
pdf(file="SC_RNA_interesting_genes_5_umap.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("PTPRZ1"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
p3 = FeaturePlot(
      object = combined_seurat,
      features = c("FABP7"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
p4 = FeaturePlot(
      object = combined_seurat,
      features = c("LGALS1"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
CombinePlots(list(p1,p2,p3,p4))
dev.off()
pdf(file="SC_RNA_interesting_genes_5_balls.pdf", width=6, height=6)
features = c("PTPRZ1","FABP7","LGALS1")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()

#Plot Genes from OPTION _6
pdf(file="SC_RNA_interesting_genes_6_umap.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("HES5"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
p3 = FeaturePlot(
      object = combined_seurat,
      features = c("PROM1"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
CombinePlots(list(p1,p2,p3))
dev.off()
pdf(file="SC_RNA_interesting_genes_6_balls.pdf", width=6, height=6)
features = c("HES5","PROM1")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()

#Plot Genes from  OPTION _7
pdf(file="SC_RNA_interesting_genes_7_umap.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("PROX1"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
p3 = FeaturePlot(
      object = combined_seurat,
      features = c("VIM"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
p4 = FeaturePlot(
      object = combined_seurat,
      features = c("DUSP6"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
CombinePlots(list(p1,p2,p3,p4))
dev.off()
pdf(file="SC_RNA_interesting_genes_7_clusters_balls.pdf", width=6, height=6)
features = c("PROX1","VIM","DUSP6")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()

pdf(file="SC_RNA_interesting_genes_8_umap.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("NXN"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
CombinePlots(list(p1,p2,p3))
dev.off()
pdf(file="SC_RNA_interesting_genes_8_balls.pdf", width=6, height=6)
features = c("NXN")
DotPlot(combined_seurat, features = features, cols = c("blue", "yellow")) + RotatedAxis()
dev.off()

#Plot cell_cycle proportions per cluster 
meta.data <- marrow@meta.data
pdf(file="CELLCYCLE_proportions_by_samples_time.pdf", width=18, height=6)
ggplot(meta.data, aes(x=protocol, fill=Phase)) + geom_bar(position = "fill")
dev.off()


###########################################
###########################################
###########################################
#Transiciones
#Heatmap de los logFC de los contrastes de controlVStratamiento
vector_days = c()
vector_days = append(vector_days,(rep("20_days", length(grep("20",(combined_seurat@meta.data$protocol))))))
vector_days = append(vector_days,(rep("13_days", length(grep("13",(combined_seurat@meta.data$protocol))))))
vector_days = append(vector_days,(rep("7_days", length(grep("7",(combined_seurat@meta.data$protocol))))))
vector_days = append(vector_days,(rep("Basal_days", length(grep("basal",(combined_seurat@meta.data$protocol))))))
combined_seurat@meta.data$vector_days = vector_days

#For not modify combined_seurat done in redoing with only neuro cells
combined_seurat_original = combined_seurat
#######
#REPEAT BUT ONLY WITH NEURONS FLOW CELLS 
neuro_clusters = c("0","1","2","3","4","6","7","8","13","15")
Idents(combined_seurat) <- "integrated_snn_res.0.6"
combined_seurat = subset(x = combined_seurat, idents = neuro_clusters))
#######

marker_genes = c()

Idents(combined_seurat) <- "vector_days"
for (n in 1:(length(levels(Idents(combined_seurat)))-1)) { # for each time_point 20/13/7
  cells_time_point_n = subset(x = combined_seurat, idents=c(levels(Idents(combined_seurat))[n]))
  for (m in neuro_clusters) { # for each cluster 0-17
    Idents(cells_time_point_n) <- "seurat_clusters"
    cells_cluster_n = subset(x = cells_time_point_n, idents=c(m))
    vector_condition = c()
    vector_condition = append(vector_condition,(rep("control", length(grep("control",(cells_cluster_n@meta.data$protocol))))))
    vector_condition = append(vector_condition,(rep("ba", length(grep("ba",(cells_cluster_n@meta.data$protocol))))))
    cells_cluster_n@meta.data$vector_condition = vector_condition
    if (length(cells_cluster_n@meta.data$vector_condition) > 0){
      if (table(cells_cluster_n@meta.data$vector_condition)[1] > 3 & table(cells_cluster_n@meta.data$vector_condition)[2] > 3){
        Idents(cells_cluster_n) <- "vector_condition"
        resulting.markers = FindMarkers(object = cells_cluster_n, ident.1 = "control", ident.2 = "ba", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.125)
        resulting.markers = subset(resulting.markers,(resulting.markers['p_val_adj'] <= 0.01 ))
        if (nrow(resulting.markers) > 0){
          marker_genes = append(marker_genes, rownames(resulting.markers))
        }
      }
    }
  }
}
all_genes_diff_on_all_contrasts = unique(marker_genes)
all_genes_diff_on_all_contrasts_final = all_genes_diff_on_all_contrasts

marker_table =c()
Idents(combined_seurat) <- "vector_days"
for (n in 1:(length(levels(Idents(combined_seurat)))-1)) { # for each time_point 20/13/7
  cells_time_point_n = subset(x = combined_seurat, idents=c(levels(Idents(combined_seurat))[n]))
  for (m in neuro_clusters) { # for each cluster 0-17
    Idents(cells_time_point_n) <- "seurat_clusters"
    cells_cluster_n = subset(x = cells_time_point_n, idents=c(m))
    vector_condition = c()
    vector_condition = append(vector_condition,(rep("control", length(grep("control",(cells_cluster_n@meta.data$protocol))))))
    vector_condition = append(vector_condition,(rep("ba", length(grep("ba",(cells_cluster_n@meta.data$protocol))))))
    cells_cluster_n@meta.data$vector_condition = vector_condition
    if (length(cells_cluster_n@meta.data$vector_condition) > 0){
      if (table(cells_cluster_n@meta.data$vector_condition)[1] > 3 & table(cells_cluster_n@meta.data$vector_condition)[2] > 3){
        Idents(cells_cluster_n) <- "vector_condition"
        resulting.markers = FindMarkers(object = cells_cluster_n, ident.1 = "control", ident.2 = "ba", only.pos = FALSE, min.pct = 0, logfc.threshold = 0, features = all_genes_diff_on_all_contrasts_final, return.thresh = 1.01)
        if (nrow(resulting.markers) > 0){
          marker_genes_logfc = resulting.markers$avg_log2FC
            if (length(marker_genes_logfc) < length(all_genes_diff_on_all_contrasts_final)){
              all_genes_diff_on_all_contrasts_final = all_genes_diff_on_all_contrasts_final[all_genes_diff_on_all_contrasts_final %in% rownames(resulting.markers)]
              resulting.markers = FindMarkers(object = cells_cluster_n, ident.1 = "control", ident.2 = "ba", only.pos = FALSE, min.pct = 0, logfc.threshold = 0, features = all_genes_diff_on_all_contrasts_final, return.thresh = 1.01)
              marker_table = append(marker_table, list(data.frame(genes=all_genes_diff_on_all_contrasts_final, logfc=marker_genes_logfc)))
              print(paste0("Bucle_tiempo_",n,"_cluster_",m))
              print(dim(data.frame(genes=all_genes_diff_on_all_contrasts_final, logfc=marker_genes_logfc))) 
            } else {
              marker_table = append(marker_table, list(data.frame(genes=all_genes_diff_on_all_contrasts_final, logfc=marker_genes_logfc)))
              print(paste0("Bucle_tiempo_",n,"_cluster_",m))
              print(dim(data.frame(genes=all_genes_diff_on_all_contrasts_final, logfc=marker_genes_logfc))) 
          }
        }
      }
    }
  }
}
true_list_of_genes = as.character(marker_table[[length(marker_table)]]$genes)

for (n in 1:(length(marker_table))) {
  marker_table[[n]] = marker_table[[n]][marker_table[[n]]$genes %in% true_list_of_genes, ]
}

#naw we have all contracst with same genes, lets create a merged table with all of them
tabla_final <- data.frame(genes = true_list_of_genes)
tiempos = c("Day20","Day13","Day7")
for (n in 1:(3)) {
  for (m in 1:(length(marker_table)/3)) {
    new = marker_table[[m*n]]$logfc
    tabla_final[ , ncol(tabla_final) + 1] = new
    colnames(tabla_final)[ncol(tabla_final)] <- paste0(tiempos[n] ,"_cluster_", m-1) #neuro_clusters
  }
}

#############
#Plot tabla_final in a heatmap
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize) # for the colorRamp2() function

tabla_final
rownames(tabla_final)= tabla_final$genes
tabla_final$genes = NULL

reorder_names = c()
tiempos = c("Day7","Day13","Day20")
for (n in 1:(3)) {
  for (m in 1:(length(marker_table)/3)) {
    reorder_names = append(reorder_names, print(paste0(tiempos[n] ,"_cluster_", m-1))) #neuro_clusters
  }
}
tabla_final = tabla_final[,reorder_names]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Day7", 16)
tiempos=append(tiempos,rep("Day13", 16))
tiempos=append(tiempos,rep("Day20", 16))
#tiempos=c("T2","T4","T4","T3","T2","T2","T4","T3","T4","T1","T1","T1","T2","T3","T2","T4")
#tiempos =c("T1","T3","T3","T2","T1","T3","T2","T3","T2","T2")
meta_data = data.frame(Control_vs_Treatment=colnames(tabla_final), Time=tiempos)
rownames(meta_data)=meta_data$Control_vs_Treatment
meta_data$Control_vs_Treatment = NULL
# order of annotations/colors are defined here
#ordered_meta_data <- meta_data[order(meta_data$numeros), ]

# OPTIONAL: YOU CAN PICK COLORS FOR EACH LEVEL OF ANNOTATION
# HERE I PROVIDE COLORS FOR ONLY TWO OF THE THREE LEVELS, COLORS FOR THE
# REMAINING LEVEL IS TAKEN CARE OF BY THE PACKAGE
annotation_colors <- list("Time"= c(     "Day7" = "blue",
                                         "Day13" = "orange",
                                         "Day20" = "green")
                          )
ha = HeatmapAnnotation(df = meta_data,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# Expression data
genes_to_use <- true_list_of_genes

#seurat_object <- ScaleData(combined_seurat, genes.use = genes_to_use)

my_data <- tabla_final #seurat_object@assays$RNA@scale.data

# COLUMN ORDER OF THE EXPRESSION DATA SHOULD MATCH THE ROW ORDER OF THE
# ANNOTATION TABLE
my_data <- my_data[, rownames(meta_data)]
my_data_max=max(my_data)
my_data_min=min(my_data)
# Heatmap
col_fun = colorRamp2(c(my_data_min/2, 0, my_data_max/2), c("blue", "white", "red"))
pdf(file = "example_heatmap.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Control vs Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha
)
dev.off()

#Binarized version of table and plot
my_data_binarized = my_data
my_data_binarized[my_data_binarized < 0] <- -1
my_data_binarized[my_data_binarized > 0] <- 1
# Heatmap
col_fun = colorRamp2(c(-1, 1), c("lightyellow", "purple"))
pdf(file = "example_heatmap_2.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data_binarized,
  heatmap_legend_param = list(title = "Control vs Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha
)
dev.off()




#################################
#################################
library(qvalue)
library(jaccard)
#Binarize my_data into 0/1 been 1 its overexpressed
#Jaccard_distance <- 1 - jaccard(Binary_A,Binary_B) #Jacard distance
jacard_similarity = matrix(, nrow = 16, ncol = 16)
colnames(jacard_similarity) = colnames(my_data_binarized)
rownames(jacard_similarity) = colnames(my_data_binarized)
for (n in 1:(ncol(my_data_binarized))) { 
  for (m in 1:(ncol(my_data_binarized))) { 
    Binary_A <- my_data_binarized[,n]
    Binary_B <- my_data_binarized[,m]
    print(jaccard(Binary_A,Binary_B))
    pair_similarity=jaccard(Binary_A,Binary_B)
    jacard_similarity[n,m] = pair_similarity
  }
}
# Heatmap
col_fun = colorRamp2(c(0, 1), c("lightyellow", "purple"))
pdf(file = "Jaccard_similarity_downregulated_heatmap.pdf",
    width = 10,
    height = 10)
Heatmap(
  jacard_similarity,
  heatmap_legend_param = list(title = "Jaccard similarity - Ctr vs T"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  #top_annotation = ha
)
dev.off()

#Corr matrix
library("Hmisc")
corr_mat=rcorr(as.matrix(my_data), type = c("spearman"))
# Heatmap
col_fun = colorRamp2(c(0, 1), c("lightyellow", "purple"))
pdf(file = "DIFFexpression_correlation_heatmap.pdf",
    width = 10,
    height = 10)
Heatmap(
  corr_mat$r,
  heatmap_legend_param = list(title = "Spearman correlation - Ctr vs T"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  #top_annotation = ha
)
dev.off()


########
########
#DExtract contrasts genes from the 4 big groups we are interested_in
level_1 = c(0,9,4,1)
level_2 = c(3,5,7,10,11)
level_3 = c(2,13,6,15,8,12)
level_4 = c(14)

0, 9 y 4
3, 5, 7, 10 y 11
2, 13, 6, 15, 8 (y pensar�a en sacar el 12 si cambia mucho haci�ndolo o incluy�ndolo)
1
y 14


# Adding column based on other column:
combined_seurat@meta.data$velocity_level = with(combined_seurat@meta.data, ifelse(seurat_clusters == c(0), 'level_1',
                                                                           ifelse(seurat_clusters == c(9), 'level_1',
                                                                           ifelse(seurat_clusters == c(4), 'level_1',
                                                                           ifelse(seurat_clusters == c(1), 'level_1',
                                                                           ifelse(seurat_clusters == c(3), 'level_2',
                                                                           ifelse(seurat_clusters == c(5), 'level_2',
                                                                           ifelse(seurat_clusters == c(7), 'level_2',
                                                                           ifelse(seurat_clusters == c(10), 'level_2',
                                                                           ifelse(seurat_clusters == c(11), 'level_2',
                                                                           ifelse(seurat_clusters == c(2), 'level_3',
                                                                           ifelse(seurat_clusters == c(13), 'level_3',
                                                                           ifelse(seurat_clusters == c(6), 'level_3',
                                                                           ifelse(seurat_clusters == c(15), 'level_3',
                                                                           ifelse(seurat_clusters == c(8), 'level_3',
                                                                           ifelse(seurat_clusters == c(12), 'level_3',
                                                                           ifelse(seurat_clusters == c(14), 'level_4', 'ERROR')))))))))))))))))

combined_seurat@meta.data$velocity_level = paste(combined_seurat@meta.data$velocity_level, combined_seurat@meta.data$protocol_2, sep="_")
Idents(combined_seurat) <- "velocity_level"

#combined_seurat <- subset(combined_seurat, subset=integrated_snn_res.0.6 %in% c(1, 2, 6, 8, 13, 15))
#Idents(combined_seurat) <- "protocol_2"
#markers_for_Dani = FindMarkers(object = combined_seurat, ident.1 = "ba", ident.2 = "control", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)

#markers_for_Dani_final_subset = (subset(markers_for_Dani_final,(markers_for_Dani_final['p_val_adj'] <= 0.01 )))


#Do contrasts
#level_1
velocity_level_markers = FindMarkers(object = combined_seurat, ident.1 = "level_1_control", ident.2 = "level_1_ba", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(velocity_level_markers, file = paste0("Velocity_level_1_markers_clusters_0_9_4_1_.csv"))
#level_2
velocity_level_markers = FindMarkers(object = combined_seurat, ident.1 = "level_2_control", ident.2 = "level_2_ba", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(velocity_level_markers, file = paste0("Velocity_level_2_markers_clusters_3_5_7_10_11_.csv"))
#level_3
velocity_level_markers = FindMarkers(object = combined_seurat, ident.1 = "level_3_control", ident.2 = "level_3_ba", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(velocity_level_markers, file = paste0("Velocity_level_3_markers_clusters_2_13_6_15_.csv"))
#level_4
velocity_level_markers = FindMarkers(object = combined_seurat, ident.1 = "level_4_control", ident.2 = "level_4_ba", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
write.csv(velocity_level_markers, file = paste0("Velocity_level_4_markers_clusters_14_.csv"))


#GSEA for velocity_levels information ranking base on control vs BA contrast
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/velocity_levels_markers/GSEA/GO/"

#Change names to entrez code
name_changer = function(especie_database, especie_name, x){
  ensembl <- useMart(especie_database, dataset = especie_name)
  annotation = getBM(attributes=c("hgnc_symbol","entrezgene_id"), values=x, mart=ensembl)
  return(annotation)
}

#Create the geneList
geneList <- velocity_level_markers$avg_log2FC
names(geneList) <- as.character(rownames(velocity_level_markers))
geneList = sort(geneList, decreasing = TRUE)

annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
names(geneList) = as.character(entrezgene)
geneList = geneList[!is.na(names(geneList))]

#ahora estoy con el lvl2

ont = c("CC", "MF", "BP")
for (n in 1:length(ont)) {
  ego <- gseGO(geneList    = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = ont[n],
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
  ego_table <- as.data.frame(ego) 
  write.csv(ego_table,paste0(OutPath,"ego_",ont[n],"_cluster_",m,".csv"))
  #plots
  if (nrow(ego_table) > 0){
    #ego <- pairwise_termsim(ego) 
    pdf(file=paste0(OutPath,"ego__barplot_",ont[n],"_cluster_",m,".pdf"), width=8, height=8)  
    plot(barplot(ego, showCategory=20))  
    dev.off()
    pdf(file=paste0(OutPath,"ego_cnetplot_",ont[n],"_cluster_",m,".pdf"), width=8, height=8)  
    plot(cnetplot(ego, categorySize="pvalue", foldChange=geneList))
    dev.off()  
  }  
}  


############################################################################################################
##################### NEW VERSION OF DE ANALYSIS
############################################################################################################

#Only control cells 
Idents(combined_seurat) <- "protocol_2"
only_control_cells = subset(x = combined_seurat, idents=c("control"))
Idents(only_control_cells) <- "integrated_snn_res.0.6"
only_control_cells_subset = subset(x = only_control_cells, idents=c("0","1","2","3","4","6","7","8","13","15"))

############################################################################################################
###########################################
#Transiciones
#Heatmap de los logFC de los contrastes de controlVStratamiento
############################################################################################################

#Contrasts para obtener los genes diferenciales en los diferentes clusters
Idents(only_control_cells_subset) = "integrated_snn_res.0.6"
control_all_markers = FindAllMarkers(object = only_control_cells_subset, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.125)
#Get all agregated unique gene names from each cluster 
de_genes_any_cluster = c()
for (n in levels(control_all_markers$cluster)){
  print(n)
  de_genes_any_cluster = append(de_genes_any_cluster, control_all_markers[control_all_markers$cluster == 0,]$gene[control_all_markers[control_all_markers$cluster == 0,]$gene %in% control_all_markers[control_all_markers$cluster == n,]$gene]) 
}
#total of 2275
de_genes_any_cluster = unique(de_genes_any_cluster)
#Realizamos un an�lisis diferencial solo para los genes diferenciales en todos los clusters para generar la tabla para calcualr correlaciones
control_all_markers = FindAllMarkers(object = only_control_cells_subset, only.pos = FALSE, min.pct = 0, logfc.threshold = 0, features = de_genes_any_cluster, return.thresh = 1.01)

#Generate matrix with fold changes
avg_fold_2fc_matrix = data.frame(matrix(nrow=2275, ncol=0))
for (n in levels(control_all_markers$cluster)){
  avg_fold_2fc_matrix[,n] = control_all_markers[control_all_markers$cluster == n,]$avg_log2FC
}

my_data = avg_fold_2fc_matrix
my_data_genes = unique(control_all_markers$gene)
#Cor plot for each cluster 
library("Hmisc")
corr_mat=rcorr(as.matrix(my_data), type = c("spearman"))
# Heatmap
col_fun = colorRamp2(c(0, 1), c("lightyellow", "purple"))
pdf(file = "Markers_correlation_heatmap_using_top_marker_genes.pdf",
    width = 10,
    height = 10)
Heatmap(
  corr_mat$r,
  heatmap_legend_param = list(title = "Spearman correlation - Markers"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  #top_annotation = ha
)
dev.off()


save.image("08_06_2022.RData")


###########################################
#DESCRIPTION QC
#Over only control cells
###########################################
meta.data <- only_control_cells_subset@meta.data

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(count = n())
levels(counts$protocol) = c("data_rna_control_7_days","data_rna_control_13_days","data_rna_control_20_days")
counts$protocol<-factor(counts$protocol,levels=levels(counts$protocol)[1:3]) 
pdf(file="Number_Cells_distribution.pdf", width=6, height=6)
ggplot(counts, aes(seurat_clusters, count, fill = protocol)) +
  geom_bar(stat = 'identity', position = "fill")
dev.off()

counts$count = (counts$count / nrow(meta.data))*100
pdf(file="Number_Cells_percentage_distribution.pdf", width=6, height=6)
ggplot(counts, aes(seurat_clusters, count, fill = protocol)) +
  geom_bar(stat = 'identity', position = "fill")
dev.off()

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(nCount_RNA = mean(nCount_RNA))
levels(counts$protocol) = c("data_rna_control_7_days","data_rna_control_13_days","data_rna_control_20_days")
counts$protocol<-factor(counts$protocol,levels=levels(counts$protocol)[1:3]) 
pdf(file="Counts_rna_mean_distribution.pdf", width=6, height=6)
ggplot(counts, aes(seurat_clusters, nCount_RNA, fill = protocol)) +
  geom_bar(stat = 'identity', position = "fill")
dev.off()

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(percent_mt = mean(percent.mt))
levels(counts$protocol) = c("data_rna_control_7_days","data_rna_control_13_days","data_rna_control_20_days")
counts$protocol<-factor(counts$protocol,levels=levels(counts$protocol)[1:3]) 
pdf(file="MT_percentage_mean_distribution.pdf", width=6, height=6)
ggplot(counts, aes(seurat_clusters, percent_mt, fill = protocol)) +
  geom_bar(stat = 'identity', position = "fill")
dev.off()

counts <- group_by(meta.data, protocol, seurat_clusters) %>% summarise(nFeature_RNA = mean(nFeature_RNA))
levels(counts$protocol) = c("data_rna_control_7_days","data_rna_control_13_days","data_rna_control_20_days")
counts$protocol<-factor(counts$protocol,levels=levels(counts$protocol)[1:3]) 
pdf(file="Number_Features_mean_distribution.pdf", width=6, height=6)
ggplot(counts, aes(seurat_clusters, nFeature_RNA, fill = protocol)) +
  geom_bar(stat = 'identity', position = "fill")
dev.off()

#Ver como cambia la expresion en base a los knn vecinos para un gen con respecto a su expresion en el cluster
only_control_cells_subset_cluster_0 = subset(x = only_control_cells_subset, idents=c("1"))

#original
mean(only_control_cells_subset_cluster_0@assays$RNA@data['DCX',])

#KNN based
#20 neighbour cells for each cell
near_data <- get.knn(only_control_cells_subset_cluster_0@reductions$umap@cell.embeddings[, 1:2], k = 20)
near_data=data.frame(near_data)
rownames(near_data) = rownames(only_control_cells_subset_cluster_0@reductions$umap@cell.embeddings)
near_data_index = near_data[,1:20]

#List of vectors with near 20 cells names for each cell
near_cells_for_each_cell = c()
for (m in 1:nrow(near_data_index)){
  near_cells_for_each_cell = append(near_cells_for_each_cell, list(rownames(near_data_index)[unlist(near_data_index[m,])]))
}

#Original RNA data
near_cells_for_each_cell_expression_each_protein = c()
near_cells_for_each_cell_expression_all = c()

#para todos los genes
for(n in 1:23350){
  near_cells_for_each_cell_expression_each_protein = c()
  for (m in 1:nrow(near_data_index)){
    near_cells_for_each_cell_expression_each_protein = append(near_cells_for_each_cell_expression_each_protein, mean(only_control_cells_subset_cluster_0@assays$RNA@data[n,][near_cells_for_each_cell[[m]]]))
  }
  near_cells_for_each_cell_expression_all = append(near_cells_for_each_cell_expression_all, list(near_cells_for_each_cell_expression_each_protein))
}
near_cells_for_each_cell_expression_all_Original = near_cells_for_each_cell_expression_all

neighbours_vs_full_cells = c()
for(n in 1:17220){
  neighbours_vs_full_cells = append(neighbours_vs_full_cells, cor(only_control_cells_subset_cluster_0@assays$RNA@data[n,], near_cells_for_each_cell_expression_all[[n]], method = "pearson"))
}
#for the position of DCX check the expresion for using only the k=20 nn cells average expression base on umap space

#for a subset of genes
#match("DCX",rownames(only_control_cells_subset_cluster_0@assays$RNA@data))
lista_interesting_genes = c(8867, 11585, 8969, 4471, 10936)
for(n in lista_interesting_genes){
  near_cells_for_each_cell_expression_each_protein = c()
  for (m in 1:nrow(near_data_index)){
    near_cells_for_each_cell_expression_each_protein = append(near_cells_for_each_cell_expression_each_protein, mean(only_control_cells_subset_cluster_0@assays$RNA@data[n,][near_cells_for_each_cell[[m]]]))
  }
  near_cells_for_each_cell_expression_all = append(near_cells_for_each_cell_expression_all, list(near_cells_for_each_cell_expression_each_protein))
}

#umap
neighbours_vs_full_cells = c()
means_origianl_data = c()
means_neighbours = c()
for(n in 1:length(lista_interesting_genes)){
  neighbours_vs_full_cells = append(neighbours_vs_full_cells, cor(only_control_cells_subset_cluster_0@assays$RNA@data[lista_interesting_genes[n],], near_cells_for_each_cell_expression_all[[n]], method = "pearson"))
  means_origianl_data = append(means_origianl_data, mean(only_control_cells_subset_cluster_0@assays$RNA@data[lista_interesting_genes[n],]))
  means_neighbours = append(means_neighbours, mean(near_cells_for_each_cell_expression_all[[n]]))
}

#pca
neighbours_vs_full_cells_2 = c()
means_origianl_data_2 = c()
means_neighbours_2 = c()
for(n in 1:length(lista_interesting_genes)){
  neighbours_vs_full_cells_2 = append(neighbours_vs_full_cells_2, cor(only_control_cells_subset_cluster_0@assays$RNA@data[lista_interesting_genes[n],], near_cells_for_each_cell_expression_all[[n]], method = "pearson"))
  means_origianl_data_2 = append(means_origianl_data_2, mean(only_control_cells_subset_cluster_0@assays$RNA@data[lista_interesting_genes[n],]))
  means_neighbours_2 = append(means_neighbours_2, mean(near_cells_for_each_cell_expression_all[[n]]))
}

test2 <- rbind(neighbours_vs_full_cells, neighbours_vs_full_cells_2)
pdf(file=paste0("Genes_pearson_umap_pca_c1",".pdf"), width=12, height=6)
barplot(test2,
main = "kNN 20",
col = c("#60be00","#bebe00"),
beside = TRUE,
names.arg = rownames(only_control_cells_subset_cluster_0@assays$RNA@data)[c(8867,11585,8969,4471,10936)], las=2, cex.names=0.75, legend = TRUE, args.legend = list(bty = "n", x = "top", ncol = 2, inset = -0.075, legend=c("UMAP_nn","PCA_nn")))
dev.off()

test3 <- rbind(means_origianl_data, means_neighbours, means_neighbours_2)
pdf(file=paste0("Genes_mean_expression_umap_pca_c1",".pdf"), width=12, height=6)
barplot(test3,
main = "kNN 20",
col = c("#60be00","#bebe00","#00bebe"),
beside = TRUE,
names.arg = rownames(only_control_cells_subset_cluster_0@assays$RNA@data)[c(8867,11585,8969,4471,10936)], las=2, cex.names=0.75, legend = TRUE, args.legend = list(bty = "n", x = "top", ncol = 3, inset = -0.075, legend=c("Original_mean_expr", "UMAP_nn_mean_exp", "PCA_nn_mean_exp")))
dev.off()




#NOT RUNNED
only_control_cells_subset_graph= FindVariableFeatures(only_control_cells_subset)
only_control_cells_subset_graph = ScaleData(object = only_control_cells_subset_graph, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"))
only_control_cells_subset_graph = RunPCA(object = only_control_cells_subset_graph, npcs = 30,verbose = FALSE)
only_control_cells_subset_graph = RunTSNE(object = only_control_cells_subset_graph, reduction = "pca", dims = 1:cc_number)
pdf(file="Graph_tsne_only_control_cells_2.pdf", width=6, height=6)
only_control_cells_subset_graph = FindNeighbors(object = only_control_cells_subset_graph, dims = 1:cc_number, do.plot = TRUE, k.param = 15, reduction = "pca")
dev.off()

#How to add colors to graph plot?
as.matrix(ceiling(as.matrix(only_control_cells_subset_graph@graphs$integrated_snn)))















#####################################################################################################################
#
#####################################################################################################################
# STEP 1 - NEURONS DEVELOPMENT SUBSET CONTROL CELLS - CONTRASTS
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/")
load("02_08_2022.RData")
#Generate new metadata grouping for the overall 5 groups of trajectory
only_control_cells_subset@meta.data$group_A = only_control_cells_subset@meta.data$integrated_snn_res.0.6
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "4"] <- "A"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "0"] <- "B"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "7"] <- "C"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "3"] <- "D"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "13"] <- "E"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "15"] <- "E"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "1"] <- "E"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "8"] <- "E"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "2"] <- "E"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "6"] <- "E"

Idents(only_control_cells_subset) = "group_A"

potential_markers_BA = FindMarkers(object = only_control_cells_subset, ident.1 = "B", ident.2 = "A", only.pos = FALSE) 
potential_markers_CB = FindMarkers(object = only_control_cells_subset, ident.1 = "C", ident.2 = "B", only.pos = FALSE)
potential_markers_DC = FindMarkers(object = only_control_cells_subset, ident.1 = "D", ident.2 = "C", only.pos = FALSE)
potential_markers_ED = FindMarkers(object = only_control_cells_subset, ident.1 = "E", ident.2 = "D", only.pos = FALSE)

potential_markers_BA = potential_markers_BA[abs(potential_markers_BA$p_val_adj) < 0.05, ]
potential_markers_CB = potential_markers_CB[abs(potential_markers_CB$p_val_adj) < 0.05, ]
potential_markers_DC = potential_markers_DC[abs(potential_markers_DC$p_val_adj) < 0.05, ]
potential_markers_ED = potential_markers_ED[abs(potential_markers_ED$p_val_adj) < 0.05, ]

#Overall developmental A-B-C-D-E
markers_BA = FindMarkers(object = only_control_cells_subset, ident.1 = "B", ident.2 = "A", only.pos = FALSE, min.pct = 0, logfc.threshold = 0) #0.25 0.125
markers_CB = FindMarkers(object = only_control_cells_subset, ident.1 = "C", ident.2 = "B", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_DC = FindMarkers(object = only_control_cells_subset, ident.1 = "D", ident.2 = "C", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_ED = FindMarkers(object = only_control_cells_subset, ident.1 = "E", ident.2 = "D", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)

#ALL to ALL in E for all clusters
Idents(only_control_cells_subset) = "integrated_snn_res.0.6"
#8,1,6,15,2,13
#8
markers_E_1_8 = FindMarkers(object = only_control_cells_subset, ident.1 = "1", ident.2 = "8", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_6_8 = FindMarkers(object = only_control_cells_subset, ident.1 = "6", ident.2 = "8", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_15_8 = FindMarkers(object = only_control_cells_subset, ident.1 = "15", ident.2 = "8", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_2_8 = FindMarkers(object = only_control_cells_subset, ident.1 = "2", ident.2 = "8", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_13_8 = FindMarkers(object = only_control_cells_subset, ident.1 = "13", ident.2 = "8", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
#1
markers_E_8_1 = FindMarkers(object = only_control_cells_subset, ident.1 = "8", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_6_1 = FindMarkers(object = only_control_cells_subset, ident.1 = "6", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_15_1 = FindMarkers(object = only_control_cells_subset, ident.1 = "15", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_2_1 = FindMarkers(object = only_control_cells_subset, ident.1 = "2", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_13_1 = FindMarkers(object = only_control_cells_subset, ident.1 = "13", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
#6
markers_E_1_6 = FindMarkers(object = only_control_cells_subset, ident.1 = "1", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_8_6 = FindMarkers(object = only_control_cells_subset, ident.1 = "8", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_15_6 = FindMarkers(object = only_control_cells_subset, ident.1 = "15", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_2_6 = FindMarkers(object = only_control_cells_subset, ident.1 = "2", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_13_6 = FindMarkers(object = only_control_cells_subset, ident.1 = "13", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
#15
markers_E_1_15 = FindMarkers(object = only_control_cells_subset, ident.1 = "1", ident.2 = "15", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_6_15 = FindMarkers(object = only_control_cells_subset, ident.1 = "6", ident.2 = "15", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_8_15 = FindMarkers(object = only_control_cells_subset, ident.1 = "8", ident.2 = "15", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_2_15 = FindMarkers(object = only_control_cells_subset, ident.1 = "2", ident.2 = "15", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_13_15 = FindMarkers(object = only_control_cells_subset, ident.1 = "13", ident.2 = "15", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
#2
markers_E_1_2 = FindMarkers(object = only_control_cells_subset, ident.1 = "1", ident.2 = "2", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_6_2 = FindMarkers(object = only_control_cells_subset, ident.1 = "6", ident.2 = "2", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_15_2 = FindMarkers(object = only_control_cells_subset, ident.1 = "15", ident.2 = "2", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_8_2 = FindMarkers(object = only_control_cells_subset, ident.1 = "8", ident.2 = "2", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_13_2 = FindMarkers(object = only_control_cells_subset, ident.1 = "13", ident.2 = "2", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
#13
markers_E_1_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "1", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_6_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "6", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_15_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "15", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_2_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "2", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_8_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "8", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)

list_of_contrasts = list(markers_BA,markers_CB,markers_DC,markers_ED,
                      markers_E_1_8,markers_E_6_8,markers_E_15_8,markers_E_2_8,markers_E_13_8,
                      markers_E_8_1,markers_E_6_1,markers_E_15_1,markers_E_2_1,markers_E_13_1,
                      markers_E_1_6,markers_E_8_6,markers_E_15_6,markers_E_2_6,markers_E_13_6,
                      markers_E_1_15,markers_E_6_15,markers_E_8_15,markers_E_2_15,markers_E_13_15,
                      markers_E_1_2,markers_E_6_2,markers_E_15_2,markers_E_8_2,markers_E_13_2,
                      markers_E_1_13,markers_E_6_13,markers_E_15_13,markers_E_2_13,markers_E_8_13)
                      
names(list_of_contrasts) = c("markers_BA","markers_CB","markers_DC","markers_ED",
                      "markers_E_1_8","markers_E_6_8","markers_E_15_8","markers_E_2_8","markers_E_13_8",
                      "markers_E_8_1","markers_E_6_1","markers_E_15_1","markers_E_2_1","markers_E_13_1",
                      "markers_E_1_6","markers_E_8_6","markers_E_15_6","markers_E_2_6","markers_E_13_6",
                      "markers_E_1_15","markers_E_6_15","markers_E_8_15","markers_E_2_15","markers_E_13_15",
                      "markers_E_1_2","markers_E_6_2","markers_E_15_2","markers_E_8_2","markers_E_13_2",
                      "markers_E_1_13","markers_E_6_13","markers_E_15_13","markers_E_2_13","markers_E_8_13")
######################################################################################################################################################
#CLUSTERPROFILER
#GSEA for velocity_levels information ranking base on control vs BA contrast
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/GO/"

#Change names to entrez code
name_changer = function(especie_database, especie_name, x){
  ensembl <- useMart(especie_database, dataset = especie_name)
  annotation = getBM(attributes=c("hgnc_symbol","entrezgene_id"), values=x, mart=ensembl)
  return(annotation)
}

library("enrichplot")
library("clusterProfiler")
library("org.Hs.eg.db")
list_of_gsea = list()
for (n in 1:length(list_of_contrasts)){
  print(n)
  #Create the geneList
  geneList <- list_of_contrasts[[n]]$avg_log2FC
  names(geneList) <- as.character(rownames(list_of_contrasts[[n]]))
  geneList = sort(geneList, decreasing = TRUE)
   
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]
  
  #GSEA
  ego <- gseGO(geneList    = geneList,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 100,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = FALSE)
  ego_table <- as.data.frame(ego) 
  write.csv(ego_table,paste0(OutPath,"gsea_","BP","_contrast_",names(list_of_contrasts[n]),".csv"))
  list_of_gsea = append(list_of_gsea, list(ego_table))
  #plots
  if (nrow(ego_table) > 1){
    pdf(file=paste0(OutPath,"gsea_dotplot_","BP","_contrast_",names(list_of_contrasts[n]),".pdf"), width=8, height=12)  
    plot(dotplot(ego, showCategory=30))
    dev.off()
    pdf(file=paste0(OutPath,"gsea_gseaPlot_","BP","_contrast_",names(list_of_contrasts[n]),".pdf"), width=12, height=12)  
    plot(gseaplot2(ego, geneSetID = 1:3))
    dev.off()    
  }
}

######################################################################################################################################################


######################################################################################################################################################
#GeneSetCluster summarization

##### GeneSetClust load
library(GeneSetCluster)
library(ggplot2)
set.seed(1234567)
##set work directory
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/GeneSetCluster_step2/")

names(list_of_gsea) = c("markers_BA","markers_CB","markers_DC","markers_ED",
                      "markers_E_1_8","markers_E_6_8","markers_E_15_8","markers_E_2_8","markers_E_13_8",
                      "markers_E_8_1","markers_E_6_1","markers_E_15_1","markers_E_2_1","markers_E_13_1",
                      "markers_E_1_6","markers_E_8_6","markers_E_15_6","markers_E_2_6","markers_E_13_6",
                      "markers_E_1_15","markers_E_6_15","markers_E_8_15","markers_E_2_15","markers_E_13_15",
                      "markers_E_1_2","markers_E_6_2","markers_E_15_2","markers_E_8_2","markers_E_13_2",
                      "markers_E_1_13","markers_E_6_13","markers_E_15_13","markers_E_2_13","markers_E_8_13")

# STEP 1: with all samples included
pathways_markers_BA = list_of_gsea[["markers_BA"]]
pathways_markers_CB = list_of_gsea[["markers_CB"]]
pathways_markers_DC = list_of_gsea[["markers_DC"]]
pathways_markers_ED = list_of_gsea[["markers_ED"]]

pathways_markers_E_1_8 = list_of_gsea[["markers_E_1_8"]]
pathways_markers_E_6_8 = list_of_gsea[["markers_E_6_8"]]
pathways_markers_E_15_8 = list_of_gsea[["markers_E_15_8"]]
pathways_markers_E_2_8 = list_of_gsea[["markers_E_2_8"]]
pathways_markers_E_13_8 = list_of_gsea[["markers_E_13_8"]]

pathways_markers_E_8_1 = list_of_gsea[["markers_E_8_1"]]
pathways_markers_E_6_1 = list_of_gsea[["markers_E_6_1"]]
pathways_markers_E_15_1 = list_of_gsea[["markers_E_15_1"]]
pathways_markers_E_2_1 = list_of_gsea[["markers_E_2_1"]]
pathways_markers_E_13_1 = list_of_gsea[["markers_E_13_1"]]

pathways_markers_E_1_6 = list_of_gsea[["markers_E_1_6"]]
pathways_markers_E_8_6 = list_of_gsea[["markers_E_8_6"]]
pathways_markers_E_15_6 = list_of_gsea[["markers_E_15_6"]]
pathways_markers_E_2_6 = list_of_gsea[["markers_E_2_6"]]
pathways_markers_E_13_6 = list_of_gsea[["markers_E_13_6"]]

pathways_markers_E_1_15 = list_of_gsea[["markers_E_1_15"]]
pathways_markers_E_6_15 = list_of_gsea[["markers_E_6_15"]]
pathways_markers_E_8_15 = list_of_gsea[["markers_E_8_15"]]
pathways_markers_E_2_15 = list_of_gsea[["markers_E_2_15"]]
pathways_markers_E_13_15 = list_of_gsea[["markers_E_13_15"]]

pathways_markers_E_1_2 = list_of_gsea[["markers_E_1_2"]]
pathways_markers_E_6_2 = list_of_gsea[["markers_E_6_2"]]
pathways_markers_E_15_2 = list_of_gsea[["markers_E_15_2"]]
pathways_markers_E_8_2 = list_of_gsea[["markers_E_8_2"]]
pathways_markers_E_13_2 = list_of_gsea[["markers_E_13_2"]]

pathways_markers_E_1_13 = list_of_gsea[["markers_E_1_13"]]
pathways_markers_E_6_13 = list_of_gsea[["markers_E_6_13"]]
pathways_markers_E_15_13 = list_of_gsea[["markers_E_15_13"]]
pathways_markers_E_2_13 = list_of_gsea[["markers_E_2_13"]]
pathways_markers_E_8_13 = list_of_gsea[["markers_E_8_13"]]

#get significant pathways

pathways_markers_BA_sig<-pathways_markers_BA[pathways_markers_BA$p.adjust<0.1,]
pathways_markers_CB_sig<-pathways_markers_CB[pathways_markers_CB$p.adjust<0.1,]
pathways_markers_DC_sig<-pathways_markers_DC[pathways_markers_DC$p.adjust<0.1,]
pathways_markers_ED_sig<-pathways_markers_ED[pathways_markers_ED$p.adjust<0.1,]

pathways_markers_E_1_8_sig<-pathways_markers_E_1_8[pathways_markers_E_1_8$p.adjust<0.1,]
pathways_markers_E_6_8_sig<-pathways_markers_E_6_8[pathways_markers_E_6_8$p.adjust<0.1,]
pathways_markers_E_15_8_sig<-pathways_markers_E_15_8[pathways_markers_E_15_8$p.adjust<0.1,]
pathways_markers_E_2_8_sig<-pathways_markers_E_2_8[pathways_markers_E_2_8$p.adjust<0.1,]
pathways_markers_E_13_8_sig<-pathways_markers_E_13_8[pathways_markers_E_13_8$p.adjust<0.1,]

pathways_markers_E_8_1_sig<-pathways_markers_E_8_1[pathways_markers_E_8_1$p.adjust<0.1,]
pathways_markers_E_6_1_sig<-pathways_markers_E_6_1[pathways_markers_E_6_1$p.adjust<0.1,]
pathways_markers_E_15_1_sig<-pathways_markers_E_15_1[pathways_markers_E_15_1$p.adjust<0.1,]
pathways_markers_E_2_1_sig<-pathways_markers_E_2_1[pathways_markers_E_2_1$p.adjust<0.1,]
pathways_markers_E_13_1_sig<-pathways_markers_E_13_1[pathways_markers_E_13_1$p.adjust<0.1,]

pathways_markers_E_1_6_sig<-pathways_markers_E_1_6[pathways_markers_E_1_6$p.adjust<0.1,]
pathways_markers_E_8_6_sig<-pathways_markers_E_8_6[pathways_markers_E_8_6$p.adjust<0.1,]
pathways_markers_E_15_6_sig<-pathways_markers_E_15_6[pathways_markers_E_15_6$p.adjust<0.1,]
pathways_markers_E_2_6_sig<-pathways_markers_E_2_6[pathways_markers_E_2_6$p.adjust<0.1,]
pathways_markers_E_13_6_sig<-pathways_markers_E_13_6[pathways_markers_E_13_6$p.adjust<0.1,]

pathways_markers_E_1_15_sig<-pathways_markers_E_1_15[pathways_markers_E_1_15$p.adjust<0.1,]
pathways_markers_E_6_15_sig<-pathways_markers_E_6_15[pathways_markers_E_6_15$p.adjust<0.1,]
pathways_markers_E_8_15_sig<-pathways_markers_E_8_15[pathways_markers_E_8_15$p.adjust<0.1,]
pathways_markers_E_2_15_sig<-pathways_markers_E_2_15[pathways_markers_E_2_15$p.adjust<0.1,]
pathways_markers_E_13_15_sig<-pathways_markers_E_13_15[pathways_markers_E_13_15$p.adjust<0.1,]

pathways_markers_E_1_2_sig<-pathways_markers_E_1_2[pathways_markers_E_1_2$p.adjust<0.1,]
pathways_markers_E_6_2_sig<-pathways_markers_E_6_2[pathways_markers_E_6_2$p.adjust<0.1,]
pathways_markers_E_15_2_sig<-pathways_markers_E_15_2[pathways_markers_E_15_2$p.adjust<0.1,]
pathways_markers_E_8_2_sig<-pathways_markers_E_8_2[pathways_markers_E_8_2$p.adjust<0.1,]
pathways_markers_E_13_2_sig<-pathways_markers_E_13_2[pathways_markers_E_13_2$p.adjust<0.1,]

pathways_markers_E_1_13_sig<-pathways_markers_E_1_13[pathways_markers_E_1_13$p.adjust<0.1,]
pathways_markers_E_6_13_sig<-pathways_markers_E_6_13[pathways_markers_E_6_13$p.adjust<0.1,]
pathways_markers_E_15_13_sig<-pathways_markers_E_15_13[pathways_markers_E_15_13$p.adjust<0.1,]
pathways_markers_E_2_13_sig<-pathways_markers_E_2_13[pathways_markers_E_2_13$p.adjust<0.1,]
pathways_markers_E_8_13_sig<-pathways_markers_E_8_13[pathways_markers_E_8_13$p.adjust<0.1,]

# STEP 2: CREATE PATHWAY OBJECT
object <- ObjectCreator(Pathways = c(as.character(pathways_markers_BA_sig$ID),
                                     as.character(pathways_markers_CB_sig$ID),
                                     as.character(pathways_markers_DC_sig$ID),
                                     as.character(pathways_markers_ED_sig$ID),
                                     as.character(pathways_markers_E_1_8_sig$ID), 
                                     as.character(pathways_markers_E_6_8_sig$ID), 
                                     as.character(pathways_markers_E_15_8_sig$ID), 
                                     as.character(pathways_markers_E_2_8_sig$ID), 
                                     as.character(pathways_markers_E_13_8_sig$ID), 
                                     as.character(pathways_markers_E_8_1_sig$ID), 
                                     as.character(pathways_markers_E_6_1_sig$ID), 
                                     as.character(pathways_markers_E_15_1_sig$ID), 
                                     as.character(pathways_markers_E_2_1_sig$ID), 
                                     as.character(pathways_markers_E_13_1_sig$ID), 
                                     as.character(pathways_markers_E_1_6_sig$ID), 
                                     as.character(pathways_markers_E_8_6_sig$ID), 
                                     as.character(pathways_markers_E_15_6_sig$ID), 
                                     as.character(pathways_markers_E_2_6_sig$ID), 
                                     as.character(pathways_markers_E_13_6_sig$ID), 
                                     as.character(pathways_markers_E_1_15_sig$ID), 
                                     as.character(pathways_markers_E_6_15_sig$ID), 
                                     as.character(pathways_markers_E_8_15_sig$ID), 
                                     as.character(pathways_markers_E_2_15_sig$ID), 
                                     as.character(pathways_markers_E_13_15_sig$ID), 
                                     as.character(pathways_markers_E_1_2_sig$ID), 
                                     as.character(pathways_markers_E_6_2_sig$ID), 
                                     as.character(pathways_markers_E_15_2_sig$ID), 
                                     as.character(pathways_markers_E_8_2_sig$ID), 
                                     as.character(pathways_markers_E_13_2_sig$ID), 
                                     as.character(pathways_markers_E_1_13_sig$ID), 
                                     as.character(pathways_markers_E_6_13_sig$ID), 
                                     as.character(pathways_markers_E_15_13_sig$ID), 
                                     as.character(pathways_markers_E_2_13_sig$ID), 
                                     as.character(pathways_markers_E_8_13_sig$ID)), 
                        Molecules =c(as.character(pathways_markers_BA_sig$core_enrichment),
                                     as.character(pathways_markers_CB_sig$core_enrichment),
                                     as.character(pathways_markers_DC_sig$core_enrichment),
                                     as.character(pathways_markers_ED_sig$core_enrichment),
                                     as.character(pathways_markers_E_1_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_13_sig$core_enrichment)),
                        Groups=c(rep("BA",times=nrow(pathways_markers_BA_sig)),
                                 rep("CB",times=nrow(pathways_markers_CB_sig)),
                                 rep("DC",times=nrow(pathways_markers_DC_sig)),
                                 rep("ED",times=nrow(pathways_markers_ED_sig)),
                                 rep("E_1_8",times=nrow(pathways_markers_E_1_8_sig)),
                                 rep("E_6_8",times=nrow(pathways_markers_E_6_8_sig)),
                                 rep("E_15_8",times=nrow(pathways_markers_E_15_8_sig)),
                                 rep("E_2_8",times=nrow(pathways_markers_E_2_8_sig)),
                                 rep("E_13_8",times=nrow(pathways_markers_E_13_8_sig)),
                                 rep("E_8_1",times=nrow(pathways_markers_E_8_1_sig)),
                                 rep("E_6_1",times=nrow(pathways_markers_E_6_1_sig)),
                                 rep("E_15_1",times=nrow(pathways_markers_E_15_1_sig)),
                                 rep("E_2_1",times=nrow(pathways_markers_E_2_1_sig)),
                                 rep("E_13_1",times=nrow(pathways_markers_E_13_1_sig)),
                                 rep("E_1_6",times=nrow(pathways_markers_E_1_6_sig)),
                                 rep("E_8_6",times=nrow(pathways_markers_E_8_6_sig)),
                                 rep("E_15_6",times=nrow(pathways_markers_E_15_6_sig)),
                                 rep("E_2_6",times=nrow(pathways_markers_E_2_6_sig)),
                                 rep("E_13_6",times=nrow(pathways_markers_E_13_6_sig)),
                                 rep("E_1_15",times=nrow(pathways_markers_E_1_15_sig)),
                                 rep("E_6_15",times=nrow(pathways_markers_E_6_15_sig)),
                                 rep("E_8_15",times=nrow(pathways_markers_E_8_15_sig)),
                                 rep("E_2_15",times=nrow(pathways_markers_E_2_15_sig)),
                                 rep("E_13_15",times=nrow(pathways_markers_E_13_15_sig)),
                                 rep("E_1_2",times=nrow(pathways_markers_E_1_2_sig)),
                                 rep("E_6_2",times=nrow(pathways_markers_E_6_2_sig)),
                                 rep("E_15_2",times=nrow(pathways_markers_E_15_2_sig)),
                                 rep("E_8_2",times=nrow(pathways_markers_E_8_2_sig)),
                                 rep("E_13_2",times=nrow(pathways_markers_E_13_2_sig)),
                                 rep("E_1_13",times=nrow(pathways_markers_E_1_13_sig)),
                                 rep("E_6_13",times=nrow(pathways_markers_E_6_13_sig)),
                                 rep("E_15_13",times=nrow(pathways_markers_E_15_13_sig)),
                                 rep("E_2_13",times=nrow(pathways_markers_E_2_13_sig)),
                                 rep("E_8_13",times=nrow(pathways_markers_E_8_13_sig))),
                        structure = "SYMBOL", Type="", sep="/", Source="GSEA", organism="org.Hs.eg.db")
                        
# STEP3: COMBINE
combine<- CombineGeneSets(object, display="Expanded")

####
OptimalGeneSets(combine, method = "gap",cluster_method = "kmeans", max_cluster = 20, main="combine")

combine2 <- combine
combine2@metadata$Groups <- c("BA","CB","DC","ED",
                      "E_1_8","E_6_8","E_15_8","E_2_8","E_13_8",
                      "E_8_1","E_6_1","E_15_1","E_2_1","E_13_1",
                      "E_1_6","E_8_6","E_15_6","E_2_6","E_13_6",
                      "E_1_15","E_6_15","E_8_15","E_2_15","E_13_15",
                      "E_1_2","E_6_2","E_15_2","E_8_2","E_13_2",
                      "E_1_13","E_6_13","E_15_13","E_2_13","E_8_13")
combine2 <- ClusterGeneSets(combine2, clusters = 6, method = "kmeans",order = "cluster")

x <- combine2@plot$aka2
x <- x[,!colnames(x) == "Group"]
combine2@plot$aka2 <- x

x <- combine2@plot$aka3
x <- x[2:length(x)]
combine2@plot$aka3 <- x

pdf("genesetcluster_3_new.pdf", width = 12, height = 15)
PlotGeneSets(Object = combine2, main = "GO_BP_Neurons", RR.max = 20, annotation.mol=F)
dev.off()

####Try diferent onfigurations subsetting some groups
#object2 <- ObjectCreator(Pathways = c(as.character(pathways_markers_BA_sig$ID),
#                                     as.character(pathways_markers_CB_sig$ID),
#                                     as.character(pathways_markers_DC_sig$ID),
#                                     as.character(pathways_markers_ED_sig$ID)),
#                        Molecules =c(as.character(pathways_markers_BA_sig$core_enrichment),
#                                     as.character(pathways_markers_DC_sig$core_enrichment),
#                                     as.character(pathways_markers_ED_sig$core_enrichment)),
#                        Groups=c(rep("BA",times=nrow(pathways_markers_BA_sig)),
#                                 rep("DC",times=nrow(pathways_markers_DC_sig)),
#                                 rep("ED",times=nrow(pathways_markers_ED_sig))),
#                        structure = "SYMBOL", Type="", sep="/", Source="GSEA", organism="org.Hs.eg.db")
#combine3<- CombineGeneSets(object2, display="Expanded")
#OptimalGeneSets(combine3, method = "gap",cluster_method = "kmeans", max_cluster = 20, main="combine")
#combine3@metadata$Groups <- c("BA","DC","ED")
#combine3 <- ClusterGeneSets(combine3, clusters = 6, method = "kmeans",order = "cluster")

x <- combine3@plot$aka2
x <- x[,!colnames(x) == "Group"]
combine3@plot$aka2 <- x

x <- combine3@plot$aka3
x <- x[2:length(x)]
combine3@plot$aka3 <- x

png("genesetcluster_3_new_combine_3.png", width = 1200, height = 1500)
PlotGeneSets(Object = combine3, main = "GO_BP_Neurons", RR.max = 20, annotation.mol=F)
dev.off()

# STEP4: Export tables
clusters<-combine3@Data[[1]]
write.table(clusters, file="Clusters_pathways_descriptions_combine_3.csv", sep=";", row.names = FALSE)

annotation_pathways = bind_rows(list_of_gsea, .id = rownames(list_of_gsea))

#merge<-rbind(pathways_Monocytes_sig, pathways_DC_sig)
#merge<-merge[,c(1:2)]
all<-merge(clusters, annotation_pathways, by.x="Pathways", by.y="ID", all=FALSE)
write.table(all, file="Clusters_pathways_descriptions_combine_2.csv", sep=";", row.names = FALSE)


####DO limma (Ta-Tb) - (Ca-Cb)
library(limma)


#solo Tratamiento
Idents(combined_seurat) <- "protocol_2"
T_group_A = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 4 & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_B = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 0 & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_C = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 7 & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_D = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 3 & combined_seurat@meta.data$protocol_2 == "ba"]

T_group_E_8 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(8) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_1 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(1) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_6 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(6) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_15 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(15) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_13 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(13) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_2 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(2) & combined_seurat@meta.data$protocol_2 == "ba"]

T_group_E = c(T_group_E_8,T_group_E_1,T_group_E_6,T_group_E_15,T_group_E_13, T_group_E_2)

#solo Control
C_group_A = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 4 & combined_seurat@meta.data$protocol_2 == "control"]
C_group_B = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 0 & combined_seurat@meta.data$protocol_2 == "control"]
C_group_C = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 7 & combined_seurat@meta.data$protocol_2 == "control"]
C_group_D = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == 3 & combined_seurat@meta.data$protocol_2 == "control"]

C_group_E_8 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(8) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_1 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(1) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_6 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(6) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_15 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(15) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_13 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(13) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_2 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(2) & combined_seurat@meta.data$protocol_2 == "control"]

C_group_E = c(C_group_E_8,C_group_E_1,C_group_E_6,C_group_E_15,C_group_E_13,C_group_E_2)

combined_seurat@meta.data$contrasts_groups = "Rest"
combined_seurat@meta.data[T_group_A,]$contrasts_groups = "T_group_A"
combined_seurat@meta.data[T_group_B,]$contrasts_groups = "T_group_B"
combined_seurat@meta.data[T_group_C,]$contrasts_groups = "T_group_C"
combined_seurat@meta.data[T_group_D,]$contrasts_groups = "T_group_D"
combined_seurat@meta.data[T_group_E,]$contrasts_groups = "T_group_E"

combined_seurat@meta.data[C_group_A,]$contrasts_groups = "C_group_A"
combined_seurat@meta.data[C_group_B,]$contrasts_groups = "C_group_B"
combined_seurat@meta.data[C_group_C,]$contrasts_groups = "C_group_C"
combined_seurat@meta.data[C_group_D,]$contrasts_groups = "C_group_D"
combined_seurat@meta.data[C_group_E,]$contrasts_groups = "C_group_E"


#g <- factor(paste0(sex, drug, intervention))
#b <- factor(batch)
#design <- model.matrix(~0 + g + b) # plus other confounding covariates

g <- factor(combined_seurat@meta.data$contrasts_groups)
b <- factor(combined_seurat@meta.data$vector_days)
design <- model.matrix(~0 + g + b)
design <- design[,setdiff(colnames(design), "b4")] # get to full rank

exp_matrix = combined_seurat@assays$RNA@data
y <- voom(exp_matrix, design, plot = T)
#corfit <- duplicateCorrelation(exp_matrix, design, block = combined_seurat@meta.data$vector_days)

fit <- lmFit(exp_matrix, design)#, block = combined_seurat@meta.data$vector_days, correlation = corfit$consensus)

contrast.matrix <- makeContrasts((gT_group_E - gT_group_D) - (gC_group_E - gC_group_D), levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)

summary(decideTests(fit2))
diff_expression_results <- topTreat(fit2, coef=1, n=Inf)

# STEP2 - Diff and significative changes of genes in contrast for Tx-Ty treatment strict derived changes
potential_genes_BA = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_CB = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_DC = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_ED = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
write.csv(potential_genes_ED, file = "Potential_markers_step_2_ED.csv")

diff_expression_results_BA = diff_expression_results
diff_expression_results_CB = diff_expression_results
diff_expression_results_DC = diff_expression_results
diff_expression_results_ED = diff_expression_results

#PREPARING gsea input
list_of_contrasts_step2 = list(diff_expression_results_BA,diff_expression_results_CB,diff_expression_results_DC,diff_expression_results_ED)
names(list_of_contrasts_step2) = c("diff_expression_results_BA","diff_expression_results_CB","diff_expression_results_DC","diff_expression_results_ED")

OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/GO_step2/"

#Change names to entrez code
name_changer = function(especie_database, especie_name, x){
  ensembl <- useMart(especie_database, dataset = especie_name)
  annotation = getBM(attributes=c("hgnc_symbol","entrezgene_id"), values=x, mart=ensembl)
  return(annotation)
}

library("enrichplot")
library("clusterProfiler")
library("org.Hs.eg.db")
list_of_gsea = list()
for (n in 1:length(list_of_contrasts_step2)){
  print(n)
  #Create the geneList
  geneList <- list_of_contrasts_step2[[n]]$logFC
  names(geneList) <- as.character(rownames(list_of_contrasts_step2[[n]]))
  geneList = sort(geneList, decreasing = TRUE)
   
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]
  
  #GSEA
  ego <- gseGO(geneList    = geneList,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 100,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = FALSE)
  ego_table <- as.data.frame(ego) 
  write.csv(ego_table,paste0(OutPath,"gsea_","BP","_contrast_",names(list_of_contrasts_step2[n]),".csv"))
  list_of_gsea = append(list_of_gsea, list(ego_table))
  #plots
  if (nrow(ego_table) > 1){
    pdf(file=paste0(OutPath,"gsea_dotplot_","BP","_contrast_",names(list_of_contrasts_step2[n]),".pdf"), width=8, height=12)  
    plot(dotplot(ego, showCategory=30))
    dev.off()
    pdf(file=paste0(OutPath,"gsea_gseaPlot_","BP","_contrast_",names(list_of_contrasts_step2[n]),".pdf"), width=12, height=12)  
    plot(gseaplot2(ego, geneSetID = 1:3))
    dev.off()    
  }
}

#Preparing genesetcluster for step2
names(list_of_gsea) = c("markers_BA","markers_CB","markers_DC","markers_ED")

# STEP 1: with all samples included
pathways_markers_BA = list_of_gsea[["markers_BA"]]
pathways_markers_CB = list_of_gsea[["markers_CB"]]
pathways_markers_DC = list_of_gsea[["markers_DC"]]
pathways_markers_ED = list_of_gsea[["markers_ED"]]

#get significant pathways

pathways_markers_BA_sig<-pathways_markers_BA[pathways_markers_BA$p.adjust<0.1,]
pathways_markers_CB_sig<-pathways_markers_CB[pathways_markers_CB$p.adjust<0.1,]
pathways_markers_DC_sig<-pathways_markers_DC[pathways_markers_DC$p.adjust<0.1,]
pathways_markers_ED_sig<-pathways_markers_ED[pathways_markers_ED$p.adjust<0.1,]

# STEP 2: CREATE PATHWAY OBJECT
object <- ObjectCreator(Pathways = c(as.character(pathways_markers_BA_sig$ID),
                                     as.character(pathways_markers_CB_sig$ID),
                                     as.character(pathways_markers_DC_sig$ID),
                                     as.character(pathways_markers_ED_sig$ID)),
                        Molecules =c(as.character(pathways_markers_BA_sig$core_enrichment),
                                     as.character(pathways_markers_CB_sig$core_enrichment),
                                     as.character(pathways_markers_DC_sig$core_enrichment),
                                     as.character(pathways_markers_ED_sig$core_enrichment)),
                        Groups=c(rep("BA",times=nrow(pathways_markers_BA_sig)),
                                 rep("CB",times=nrow(pathways_markers_CB_sig)),
                                 rep("DC",times=nrow(pathways_markers_DC_sig)),
                                 rep("ED",times=nrow(pathways_markers_ED_sig))),
                        structure = "SYMBOL", Type="", sep="/", Source="GSEA", organism="org.Hs.eg.db")
                        
# STEP3: COMBINE
combine<- CombineGeneSets(object, display="Expanded")

####
OptimalGeneSets(combine, method = "gap",cluster_method = "kmeans", max_cluster = 20, main="combine")

combine2 <- combine
combine2@metadata$Groups <- c("BA","DC","ED")
combine2 <- ClusterGeneSets(combine2, clusters = 6, method = "kmeans",order = "cluster")

x <- combine2@plot$aka2
x <- x[,!colnames(x) == "Group"]
combine2@plot$aka2 <- x

x <- combine2@plot$aka3
x <- x[2:length(x)]
combine2@plot$aka3 <- x

pdf("genesetcluster_3_new.pdf", width = 12, height = 15)
PlotGeneSets(Object = combine2, main = "GO_BP_Neurons", RR.max = 20, annotation.mol=F)
dev.off()

# STEP4: Export tables
clusters<-combine2@Data[[1]]
write.table(clusters, file="Clusters_pathways_descriptions.csv", sep=";", row.names = FALSE)

annotation_pathways = bind_rows(list_of_gsea, .id = rownames(list_of_gsea))

all<-merge(clusters, annotation_pathways, by.x="Pathways", by.y="ID", all=FALSE)
write.table(all, file="Clusters_pathways_descriptions.csv", sep=";", row.names = FALSE)


######################################################################################################################################################
# STEP 3 - create matrix with pertainti of genes on each category

names_sig_C_markers_BA = rownames(potential_markers_BA)
names_sig_C_markers_CB = rownames(potential_markers_CB)
names_sig_C_markers_DC = rownames(potential_markers_DC)
names_sig_C_markers_ED = rownames(potential_markers_ED)
write.table(potential_markers_ED, file="Potential_markers_step_1_ED.csv", sep=";", row.names = TRUE)

names_sig_T_potential_BA = rownames(potential_genes_BA)
names_sig_T_potential_CB = rownames(potential_genes_CB)
names_sig_T_potential_DC = rownames(potential_genes_DC)
names_sig_T_potential_ED = rownames(potential_genes_ED)


out <- table(stack(mget(ls(pattern = "^names_sig_"))))
out <- out[, c("names_sig_C_markers_BA", "names_sig_C_markers_CB", "names_sig_C_markers_DC", "names_sig_C_markers_ED", "names_sig_T_potential_BA", "names_sig_T_potential_CB", "names_sig_T_potential_DC", "names_sig_T_potential_ED")]

out_table = data.frame(matrix(data = out, ncol = ncol(out), nrow = nrow(out)))
rownames(out_table) = rownames(out)
colnames(out_table) = colnames(out)
out_table$Events = rowSums(out_table)
write.table(out_table, file="genes_of_interest_per_contrast.csv", sep=";", row.names = TRUE)

out_summary = t(out) %*% out
write.table(out_summary, file="genes_of_interest_per_contrast_summary.csv", sep=";", row.names = TRUE)

#GENES IN BOTH CONTROL + TREATMENT
lista_gene_names=Reduce(intersect, list(names_sig_T_potential_ED,names_sig_C_markers_CB))

potential_markers_CB[lista_gene_names,]
potential_genes_ED[lista_gene_names,]

#GENES ONLY IN CONTROL
lista_gene_names = names_sig_C_markers_ED[!(names_sig_C_markers_CB %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
potential_markers_ED[lista_gene_names,]
#GENES ONLY IN TREATMENT
lista_gene_names = names_sig_T_potential_ED[!(names_sig_T_potential_ED %in% c(names_sig_C_markers_BA,names_sig_C_markers_CB,names_sig_C_markers_DC,names_sig_C_markers_ED))]
potential_genes_ED[lista_gene_names,]


#For heatmap rownames colors
both_conditions_all=Reduce(intersect, list(names_sig_T_potential_BA,names_sig_C_markers_BA))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_BA,names_sig_C_markers_CB)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_BA,names_sig_C_markers_DC)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_BA,names_sig_C_markers_ED)))

both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_CB,names_sig_C_markers_BA)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_CB,names_sig_C_markers_CB)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_CB,names_sig_C_markers_DC)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_CB,names_sig_C_markers_ED)))

both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_DC,names_sig_C_markers_BA)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_DC,names_sig_C_markers_CB)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_DC,names_sig_C_markers_DC)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_DC,names_sig_C_markers_ED)))

both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_ED,names_sig_C_markers_BA)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_ED,names_sig_C_markers_CB)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_ED,names_sig_C_markers_DC)))
both_conditions_all=append(both_conditions_all, Reduce(intersect, list(names_sig_T_potential_ED,names_sig_C_markers_ED)))
both_conditions_all=unique(both_conditions_all)

only_T_conditions_all=c("NOL7","CLUAP1","SEC61G","SRP14","SLIRP","SEC61G","SNRPF","TMEM59","ROMO1","UBE2E3","STRAP","ENAH")

only_C_conditions_all_BA = names_sig_C_markers_BA[!(names_sig_C_markers_BA %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_conditions_all_CB = names_sig_C_markers_CB[!(names_sig_C_markers_CB %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_conditions_all_DC = names_sig_C_markers_DC[!(names_sig_C_markers_DC %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_conditions_all_ED = names_sig_C_markers_ED[!(names_sig_C_markers_ED %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_condition_all = unique(c(only_C_conditions_all_BA,only_C_conditions_all_CB,only_C_conditions_all_DC,only_C_conditions_all_ED))

######################################################################################################################################################
#HEATMAP PLOT OF BINARIZED AND EXPRESION CASCADES OF ACTIVATION FOR ALL GENES SIGNIFICANT CHANGING
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize) # for the colorRamp2() function

#Table of binarized cascade activation + event(sum of cascade activations)
out_table_plot = out_table
out_table_plot$Events = NULL

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_BA","C_CB","C_DC","C_ED","T_BA","T_CB","T_DC","T_ED")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("BA","CB","DC","ED","BA","CB","DC","ED")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("BA" = "#d01c8b",
                                           "CB" = "#f1b6da",
                                           "DC" = "#b8e186",
                                           "ED" = "#4dac26")
                          )
                          )
ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot
colnames(my_data) = rownames(meta_data_plot)
#Binarized version of table and plot
# Heatmap
col_fun = colorRamp2(c(0, 1), c("#f7f7f7", "#998ec3"))

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))


pdf(file = "Binary_cascade_of_activation_transition_lvl_1.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

######################################################################################################################################################
#REPEAT HEATMAP WITH EXPRESSION VALUES
#Table of binarized cascade activation + event(sum of cascade activations)
out_table_plot = out_table
out_table_plot$Events = NULL

#Fom binary to logFC matrix
out_table_plot[,1] = replace(out_table_plot[,1], out_table_plot[,1]==1, potential_markers_BA[rownames(out_table_plot)[out_table_plot[,1]==1],]$avg_log2FC)
out_table_plot[,2] = replace(out_table_plot[,2], out_table_plot[,2]==1, potential_markers_CB[rownames(out_table_plot)[out_table_plot[,2]==1],]$avg_log2FC)
out_table_plot[,3] = replace(out_table_plot[,3], out_table_plot[,3]==1, potential_markers_DC[rownames(out_table_plot)[out_table_plot[,3]==1],]$avg_log2FC)
out_table_plot[,4] = replace(out_table_plot[,4], out_table_plot[,4]==1, potential_markers_ED[rownames(out_table_plot)[out_table_plot[,4]==1],]$avg_log2FC)

out_table_plot[,5] = replace(out_table_plot[,5], out_table_plot[,5]==1, potential_genes_BA[rownames(out_table_plot)[out_table_plot[,5]==1],]$logFC)
out_table_plot[,6] = replace(out_table_plot[,6], out_table_plot[,6]==1, potential_genes_CB[rownames(out_table_plot)[out_table_plot[,6]==1],]$logFC)
out_table_plot[,7] = replace(out_table_plot[,7], out_table_plot[,7]==1, potential_genes_DC[rownames(out_table_plot)[out_table_plot[,7]==1],]$logFC)
out_table_plot[,8] = replace(out_table_plot[,8], out_table_plot[,8]==1, potential_genes_ED[rownames(out_table_plot)[out_table_plot[,8]==1],]$logFC)

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_BA","C_CB","C_DC","C_ED","T_BA","T_CB","T_DC","T_ED")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("BA","CB","DC","ED","BA","CB","DC","ED")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("BA" = "#d01c8b",
                                           "CB" = "#f1b6da",
                                           "DC" = "#b8e186",
                                           "ED" = "#4dac26")
                          )
                          
ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_1.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

col_fun = colorRamp2(c(my_data_min/100, 0, my_data_max/100), c("blue", "white", "red"))
pdf(file = "Expression_2 colors_cascade_of_activation_transition_lvl_1.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment", at = c(-1, 0, 1), labels = c("Sig.Down.regulated", "No.significant", "Sig.Up.regulated")),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

#GET 285 TF from JASPAR2020 285
ocr_tf_matrix=read.table("ocr_tf_matrix.txt", sep = "")
names(ocr_tf_matrix) 

library(JASPAR2022)
getJasparMotifs = function(species = "Homo sapiens", collection = "CORE"){
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- collection
    out <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
    if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
        names(out) <- paste(names(out), TFBSTools::name(out),
            sep = "_")
    return(out)
}
pfm <- getJasparMotifs()
sub(".*_", "", names(pfm)

#191 in common from 285 vs 692
sum(names(ocr_tf_matrix) %in% sub(".*_", "", names(pfm)))

#How many TF do we have in out 1728? 
sum(rownames(my_data) %in% names(ocr_tf_matrix)) #32
sum(rownames(my_data) %in% sub(".*_", "", names(pfm))) #67 USE JASPAR 2022
tf_in_our_contrasts= rownames(my_data)[rownames(my_data) %in% sub(".*_", "", names(pfm))]

write.table(out_table[tf_in_our_contrasts,], file="TF_binary_of_interest_per_contrast.csv", sep=";", row.names = TRUE)
write.table(out_table_plot[tf_in_our_contrasts,], file="TF_of_interest_per_contrast.csv", sep=";", row.names = TRUE)

##################################################################################################################################################
##################################################################################################################################################
#HEATMAPS SUMMARY 1-c+t+only_T AND 2-only_C
genelist_both_and_T_only = c(both_conditions_all,only_T_conditions_all)

out_table_plot_genelist_both_and_T_only= out_table_plot[genelist_both_and_T_only,]
# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_BA","C_CB","C_DC","C_ED","T_BA","T_CB","T_DC","T_ED")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("BA","CB","DC","ED","BA","CB","DC","ED")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("BA" = "#d01c8b",
                                           "CB" = "#f1b6da",
                                           "DC" = "#b8e186",
                                           "ED" = "#4dac26")
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_genelist_both_and_T_only
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('black', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_1_genelist_both_and_T_only.pdf",
    width = 10,
    height = 10)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

##########
#2-only_C

#Order base on logFC
temp = potential_markers_BA[order(potential_markers_BA$avg_log2FC),]
temp_BA = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_CB[order(potential_markers_CB$avg_log2FC),]
temp_CB = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_DC[order(potential_markers_DC$avg_log2FC),]
temp_DC = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_ED[order(potential_markers_ED$avg_log2FC),]
temp_ED = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp_all = unique(c(temp_BA,temp_CB,temp_DC,temp_ED)) #116

only_C_condition_all_subset = only_C_condition_all[only_C_condition_all %in% temp_all]
out_table_plot_genelist_C_only= out_table_plot[only_C_condition_all_subset,]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_BA","C_CB","C_DC","C_ED","T_BA","T_CB","T_DC","T_ED")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("BA","CB","DC","ED","BA","CB","DC","ED")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("BA" = "#d01c8b",
                                           "CB" = "#f1b6da",
                                           "DC" = "#b8e186",
                                           "ED" = "#4dac26")
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_genelist_C_only
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_C_condition_all_subset
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('black', nrow(my_data))
fontcolors[row_idx] <- 'd8b365'
fontcolors[row_idx_t] <- '#d8b365'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_1_genelist_C_only.pdf",
    width = 10,
    height = 10)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

##########
#Caracterization
genelist_caracterization = c("S100B","ANXA2","TUBB2A","SOX4","PROX1","DCX","NEUROG2","ELAVL3","VCAM1","STAT3","FYN","HES6","FOS","HES5","ASCL1","VIM","SOX3","PAX6","BAX","SOX21","LGALSL","CALB1","NES","SOX11","NCAM1","SOX9","SOX21","REST","MALAT1")
out_table_plot_genelist_caracterization= out_table_plot[genelist_caracterization,]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_BA","C_CB","C_DC","C_ED","T_BA","T_CB","T_DC","T_ED")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("BA","CB","DC","ED","BA","CB","DC","ED")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("BA" = "#d01c8b",
                                           "CB" = "#f1b6da",
                                           "DC" = "#b8e186",
                                           "ED" = "#4dac26")
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_genelist_caracterization
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_1_genelist_caracterization.pdf",
    width = 10,
    height = 5)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()








######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#NEURONS-LVL2 CLUSTER IDENTITY SETTING WITH 1 VS ALL CONTRASTS 
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
Idents(only_control_cells_subset) = "integrated_snn_res.0.6"
only_control_cells_subset_neurons = subset(x = only_control_cells_subset, idents=c("1","2","6","8","13","15"))

#########################################
#ONE vs ALL significant identity markers
control_all_markers_sub_neurons = FindAllMarkers(object = only_control_cells_subset_neurons, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.125)
write.table(control_all_markers_sub_neurons, file="Identity_markers_sub_neurons_only_Control_cells.csv", sep=";", row.names = TRUE)
#ENOUGH FOR KEEP GOING

setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/level_2_neurons/")


# STEP 0: First generate list_of_gsea with previous code withoit regular 4 neurogenesis groups

names(list_of_gsea) = c("markers_E_1_8","markers_E_6_8","markers_E_15_8","markers_E_2_8","markers_E_13_8",
                      "markers_E_8_1","markers_E_6_1","markers_E_15_1","markers_E_2_1","markers_E_13_1",
                      "markers_E_1_6","markers_E_8_6","markers_E_15_6","markers_E_2_6","markers_E_13_6",
                      "markers_E_1_15","markers_E_6_15","markers_E_8_15","markers_E_2_15","markers_E_13_15",
                      "markers_E_1_2","markers_E_6_2","markers_E_15_2","markers_E_8_2","markers_E_13_2",
                      "markers_E_1_13","markers_E_6_13","markers_E_15_13","markers_E_2_13","markers_E_8_13")

# STEP 1: with all samples included
pathways_markers_E_1_8 = list_of_gsea[["markers_E_1_8"]]
pathways_markers_E_6_8 = list_of_gsea[["markers_E_6_8"]]
pathways_markers_E_15_8 = list_of_gsea[["markers_E_15_8"]]
pathways_markers_E_2_8 = list_of_gsea[["markers_E_2_8"]]
pathways_markers_E_13_8 = list_of_gsea[["markers_E_13_8"]]

pathways_markers_E_8_1 = list_of_gsea[["markers_E_8_1"]]
pathways_markers_E_6_1 = list_of_gsea[["markers_E_6_1"]]
pathways_markers_E_15_1 = list_of_gsea[["markers_E_15_1"]]
pathways_markers_E_2_1 = list_of_gsea[["markers_E_2_1"]]
pathways_markers_E_13_1 = list_of_gsea[["markers_E_13_1"]]

pathways_markers_E_1_6 = list_of_gsea[["markers_E_1_6"]]
pathways_markers_E_8_6 = list_of_gsea[["markers_E_8_6"]]
pathways_markers_E_15_6 = list_of_gsea[["markers_E_15_6"]]
pathways_markers_E_2_6 = list_of_gsea[["markers_E_2_6"]]
pathways_markers_E_13_6 = list_of_gsea[["markers_E_13_6"]]

pathways_markers_E_1_15 = list_of_gsea[["markers_E_1_15"]]
pathways_markers_E_6_15 = list_of_gsea[["markers_E_6_15"]]
pathways_markers_E_8_15 = list_of_gsea[["markers_E_8_15"]]
pathways_markers_E_2_15 = list_of_gsea[["markers_E_2_15"]]
pathways_markers_E_13_15 = list_of_gsea[["markers_E_13_15"]]

pathways_markers_E_1_2 = list_of_gsea[["markers_E_1_2"]]
pathways_markers_E_6_2 = list_of_gsea[["markers_E_6_2"]]
pathways_markers_E_15_2 = list_of_gsea[["markers_E_15_2"]]
pathways_markers_E_8_2 = list_of_gsea[["markers_E_8_2"]]
pathways_markers_E_13_2 = list_of_gsea[["markers_E_13_2"]]

pathways_markers_E_1_13 = list_of_gsea[["markers_E_1_13"]]
pathways_markers_E_6_13 = list_of_gsea[["markers_E_6_13"]]
pathways_markers_E_15_13 = list_of_gsea[["markers_E_15_13"]]
pathways_markers_E_2_13 = list_of_gsea[["markers_E_2_13"]]
pathways_markers_E_8_13 = list_of_gsea[["markers_E_8_13"]]

#get significant pathways
pathways_markers_E_1_8_sig<-pathways_markers_E_1_8[pathways_markers_E_1_8$p.adjust<0.1,]
pathways_markers_E_6_8_sig<-pathways_markers_E_6_8[pathways_markers_E_6_8$p.adjust<0.1,]
pathways_markers_E_15_8_sig<-pathways_markers_E_15_8[pathways_markers_E_15_8$p.adjust<0.1,]
pathways_markers_E_2_8_sig<-pathways_markers_E_2_8[pathways_markers_E_2_8$p.adjust<0.1,]
pathways_markers_E_13_8_sig<-pathways_markers_E_13_8[pathways_markers_E_13_8$p.adjust<0.1,]

pathways_markers_E_8_1_sig<-pathways_markers_E_8_1[pathways_markers_E_8_1$p.adjust<0.1,]
pathways_markers_E_6_1_sig<-pathways_markers_E_6_1[pathways_markers_E_6_1$p.adjust<0.1,]
pathways_markers_E_15_1_sig<-pathways_markers_E_15_1[pathways_markers_E_15_1$p.adjust<0.1,]
pathways_markers_E_2_1_sig<-pathways_markers_E_2_1[pathways_markers_E_2_1$p.adjust<0.1,]
pathways_markers_E_13_1_sig<-pathways_markers_E_13_1[pathways_markers_E_13_1$p.adjust<0.1,]

pathways_markers_E_1_6_sig<-pathways_markers_E_1_6[pathways_markers_E_1_6$p.adjust<0.1,]
pathways_markers_E_8_6_sig<-pathways_markers_E_8_6[pathways_markers_E_8_6$p.adjust<0.1,]
pathways_markers_E_15_6_sig<-pathways_markers_E_15_6[pathways_markers_E_15_6$p.adjust<0.1,]
pathways_markers_E_2_6_sig<-pathways_markers_E_2_6[pathways_markers_E_2_6$p.adjust<0.1,]
pathways_markers_E_13_6_sig<-pathways_markers_E_13_6[pathways_markers_E_13_6$p.adjust<0.1,]

pathways_markers_E_1_15_sig<-pathways_markers_E_1_15[pathways_markers_E_1_15$p.adjust<0.1,]
pathways_markers_E_6_15_sig<-pathways_markers_E_6_15[pathways_markers_E_6_15$p.adjust<0.1,]
pathways_markers_E_8_15_sig<-pathways_markers_E_8_15[pathways_markers_E_8_15$p.adjust<0.1,]
pathways_markers_E_2_15_sig<-pathways_markers_E_2_15[pathways_markers_E_2_15$p.adjust<0.1,]
pathways_markers_E_13_15_sig<-pathways_markers_E_13_15[pathways_markers_E_13_15$p.adjust<0.1,]

pathways_markers_E_1_2_sig<-pathways_markers_E_1_2[pathways_markers_E_1_2$p.adjust<0.1,]
pathways_markers_E_6_2_sig<-pathways_markers_E_6_2[pathways_markers_E_6_2$p.adjust<0.1,]
pathways_markers_E_15_2_sig<-pathways_markers_E_15_2[pathways_markers_E_15_2$p.adjust<0.1,]
pathways_markers_E_8_2_sig<-pathways_markers_E_8_2[pathways_markers_E_8_2$p.adjust<0.1,]
pathways_markers_E_13_2_sig<-pathways_markers_E_13_2[pathways_markers_E_13_2$p.adjust<0.1,]

pathways_markers_E_1_13_sig<-pathways_markers_E_1_13[pathways_markers_E_1_13$p.adjust<0.1,]
pathways_markers_E_6_13_sig<-pathways_markers_E_6_13[pathways_markers_E_6_13$p.adjust<0.1,]
pathways_markers_E_15_13_sig<-pathways_markers_E_15_13[pathways_markers_E_15_13$p.adjust<0.1,]
pathways_markers_E_2_13_sig<-pathways_markers_E_2_13[pathways_markers_E_2_13$p.adjust<0.1,]
pathways_markers_E_8_13_sig<-pathways_markers_E_8_13[pathways_markers_E_8_13$p.adjust<0.1,]

# STEP 2: CREATE PATHWAY OBJECT
object <- ObjectCreator(Pathways = c(as.character(pathways_markers_E_1_8_sig$ID), 
                                     as.character(pathways_markers_E_6_8_sig$ID), 
                                     as.character(pathways_markers_E_15_8_sig$ID), 
                                     as.character(pathways_markers_E_2_8_sig$ID), 
                                     as.character(pathways_markers_E_13_8_sig$ID), 
                                     as.character(pathways_markers_E_8_1_sig$ID), 
                                     as.character(pathways_markers_E_6_1_sig$ID), 
                                     as.character(pathways_markers_E_15_1_sig$ID), 
                                     as.character(pathways_markers_E_2_1_sig$ID), 
                                     as.character(pathways_markers_E_13_1_sig$ID), 
                                     as.character(pathways_markers_E_1_6_sig$ID), 
                                     as.character(pathways_markers_E_8_6_sig$ID), 
                                     as.character(pathways_markers_E_15_6_sig$ID), 
                                     as.character(pathways_markers_E_2_6_sig$ID), 
                                     as.character(pathways_markers_E_13_6_sig$ID), 
                                     as.character(pathways_markers_E_1_15_sig$ID), 
                                     as.character(pathways_markers_E_6_15_sig$ID), 
                                     as.character(pathways_markers_E_8_15_sig$ID), 
                                     as.character(pathways_markers_E_2_15_sig$ID), 
                                     as.character(pathways_markers_E_13_15_sig$ID), 
                                     as.character(pathways_markers_E_1_2_sig$ID), 
                                     as.character(pathways_markers_E_6_2_sig$ID), 
                                     as.character(pathways_markers_E_15_2_sig$ID), 
                                     as.character(pathways_markers_E_8_2_sig$ID), 
                                     as.character(pathways_markers_E_13_2_sig$ID), 
                                     as.character(pathways_markers_E_1_13_sig$ID), 
                                     as.character(pathways_markers_E_6_13_sig$ID), 
                                     as.character(pathways_markers_E_15_13_sig$ID), 
                                     as.character(pathways_markers_E_2_13_sig$ID), 
                                     as.character(pathways_markers_E_8_13_sig$ID)), 
                        Molecules =c(as.character(pathways_markers_E_1_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_8_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_1_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_6_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_15_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_13_2_sig$core_enrichment), 
                                     as.character(pathways_markers_E_1_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_6_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_15_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_2_13_sig$core_enrichment), 
                                     as.character(pathways_markers_E_8_13_sig$core_enrichment)),
                        Groups=c(rep("E_1_8",times=nrow(pathways_markers_E_1_8_sig)),
                                 rep("E_6_8",times=nrow(pathways_markers_E_6_8_sig)),
                                 rep("E_15_8",times=nrow(pathways_markers_E_15_8_sig)),
                                 rep("E_2_8",times=nrow(pathways_markers_E_2_8_sig)),
                                 rep("E_13_8",times=nrow(pathways_markers_E_13_8_sig)),
                                 rep("E_8_1",times=nrow(pathways_markers_E_8_1_sig)),
                                 rep("E_6_1",times=nrow(pathways_markers_E_6_1_sig)),
                                 rep("E_15_1",times=nrow(pathways_markers_E_15_1_sig)),
                                 rep("E_2_1",times=nrow(pathways_markers_E_2_1_sig)),
                                 rep("E_13_1",times=nrow(pathways_markers_E_13_1_sig)),
                                 rep("E_1_6",times=nrow(pathways_markers_E_1_6_sig)),
                                 rep("E_8_6",times=nrow(pathways_markers_E_8_6_sig)),
                                 rep("E_15_6",times=nrow(pathways_markers_E_15_6_sig)),
                                 rep("E_2_6",times=nrow(pathways_markers_E_2_6_sig)),
                                 rep("E_13_6",times=nrow(pathways_markers_E_13_6_sig)),
                                 rep("E_1_15",times=nrow(pathways_markers_E_1_15_sig)),
                                 rep("E_6_15",times=nrow(pathways_markers_E_6_15_sig)),
                                 rep("E_8_15",times=nrow(pathways_markers_E_8_15_sig)),
                                 rep("E_2_15",times=nrow(pathways_markers_E_2_15_sig)),
                                 rep("E_13_15",times=nrow(pathways_markers_E_13_15_sig)),
                                 rep("E_1_2",times=nrow(pathways_markers_E_1_2_sig)),
                                 rep("E_6_2",times=nrow(pathways_markers_E_6_2_sig)),
                                 rep("E_15_2",times=nrow(pathways_markers_E_15_2_sig)),
                                 rep("E_8_2",times=nrow(pathways_markers_E_8_2_sig)),
                                 rep("E_13_2",times=nrow(pathways_markers_E_13_2_sig)),
                                 rep("E_1_13",times=nrow(pathways_markers_E_1_13_sig)),
                                 rep("E_6_13",times=nrow(pathways_markers_E_6_13_sig)),
                                 rep("E_15_13",times=nrow(pathways_markers_E_15_13_sig)),
                                 rep("E_2_13",times=nrow(pathways_markers_E_2_13_sig)),
                                 rep("E_8_13",times=nrow(pathways_markers_E_8_13_sig))),
                        structure = "SYMBOL", Type="", sep="/", Source="GSEA", organism="org.Hs.eg.db")
                        
# STEP3: COMBINE
combine<- CombineGeneSets(object, display="Expanded")

####
OptimalGeneSets(combine, method = "gap",cluster_method = "kmeans", max_cluster = 20, main="combine")

combine2 <- combine
combine2@metadata$Groups <- c("E_1_8","E_6_8","E_15_8","E_2_8","E_13_8",
                      "E_8_1","E_6_1","E_15_1","E_2_1","E_13_1",
                      "E_1_6","E_8_6","E_15_6","E_2_6","E_13_6",
                      "E_1_15","E_6_15","E_8_15","E_2_15","E_13_15",
                      "E_1_2","E_6_2","E_15_2","E_8_2","E_13_2",
                      "E_1_13","E_6_13","E_15_13","E_2_13","E_8_13")
combine2 <- ClusterGeneSets(combine2, clusters = 8, method = "kmeans",order = "cluster")

x <- combine2@plot$aka2
x <- x[,!colnames(x) == "Group"]
combine2@plot$aka2 <- x

x <- combine2@plot$aka3
x <- x[2:length(x)]
combine2@plot$aka3 <- x

png("genesetcluster_lvl2_neurons_8c.png", width = 1200, height = 1500)
PlotGeneSets(Object = combine2, main = "GO_BP_Neurons", RR.max = 20, annotation.mol=F)
dev.off()

# STEP4: Export tables
clusters<-combine2@Data[[1]]
write.table(clusters, file="Clusters_pathways_descriptions_combine_lvl2_neurons_8c_.csv", sep=";", row.names = FALSE)

annotation_pathways = bind_rows(list_of_gsea, .id = rownames(list_of_gsea))

#merge<-rbind(pathways_Monocytes_sig, pathways_DC_sig)
#merge<-merge[,c(1:2)]
all<-merge(clusters, annotation_pathways, by.x="Pathways", by.y="ID", all=FALSE)
write.table(all, file="Clusters_pathways_descriptions_combine_lvl2_neurons_8c_.csv", sep=";", row.names = FALSE)

#cluster 6 y 2 que son????
write.table(markers_E_6_2, file="Clusters_pathways_descriptions_combine_lvl2_neurons_8c_markers_E_6_2.csv", sep=";", row.names = FALSE)


####DO limma (Ta-Tb) - (Ca-Cb)
library(limma)


#solo Tratamiento
Idents(combined_seurat) <- "protocol_2"
T_group_E_8 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(8) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_1 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(1) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_6 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(6) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_15 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(15) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_13 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(13) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_2 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(2) & combined_seurat@meta.data$protocol_2 == "ba"]

#T_group_E = c(T_group_E_8,T_group_E_1,T_group_E_6,T_group_E_15,T_group_E_13, T_group_E_2)

#solo Control
C_group_E_8 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(8) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_1 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(1) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_6 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(6) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_15 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(15) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_13 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(13) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_2 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(2) & combined_seurat@meta.data$protocol_2 == "control"]

#C_group_E = c(C_group_E_8,C_group_E_1,C_group_E_6,C_group_E_15,C_group_E_13,C_group_E_2)

combined_seurat@meta.data$contrasts_groups = "Rest"
combined_seurat@meta.data[T_group_E_8,]$contrasts_groups = "T_group_E_8"
combined_seurat@meta.data[T_group_E_1,]$contrasts_groups = "T_group_E_1"
combined_seurat@meta.data[T_group_E_6,]$contrasts_groups = "T_group_E_6"
combined_seurat@meta.data[T_group_E_15,]$contrasts_groups = "T_group_E_15"
combined_seurat@meta.data[T_group_E_13,]$contrasts_groups = "T_group_E_13"
combined_seurat@meta.data[T_group_E_2,]$contrasts_groups = "T_group_E_2"

combined_seurat@meta.data[C_group_E_8,]$contrasts_groups = "C_group_E_8"
combined_seurat@meta.data[C_group_E_1,]$contrasts_groups = "C_group_E_1"
combined_seurat@meta.data[C_group_E_6,]$contrasts_groups = "C_group_E_6"
combined_seurat@meta.data[C_group_E_15,]$contrasts_groups = "C_group_E_15"
combined_seurat@meta.data[C_group_E_13,]$contrasts_groups = "C_group_E_13"
combined_seurat@meta.data[C_group_E_2,]$contrasts_groups = "C_group_E_2"




#g <- factor(paste0(sex, drug, intervention))
#b <- factor(batch)
#design <- model.matrix(~0 + g + b) # plus other confounding covariates

g <- factor(combined_seurat@meta.data$contrasts_groups)
b <- factor(combined_seurat@meta.data$vector_days)
design <- model.matrix(~0 + g + b)
design <- design[,setdiff(colnames(design), "b4")] # get to full rank

exp_matrix = combined_seurat@assays$RNA@data
y <- voom(exp_matrix, design, plot = T)
#corfit <- duplicateCorrelation(exp_matrix, design, block = combined_seurat@meta.data$vector_days)

fit <- lmFit(exp_matrix, design)#, block = combined_seurat@meta.data$vector_days, correlation = corfit$consensus)

contrast.matrix <- makeContrasts((gT_group_E_2 - gT_group_E_8) - (gC_group_E_2 - gC_group_E_8), levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)

summary(decideTests(fit2))
diff_expression_results <- topTreat(fit2, coef=1, n=Inf)

# STEP2 - Diff and significative changes of genes in contrast for Tx-Ty treatment strict derived changes
potential_genes_E_1_13 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_E_15_1 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_E_6_15 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_E_2_15 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_E_8_1 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_E_6_8 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_E_2_8 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
write.csv(potential_genes_E_2_8, file = "Potential_markers_step_2_E_2_8.csv")

diff_expression_results_E_1_13 = diff_expression_results
diff_expression_results_E_15_1 = diff_expression_results
diff_expression_results_E_6_15 = diff_expression_results
diff_expression_results_E_2_15 = diff_expression_results
diff_expression_results_E_8_1 = diff_expression_results
diff_expression_results_E_6_8 = diff_expression_results
diff_expression_results_E_2_8 = diff_expression_results

#PREPARING gsea input
list_of_contrasts_step2 = list(diff_expression_results_E_1_13,diff_expression_results_E_15_1,diff_expression_results_E_6_15,diff_expression_results_E_2_15,
                               diff_expression_results_E_8_1,diff_expression_results_E_6_8,diff_expression_results_E_2_8)
names(list_of_contrasts_step2) = c("diff_expression_results_E_1_13","diff_expression_results_E_15_1","diff_expression_results_E_6_15","diff_expression_results_E_2_15",
                                   "diff_expression_results_E_8_1","diff_expression_results_E_6_8","diff_expression_results_E_2_8")

OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/GO_step2/"

#Change names to entrez code
name_changer = function(especie_database, especie_name, x){
  ensembl <- useMart(especie_database, dataset = especie_name)
  annotation = getBM(attributes=c("hgnc_symbol","entrezgene_id"), values=x, mart=ensembl)
  return(annotation)
}

library("enrichplot")
library("clusterProfiler")
library("org.Hs.eg.db")
list_of_gsea = list()
for (n in 1:length(list_of_contrasts_step2)){
  print(n)
  #Create the geneList
  geneList <- list_of_contrasts_step2[[n]]$logFC
  names(geneList) <- as.character(rownames(list_of_contrasts_step2[[n]]))
  geneList = sort(geneList, decreasing = TRUE)
   
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]
  
  #GSEA
  ego <- gseGO(geneList    = geneList,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 100,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = FALSE)
  ego_table <- as.data.frame(ego) 
  write.csv(ego_table,paste0(OutPath,"gsea_","BP","_contrast_",names(list_of_contrasts_step2[n]),".csv"))
  list_of_gsea = append(list_of_gsea, list(ego_table))
  #plots
  if (nrow(ego_table) > 1){
    pdf(file=paste0(OutPath,"gsea_dotplot_","BP","_contrast_",names(list_of_contrasts_step2[n]),".pdf"), width=8, height=12)  
    plot(dotplot(ego, showCategory=30))
    dev.off()
    pdf(file=paste0(OutPath,"gsea_gseaPlot_","BP","_contrast_",names(list_of_contrasts_step2[n]),".pdf"), width=12, height=12)  
    plot(gseaplot2(ego, geneSetID = 1:3))
    dev.off()    
  }
}

#Preparing genesetcluster for step2
names(list_of_gsea) = c("diff_expression_results_E_1_13","diff_expression_results_E_15_1","diff_expression_results_E_6_15","diff_expression_results_E_2_15",
                                   "diff_expression_results_E_8_1","diff_expression_results_E_6_8","diff_expression_results_E_2_8")

# STEP 1: with all samples included
pathways_markers_E_1_13 = list_of_gsea[["diff_expression_results_E_1_13"]]
pathways_markers_E_15_1 = list_of_gsea[["diff_expression_results_E_15_1"]]
pathways_markers_E_6_15 = list_of_gsea[["diff_expression_results_E_6_15"]]
pathways_markers_E_2_15 = list_of_gsea[["diff_expression_results_E_2_15"]]
pathways_markers_E_8_1 = list_of_gsea[["diff_expression_results_E_8_1"]]
pathways_markers_E_6_8 = list_of_gsea[["diff_expression_results_E_6_8"]]
pathways_markers_E_2_8 = list_of_gsea[["diff_expression_results_E_2_8"]]

#get significant pathways
pathways_markers_E_1_13_sig<-pathways_markers_E_1_13[pathways_markers_E_1_13$p.adjust<0.1,]
pathways_markers_E_15_1_sig<-pathways_markers_E_15_1[pathways_markers_E_15_1$p.adjust<0.1,]
pathways_markers_E_6_15_sig<-pathways_markers_E_6_15[pathways_markers_E_6_15$p.adjust<0.1,]
pathways_markers_E_2_15_sig<-pathways_markers_E_2_15[pathways_markers_E_2_15$p.adjust<0.1,]
pathways_markers_E_8_1_sig<-pathways_markers_E_8_1[pathways_markers_E_8_1$p.adjust<0.1,]
pathways_markers_E_6_8_sig<-pathways_markers_E_6_8[pathways_markers_E_6_8$p.adjust<0.1,]
pathways_markers_E_2_8_sig<-pathways_markers_E_2_8[pathways_markers_E_2_8$p.adjust<0.1,]

# STEP 2: CREATE PATHWAY OBJECT
object <- ObjectCreator(Pathways = c(as.character(pathways_markers_E_1_13_sig$ID),
                                     as.character(pathways_markers_E_15_1_sig$ID),
                                     as.character(pathways_markers_E_6_15_sig$ID),
                                     as.character(pathways_markers_E_2_15_sig$ID),
                                     as.character(pathways_markers_E_8_1_sig$ID),
                                     as.character(pathways_markers_E_6_8_sig$ID),
                                     as.character(pathways_markers_E_2_8_sig$ID)),
                        Molecules =c(as.character(pathways_markers_E_1_13_sig$core_enrichment),
                                     as.character(pathways_markers_E_15_1_sig$core_enrichment),
                                     as.character(pathways_markers_E_6_15_sig$core_enrichment),
                                     as.character(pathways_markers_E_2_15_sig$core_enrichment),
                                     as.character(pathways_markers_E_8_1_sig$core_enrichment),
                                     as.character(pathways_markers_E_6_8_sig$core_enrichment),
                                     as.character(pathways_markers_E_2_8_sig$core_enrichment)),
                        Groups=c(rep("E_1_13",times=nrow(pathways_markers_E_1_13_sig)),
                                 rep("E_15_1",times=nrow(pathways_markers_E_15_1_sig)),
                                 rep("E_6_15",times=nrow(pathways_markers_E_6_15_sig)),
                                 rep("E_2_15",times=nrow(pathways_markers_E_2_15_sig)),
                                 rep("E_8_1",times=nrow(pathways_markers_E_8_1_sig)),
                                 rep("E_6_8",times=nrow(pathways_markers_E_6_8_sig)),
                                 rep("E_2_8",times=nrow(pathways_markers_E_2_8_sig))),
                        structure = "SYMBOL", Type="", sep="/", Source="GSEA", organism="org.Hs.eg.db")
                        
# STEP3: COMBINE
combine<- CombineGeneSets(object, display="Expanded")

combine2 <- combine
combine2@metadata$Groups <- c("E_1_13","E_15_1","E_2_15","E_8_1","E_6_8","E_2_8")
combine2 <- ClusterGeneSets(combine2, clusters = 8, method = "kmeans",order = "cluster")

x <- combine2@plot$aka2
x <- x[,!colnames(x) == "Group"]
combine2@plot$aka2 <- x

x <- combine2@plot$aka3
x <- x[2:length(x)]
combine2@plot$aka3 <- x

png("genesetcluster_lvl2_T_neurons_8c.png", width = 1200, height = 1500)
PlotGeneSets(Object = combine2, main = "GO_BP_Neurons", RR.max = 20, annotation.mol=F, legend = T)
dev.off()

# STEP4: Export tables
clusters<-combine2@Data[[1]]
write.table(clusters, file="Clusters_pathways_descriptions_combine_lvl2_T_neurons_8c_.csv", sep=";", row.names = FALSE)

annotation_pathways = bind_rows(list_of_gsea, .id = rownames(list_of_gsea))

all<-merge(clusters, annotation_pathways, by.x="Pathways", by.y="ID", all=FALSE)
write.table(all, file="Clusters_pathways_descriptions_combine_lvl2_T_neurons_8c_.csv", sep=";", row.names = FALSE)


######################################################################################################################################################
# STEP 3 - create matrix with pertainti of genes on each category
potential_markers_E_1_13 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "1", ident.2 = "13", only.pos = FALSE)
potential_markers_E_15_1 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "15", ident.2 = "1", only.pos = FALSE)
potential_markers_E_6_15 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "6", ident.2 = "15", only.pos = FALSE)
potential_markers_E_2_15 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "2", ident.2 = "15", only.pos = FALSE)
potential_markers_E_8_1 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "8", ident.2 = "1", only.pos = FALSE)
potential_markers_E_6_8 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "6", ident.2 = "8", only.pos = FALSE)
potential_markers_E_2_8 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "2", ident.2 = "8", only.pos = FALSE)

potential_markers_E_1_13 = potential_markers_E_1_13[abs(potential_markers_E_1_13$p_val_adj) < 0.05, ]
potential_markers_E_15_1 = potential_markers_E_15_1[abs(potential_markers_E_15_1$p_val_adj) < 0.05, ]
potential_markers_E_6_15 = potential_markers_E_6_15[abs(potential_markers_E_6_15$p_val_adj) < 0.05, ]
potential_markers_E_2_15 = potential_markers_E_2_15[abs(potential_markers_E_2_15$p_val_adj) < 0.05, ]
potential_markers_E_8_1 = potential_markers_E_8_1[abs(potential_markers_E_8_1$p_val_adj) < 0.05, ]
potential_markers_E_6_8 = potential_markers_E_6_8[abs(potential_markers_E_6_8$p_val_adj) < 0.05, ]
potential_markers_E_2_8 = potential_markers_E_2_8[abs(potential_markers_E_2_8$p_val_adj) < 0.05, ]

write.table(potential_markers_E_2_8, file="Potential_markers_step_1_E_2_8.csv", sep=";", row.names = TRUE)

names_sig_T_potential_E_1_13 = rownames(potential_genes_E_1_13)
names_sig_T_potential_E_15_1 = rownames(potential_genes_E_15_1)
names_sig_T_potential_E_6_15 = rownames(potential_genes_E_6_15)
names_sig_T_potential_E_2_15 = rownames(potential_genes_E_2_15)
names_sig_T_potential_E_8_1 = rownames(potential_genes_E_8_1)
names_sig_T_potential_E_6_8 = rownames(potential_genes_E_6_8)
names_sig_T_potential_E_2_8 = rownames(potential_genes_E_2_8)

names_sig_C_markers_E_1_13 = rownames(potential_markers_E_1_13)
names_sig_C_markers_E_15_1 = rownames(potential_markers_E_15_1)
names_sig_C_markers_E_6_15 = rownames(potential_markers_E_6_15)
names_sig_C_markers_E_2_15 = rownames(potential_markers_E_2_15)
names_sig_C_markers_E_8_1 = rownames(potential_markers_E_8_1)
names_sig_C_markers_E_6_8 = rownames(potential_markers_E_6_8)
names_sig_C_markers_E_2_8 = rownames(potential_markers_E_2_8)

out <- table(stack(mget(ls(pattern = "^names_sig_")[c(4:10,15:21)])))
out <- out[, c("names_sig_C_markers_E_1_13", "names_sig_C_markers_E_15_1", "names_sig_C_markers_E_6_15", "names_sig_C_markers_E_2_15", "names_sig_C_markers_E_8_1", "names_sig_C_markers_E_6_8", "names_sig_C_markers_E_2_8", "names_sig_T_potential_E_1_13","names_sig_T_potential_E_15_1","names_sig_T_potential_E_6_15","names_sig_T_potential_E_2_15","names_sig_T_potential_E_8_1",
"names_sig_T_potential_E_6_8","names_sig_T_potential_E_2_8")]

out_table = data.frame(matrix(data = out, ncol = ncol(out), nrow = nrow(out)))
rownames(out_table) = rownames(out)
colnames(out_table) = colnames(out)
out_table$Events = rowSums(out_table)
write.table(out_table, file="genes_of_interest_per_contrast_lvl2.csv", sep=";", row.names = TRUE)

out_summary = t(out) %*% out
write.table(out_summary, file="genes_of_interest_per_contrast_lvl2_summary.csv", sep=";", row.names = TRUE)

#GENES IN BOTH CONTROL + TREATMENT
lista_gene_names=Reduce(intersect, list(names_sig_T_potential_E_1_13,names_sig_C_markers_E_6_15))

potential_markers_E_6_15[lista_gene_names,]
potential_genes_E_1_13[lista_gene_names,]

###########################################
#Auto generate all tables for each T contrast vs all C 
aa = list(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8)
aa_2 = list(potential_genes_E_1_13,potential_genes_E_15_1,potential_genes_E_6_15,potential_genes_E_2_15,potential_genes_E_8_1,potential_genes_E_6_8,potential_genes_E_2_8)

bb = list(names_sig_C_markers_E_1_13,names_sig_C_markers_E_15_1,names_sig_C_markers_E_6_15,names_sig_C_markers_E_2_15,names_sig_C_markers_E_8_1,names_sig_C_markers_E_6_8,names_sig_C_markers_E_2_8)
bb_2 = list(potential_markers_E_1_13,potential_markers_E_15_1,potential_markers_E_6_15,potential_markers_E_2_15,potential_markers_E_8_1,potential_markers_E_6_8,potential_markers_E_2_8)

lista_all = list()
for (n in (1:length(aa))){
  lista_semi = list()
  for (m in (1:length(bb))){
    lista_gene_names=Reduce(intersect, list(aa[[n]],bb[[m]]))
    control = bb_2[[m]][lista_gene_names,]
    control = control[,c(2,5)]
    tratamiento = aa_2[[n]][lista_gene_names,]
    tratamiento = tratamiento[,c(1,5)]
    union = data.frame(logfc_c = control[,1],
                       logfc_t = tratamiento[,1],
                       adj.pval_c = control[,2],
                       adj.pval_t = tratamiento[,2])
    rownames(union) = rownames(control)
    lista_semi = append(lista_semi, list(union))
  }
  lista_all = append(lista_all,list(lista_semi))
}
###########################################

#GENES ONLY IN CONTROL
lista_gene_names = names_sig_C_markers_E_2_8[!(names_sig_C_markers_E_2_8 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
potential_markers_E_2_8[lista_gene_names,][,c(2,5)]
#GENES ONLY IN TREATMENT
lista_gene_names = names_sig_T_potential_E_2_8[!(names_sig_T_potential_E_2_8 %in% c(names_sig_C_markers_E_1_13,names_sig_C_markers_E_15_1,names_sig_C_markers_E_6_15,names_sig_C_markers_E_2_15,names_sig_C_markers_E_8_1,names_sig_C_markers_E_6_8,names_sig_C_markers_E_2_8))]
potential_genes_E_2_8[lista_gene_names,][,c(1,5)]


#For heatmap rownames colors
#BOTH
both_conditions_all = c()
for (n in (1:length(lista_all))){
  for (m in (1:length(lista_all[[n]]))){
    both_conditions_all = append(both_conditions_all, rownames(lista_all[[n]][[m]]))
  }
}
both_conditions_all = unique(both_conditions_all)

#only T
only_T = c()
for (n in (1:length(bb))){
  lista_gene_names = aa[[n]][!(aa[[n]] %in% c(names_sig_C_markers_E_1_13,names_sig_C_markers_E_15_1,names_sig_C_markers_E_6_15,names_sig_C_markers_E_2_15,names_sig_C_markers_E_8_1,names_sig_C_markers_E_6_8,names_sig_C_markers_E_2_8))]
  only_T = append(only_T, lista_gene_names)
}

only_T_conditions_all = unique(only_T)

#Only C
only_C_conditions_all_E_1_13 = names_sig_C_markers_E_1_13[!(names_sig_C_markers_E_1_13 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
only_C_conditions_all_E_15_1 = names_sig_C_markers_E_15_1[!(names_sig_C_markers_E_15_1 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
only_C_conditions_all_E_6_15 = names_sig_C_markers_E_6_15[!(names_sig_C_markers_E_6_15 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
only_C_conditions_all_E_2_15 = names_sig_C_markers_E_2_15[!(names_sig_C_markers_E_2_15 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
only_C_conditions_all_E_8_1 = names_sig_C_markers_E_8_1[!(names_sig_C_markers_E_8_1 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
only_C_conditions_all_E_6_8 = names_sig_C_markers_E_6_8[!(names_sig_C_markers_E_6_8 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]
only_C_conditions_all_E_2_8 = names_sig_C_markers_E_2_8[!(names_sig_C_markers_E_2_8 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_15_1,names_sig_T_potential_E_6_15,names_sig_T_potential_E_2_15,names_sig_T_potential_E_8_1,names_sig_T_potential_E_6_8,names_sig_T_potential_E_2_8))]

only_C_condition_all = unique(c(only_C_conditions_all_E_1_13,only_C_conditions_all_E_15_1,only_C_conditions_all_E_6_15,only_C_conditions_all_E_2_15,only_C_conditions_all_E_8_1,only_C_conditions_all_E_6_8,only_C_conditions_all_E_2_8))

######################################################################################################################################################
#PLOT RESULTS FROM LVL 2 NEURONS
######################################################################################################################################################
#HEATMAP PLOT OF BINARIZED AND EXPRESION CASCADES OF ACTIVATION FOR ALL GENES SIGNIFICANT CHANGING
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize) # for the colorRamp2() function

#Table of binarized cascade activation + event(sum of cascade activations)
out_table_plot_lvl2_neurons = out_table
out_table_plot_lvl2_neurons$Events = NULL

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 7)
tiempos=append(tiempos,rep("Treatment", 7))

#replace colnames by nice ones
table_plot_names =c("C_E_1_13","C_E_15_1","C_E_6_15","C_E_2_15","C_E_8_1","C_E_6_8","C_E_2_8","T_E_1_13","T_E_15_1","T_E_6_15","T_E_2_15","T_E_8_1","T_E_6_8","T_E_2_8")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8","E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("E_1_13" = "#c51b7d",
                                           "E_15_1" = "#e9a3c9",
                                           "E_6_15" = "#fde0ef",
                                           "E_2_15" = "#f7f7f7",
                                           "E_8_1" = "#e6f5d0",
                                           "E_6_8" = "#a1d76a",
                                           "E_2_8" = "#4d9221")
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_lvl2_neurons
colnames(my_data) = rownames(meta_data_plot)
#Binarized version of table and plot
# Heatmap
col_fun = colorRamp2(c(0, 1), c("#f7f7f7", "#998ec3"))

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))


pdf(file = "Binary_cascade_of_activation_transition_lvl_2_neurons.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

######################################################################################################################################################
#REPEAT HEATMAP WITH EXPRESSION VALUES
#Table of binarized cascade activation + event(sum of cascade activations)
out_table_plot_lvl2_neurons = out_table
out_table_plot_lvl2_neurons$Events = NULL

#Fom binary to logFC matrix
out_table_plot_lvl2_neurons[,1] = replace(out_table_plot_lvl2_neurons[,1], out_table_plot_lvl2_neurons[,1]==1, potential_markers_E_1_13[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,1]==1],]$avg_log2FC)
out_table_plot_lvl2_neurons[,2] = replace(out_table_plot_lvl2_neurons[,2], out_table_plot_lvl2_neurons[,2]==1, potential_markers_E_15_1[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,2]==1],]$avg_log2FC)
out_table_plot_lvl2_neurons[,3] = replace(out_table_plot_lvl2_neurons[,3], out_table_plot_lvl2_neurons[,3]==1, potential_markers_E_6_15[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,3]==1],]$avg_log2FC)
out_table_plot_lvl2_neurons[,4] = replace(out_table_plot_lvl2_neurons[,4], out_table_plot_lvl2_neurons[,4]==1, potential_markers_E_2_15[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,4]==1],]$avg_log2FC)
out_table_plot_lvl2_neurons[,5] = replace(out_table_plot_lvl2_neurons[,5], out_table_plot_lvl2_neurons[,5]==1, potential_markers_E_8_1[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,5]==1],]$avg_log2FC)
out_table_plot_lvl2_neurons[,6] = replace(out_table_plot_lvl2_neurons[,6], out_table_plot_lvl2_neurons[,6]==1, potential_markers_E_6_8[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,6]==1],]$avg_log2FC)
out_table_plot_lvl2_neurons[,7] = replace(out_table_plot_lvl2_neurons[,7], out_table_plot_lvl2_neurons[,7]==1, potential_markers_E_2_8[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,7]==1],]$avg_log2FC)

out_table_plot_lvl2_neurons[,8] = replace(out_table_plot_lvl2_neurons[,8], out_table_plot_lvl2_neurons[,8]==1, potential_genes_E_1_13[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,8]==1],]$logFC)
out_table_plot_lvl2_neurons[,9] = replace(out_table_plot_lvl2_neurons[,9], out_table_plot_lvl2_neurons[,9]==1, potential_genes_E_15_1[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,9]==1],]$logFC)
out_table_plot_lvl2_neurons[,10] = replace(out_table_plot_lvl2_neurons[,10], out_table_plot_lvl2_neurons[,10]==1, potential_genes_E_6_15[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,10]==1],]$logFC)
out_table_plot_lvl2_neurons[,11] = replace(out_table_plot_lvl2_neurons[,11], out_table_plot_lvl2_neurons[,11]==1, potential_genes_E_2_15[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,11]==1],]$logFC)
out_table_plot_lvl2_neurons[,12] = replace(out_table_plot_lvl2_neurons[,12], out_table_plot_lvl2_neurons[,12]==1, potential_genes_E_8_1[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,12]==1],]$logFC)
out_table_plot_lvl2_neurons[,13] = replace(out_table_plot_lvl2_neurons[,13], out_table_plot_lvl2_neurons[,13]==1, potential_genes_E_6_8[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,13]==1],]$logFC)
out_table_plot_lvl2_neurons[,14] = replace(out_table_plot_lvl2_neurons[,14], out_table_plot_lvl2_neurons[,14]==1, potential_genes_E_2_8[rownames(out_table_plot_lvl2_neurons)[out_table_plot_lvl2_neurons[,14]==1],]$logFC)

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 7)
tiempos=append(tiempos,rep("Treatment", 7))

#replace colnames by nice ones
table_plot_names =c("C_E_1_13","C_E_15_1","C_E_6_15","C_E_2_15","C_E_8_1","C_E_6_8","C_E_2_8","T_E_1_13","T_E_15_1","T_E_6_15","T_E_2_15","T_E_8_1","T_E_6_8","T_E_2_8")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8","E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("E_1_13" = "#c51b7d",
                                           "E_15_1" = "#e9a3c9",
                                           "E_6_15" = "#fde0ef",
                                           "E_2_15" = "#f7f7f7",
                                           "E_8_1" = "#e6f5d0",
                                           "E_6_8" = "#a1d76a",
                                           "E_2_8" = "#4d9221")
                          )
                          
ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_lvl2_neurons
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_2_neurons.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

col_fun = colorRamp2(c(my_data_min/100, 0, my_data_max/100), c("blue", "white", "red"))
pdf(file = "Expression_2 colors_cascade_of_activation_transition_lvl_2_neurons.pdf",
    width = 10,
    height = 110)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment", at = c(-1, 0, 1), labels = c("Sig.Down.regulated", "No.significant", "Sig.Up.regulated")),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

#GET 285 TF from JASPAR2020 285
ocr_tf_matrix=read.table("ocr_tf_matrix.txt", sep = "")
names(ocr_tf_matrix) 

library(JASPAR2022)
getJasparMotifs = function(species = "Homo sapiens", collection = "CORE"){
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- collection
    out <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
    if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
        names(out) <- paste(names(out), TFBSTools::name(out),
            sep = "_")
    return(out)
}
pfm <- getJasparMotifs()
sub(".*_", "", names(pfm))

#191 in common from 285 vs 692
sum(names(ocr_tf_matrix) %in% sub(".*_", "", names(pfm)))

#How many TF do we have in out 1728? 
sum(rownames(my_data) %in% names(ocr_tf_matrix)) #32
sum(rownames(my_data) %in% sub(".*_", "", names(pfm))) #67 USE JASPAR 2022
tf_in_our_contrasts= rownames(my_data)[rownames(my_data) %in% sub(".*_", "", names(pfm))]

write.table(out_table[tf_in_our_contrasts,], file="TF_binary_of_interest_per_contrast_lvl_2_neurons.csv", sep=";", row.names = TRUE)
write.table(out_table_plot_lvl2_neurons[tf_in_our_contrasts,], file="TF_of_interest_per_contrast_lvl_2_neurons.csv", sep=";", row.names = TRUE)

##################################################################################################################################################
##################################################################################################################################################
#HEATMAPS SUMMARY 1-c+t+only_T AND 2-only_C
genelist_both_and_T_only = c(both_conditions_all,only_T_conditions_all)

out_table_plot_lvl2_neurons_genelist_both_and_T_only= out_table_plot_lvl2_neurons[genelist_both_and_T_only,]
# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 7)
tiempos=append(tiempos,rep("Treatment", 7))

#replace colnames by nice ones
table_plot_names =c("C_E_1_13","C_E_15_1","C_E_6_15","C_E_2_15","C_E_8_1","C_E_6_8","C_E_2_8","T_E_1_13","T_E_15_1","T_E_6_15","T_E_2_15","T_E_8_1","T_E_6_8","T_E_2_8")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8","E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8")


annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("E_1_13" = "#c51b7d",
                                           "E_15_1" = "#e9a3c9",
                                           "E_6_15" = "#fde0ef",
                                           "E_2_15" = "#f7f7f7",
                                           "E_8_1" = "#e6f5d0",
                                           "E_6_8" = "#a1d76a",
                                           "E_2_8" = "#4d9221")
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_lvl2_neurons_genelist_both_and_T_only
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('black', nrow(my_data))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_2_neurons_genelist_both_and_T_only.pdf",
    width = 10,
    height = 16)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

##########
#2-only_C

#Order base on logFC
temp = potential_markers_E_1_13[order(potential_markers_E_1_13$avg_log2FC),]
temp_E_1_13 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_15_1[order(potential_markers_E_15_1$avg_log2FC),]
temp_E_15_1 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_6_15[order(potential_markers_E_6_15$avg_log2FC),]
temp_E_6_15 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_2_15[order(potential_markers_E_2_15$avg_log2FC),]
temp_E_2_15 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_8_1[order(potential_markers_E_8_1$avg_log2FC),]
temp_E_8_1 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_6_8[order(potential_markers_E_6_8$avg_log2FC),]
temp_E_6_8 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_2_8[order(potential_markers_E_2_8$avg_log2FC),]
temp_E_2_8 = c(rownames(head(temp,30)), rownames(tail(temp,30)))

temp_all = unique(c(temp_E_1_13,temp_E_15_1,temp_E_6_15,temp_E_2_15,temp_E_8_1,temp_E_6_8,temp_E_2_8)) #206

only_C_condition_all_subset = only_C_condition_all[only_C_condition_all %in% temp_all]
out_table_plot_lvl2_neurons_genelist_C_only = out_table_plot_lvl2_neurons[only_C_condition_all_subset,]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 7)
tiempos=append(tiempos,rep("Treatment", 7))

#replace colnames by nice ones
table_plot_names =c("C_E_1_13","C_E_15_1","C_E_6_15","C_E_2_15","C_E_8_1","C_E_6_8","C_E_2_8","T_E_1_13","T_E_15_1","T_E_6_15","T_E_2_15","T_E_8_1","T_E_6_8","T_E_2_8")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8","E_1_13","E_15_1","E_6_15","E_2_15","E_8_1","E_6_8","E_2_8")


annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("E_1_13" = "#c51b7d",
                                           "E_15_1" = "#e9a3c9",
                                           "E_6_15" = "#fde0ef",
                                           "E_2_15" = "#f7f7f7",
                                           "E_8_1" = "#e6f5d0",
                                           "E_6_8" = "#a1d76a",
                                           "E_2_8" = "#4d9221")
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# ANNOTATION TABLE
my_data <- out_table_plot_lvl2_neurons_genelist_C_only
colnames(my_data) = rownames(meta_data_plot)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_C_condition_all_subset
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('black', nrow(my_data))
fontcolors[row_idx] <- 'd8b365'
fontcolors[row_idx_t] <- '#d8b365'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

#Gradient for EXPESSION data
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "Expression_cascade_of_activation_transition_lvl_2_neurons_genelist_C_only.pdf",
    width = 10,
    height = 10)
Heatmap(
  my_data,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(my_data),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()


##########
#GO FROM OUTPUT GROUPS FROM HEATMAPS BASE ON SAME TENDENCIES IN NEUROGENESIS
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/GO_outcomes/"

both_and_t_group_1_lvl_1 =c("NEFM","NEFL","STMN2","THSD7A","POU2F2","SCG2","UCHL1","ONECUT2","SNCG","NSG1","KIF5C","TUBA1A","TUBB",
                            "RBP1","NCAM1","TUBB2B","TUBB4A","GUK1","FAM57B","NUDT14","DDAH2","EBF1","ZNF428","DCTN3","SEC62","NDUFB1","PRPF40A","IDH2")               
both_and_t_group_2_lvl_1 =c("APC","NOL7","SRP14","CLUAP1","UBE2E3","STRAP","ENAH","ROMO1","SNRPF","TMEM59","SLIRP","SEC61G","SEC61G.1")
both_and_t_group_3_lvl_1 =c("PCDH17","MAGED2","SKA2","MYL6","RPL37","B2M","MYL12A","CDC34","EIF4H","COPZ1","POLR1D","SIPA1L2","PPM1K","TCF7L1","SLC2A1","MIR99AHG","TIMP1","BCAN","S100B","ARL6IP5","TMEM161B-AS1","ANXA2")
both_and_t_group_4_lvl_1 =c("LIMCH1","PTN","FOS","OPRK1","ID1","ID3","CLU","SPARC","CD99","CTGF")
both_and_t_group_5_lvl_1 =c("CRABP1","HES6","SOX11","RND3","GLRX","KLHL35")
both_and_t_group_6_lvl_1 =c("COTL1","PTMA","CDC25B","PPP1R14C","DLL1","NIN","SPSB4","PRMT8","ASCL1","DIO3","GLUL","PRSS23")
both_and_t_group_7_lvl_1 =c("HOXB2","DRAXIN","IER2","CBLB","CRYBA1","SCRG1","SH3BGRL3","LMO1","IGDCC3","ZBTB16","SFRP1","GNG5","RPS27L","DLL3","CCND1","C1orf61")
both_and_t_group_8_lvl_1 =c("NCALD","RGMA","NES","HES5","IGFBP5")

both_C_and_T_lvl_1_group_sets = list(both_and_t_group_1_lvl_1,both_and_t_group_2_lvl_1,both_and_t_group_3_lvl_1,both_and_t_group_4_lvl_1,both_and_t_group_5_lvl_1,both_and_t_group_6_lvl_1,both_and_t_group_7_lvl_1,both_and_t_group_8_lvl_1)
names(both_C_and_T_lvl_1_group_sets) = c("both_and_t_group_1_lvl_1","both_and_t_group_2_lvl_1","both_and_t_group_3_lvl_1","both_and_t_group_4_lvl_1","both_and_t_group_5_lvl_1","both_and_t_group_6_lvl_1","both_and_t_group_7_lvl_1","both_and_t_group_8_lvl_1")


only_C_group_1_lvl_1 =c("CDKN1C","GADD45G","JAG1")
only_C_group_2_lvl_1 =c("SYT4","GAP43","ANK3","PCDH9","ZFHX4","MAPT","CELF4","TUBB2A","TTC9B","NRXN1","TMSB10","PCSK1N","CSRNP3","ZFHX3","INA","MEIS2","MAP1B","SCG3","STMN4","MLLT11",
"STMN1","CD24","C1QL1","CALM1","RUNX1T1","SHOX2","NTM","CALB1","TMX4","FEV","SCN2A","LINC01828","SH3BP5","PCLO")
only_C_group_3_lvl_1 =c("SMS","CNTNAP2","TUBA1B","RPL37A","CENPF","PHGDH","GYPC","EIF4EBP1","HMGA1","ISYNA1","RPS27A","EEF1B2","TKT")
only_C_group_4_lvl_1 =c("EPB41","RGMB","TBCA","ARL4D","NEUROD4","TCF12","EZR","NEAT1","RPS7","SYT1","DCX","ELAVL2","SOX4","PCBP4","CBFA2T2","KCNQ1OT1","TAGLN3","MIAT","PLXNA2","BTG1","TCF4","NHLH1","ENC1","MAGI1",
"NEUROG2","PCDH18","BTG2","TFDP2")
only_C_group_5_lvl_1 =c("ID4","SPECC1","SOX9","QKI","PTPRZ1","FBXO32","SMOC1","ZFP36L1","AL139246.5","TFPI2","PMP2","ABAT","FAM213A","CYP26B1","WNT7B","FABP5","GNG12","LIX1","ENPP2","SLC1A3","PLP1","FGFBP3","FRAS1",
"LAMA4","TPM1","FABP7","IFITM3","IDH1","HES4","S100A13","CYR61","CXCR4","METRN","LINC00461","ANXA5","SOX3","PON2","PAX6","IGFBP2","GM6B","ITGB8","NTRK2","TTYH1","SFRP2","VIM")

only_C_lvl_1_group_sets = list(only_C_group_1_lvl_1,only_C_group_2_lvl_1,only_C_group_3_lvl_1,only_C_group_4_lvl_1,only_C_group_5_lvl_1)
names(only_C_lvl_1_group_sets) = c("only_C_group_1_lvl_1","only_C_group_2_lvl_1","only_C_group_3_lvl_1","only_C_group_4_lvl_1","only_C_group_5_lvl_1")


#######
#GO - BP
#######
#Create the geneList
for (n in 1:length(only_C_lvl_1_group_sets)){
  geneList <- only_C_lvl_1_group_sets[[n]]
  ego <- enrichGO(gene       = geneList,
             universe      = rownames(combined_seurat),
             OrgDb         = org.Hs.eg.db,
             keyType       = 'SYMBOL',
             ont           = 'BP',
             pAdjustMethod = "BH",
             minGSSize     = 10,
             maxGSSize     = 500,
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.1)
  ego_table <- as.data.frame(ego) 
  write.csv(ego_table,paste0(OutPath,"_",names(only_C_lvl_1_group_sets[n]),".csv"))
}

geneList

both_conditions_all
only_T_conditions_all
only_C_condition_al

#plots
if (nrow(ego_table) > 0){
  #ego <- pairwise_termsim(ego) 
  pdf(file=paste0(OutPath,"ego_barplot_BP_cluster_group_2.pdf"), width=8, height=8)  
  plot(barplot(ego, showCategory=20))  
  dev.off()
  pdf(file=paste0(OutPath,"egocneplot__BP_cluster_group_2.pdf"), width=8, height=8)  
  plot(cnetplot(ego, categorySize="pvalue", foldChange=geneList))
  dev.off()  
} 

#######
#DO
#######
library(DOSE)

annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', rownames(combined_seurat))
entrezgene = annotation$entrezgene_id[match(rownames(combined_seurat), annotation$hgnc_symbol)]
universe_genes = as.character(entrezgene)
universe_genes = universe_genes[!is.na(universe_genes)]

#Create the geneList
for (n in 1:length(both_C_and_T_lvl_1_group_sets)){
  geneList <- both_C_and_T_lvl_1_group_sets[[n]]
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', geneList)
  entrezgene = annotation$entrezgene_id[match(geneList, annotation$hgnc_symbol)]
  geneList = as.character(entrezgene)
  
  x <- enrichDO(gene          = geneList,
                ont           = "DO",
                pAdjustMethod = "BH",
                universe      = universe_genes,
                minGSSize     = 5,
                maxGSSize     = 500,
                pvalueCutoff  = 0.5,
                qvalueCutoff  = 1,
                readable      = FALSE)
  x_table <- as.data.frame(x) 
  write.csv(x_table,paste0(OutPath,"_",names(both_C_and_T_lvl_1_group_sets[n]),"DO_.csv"))
}

write.csv(x_table,paste0(OutPath,"_both_C_and_T_lvl_2_all_genes_DO.csv.csv"))

#######
#KEGG
#######
for (n in 1:length(both_C_and_T_lvl_1_group_sets)){
  geneList <- both_C_and_T_lvl_1_group_sets[[n]]
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', geneList)
  entrezgene = annotation$entrezgene_id[match(geneList, annotation$hgnc_symbol)]
  geneList = as.character(entrezgene)
  
  kk <- enrichKEGG(gene          = geneList,
                   keyType       = "kegg",
                   organism      = 'hsa',
                   universe      = universe_genes, 
                   minGSSize     = 10, 
                   maxGSSize     = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)
  x_table <- as.data.frame(kk) 
  write.csv(x_table, paste0(OutPath,"_",names(both_C_and_T_lvl_1_group_sets[n]),"KEGG_.csv"))
}

geneList = both_conditions_all

both_conditions_all
only_T_conditions_all
only_C_condition_all

#######
#PATHVIEW -KEGG
#######
library("pathview")

hsa05224 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05224",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))






#######
#CHECK GENES TO TRY YO FIT TRAYECTORY ON LVL2
#######
#bueno "GFAP","PROX1","NEUROG2","NES","SOX11","FAM114A1","SOX4","PAX6","REST","MAP2","NCAM1","HES5","MCM2","CALB2","DCX","NEFH","STMN1","SYP","SOX2","RBFOX3","NEU1","CALB1","PCNA","CBS"
pdf(file="CHECK_some_genes_validate_directions_dots_caracterization_4only_control_cells_subset.pdf", width=10, height=6)
features = c("PCNA","SOX2","GFAP","PROX1","NEUROG2","NES","SOX11","FAM114A1","SOX4","PAX6","REST","MAP2","MAP2A","MAP2B","MAP2C","NCAM1","HES5","MCM2","CALB2","DCX","NEFH","STMN1","SYP","RBFOX3","NEU1","CALB1","CBS")
DotPlot(only_control_cells_subset, features = features, cols = c("blue", "yellow")) + RotatedAxis()  #only_control_cells_subset_neurons combined_seurat
dev.off()

Idents(only_control_cells_subset_neurons) <- factor(Idents(only_control_cells_subset_neurons), levels= c(13,15,1,2,6,8))

"GFAP","SOX11","PROX1","NEUROG2","NES","FAM114A1","PROX1","SOX4","PAX6","REST","MAP2","NCAM1","HES5","NES","MCM2","CALB2","DCX","NEFH","MAP2","MAP2A","MAP2B","MAP2C","STMN1","SYP","SOX2","NEUN","CALB1"

pdf(file="CHECK_some_genes_validate_directions_umaps.pdf", width=12, height=12)
p1 = DimPlot(combined_seurat, reduction = 'umap', group.by = "seurat_clusters",label = TRUE)
p2 = FeaturePlot(
      object = combined_seurat,
      features = c("CRABP1","NEFM","NHLH1"),
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 1
      )
CombinePlots(list(p1,p2))
dev.off()

Idents(combined_seurat) = combined_seurat@meta.data$integrated_snn_res.0.6 #contrasts_groups

#POSIBLE EARLY
"RACK1","BMP7","ISG15","TENT5A",

#POSIBLE LATE
"MAP2","CDK5R1","SEMA3A","DRAXIN","PTEN","DLX1","RTN4","DLX2","SEMA6A","SLIT1","DCC","MBP"

#from cluster 13
"CRABP1","NEFM","NHLH1","IGDCC3","KLHL35","ST18","TUBA1A","GLRX","ACTG1","TUBB2B","RPL29","RGMB","TUBB"

#from clsuter 1
"LINC00461","MALAT1","FIGN","C1orf61","FOXP2","GPM6A","TSHZ2","RUNX1T1","MARCKS","SOX4","CDH8","LHX5-AS1","EPHA4","FAT3","NR2F2","ANK3","PAX6","IRX5","FABP7","NOVA1","VCAN","MAP2","SESN3"

#from cluster 8
"FEV","GATA2","PEG10","GATA3","CRYBA2","GATA3-AS1","CALB1","SLC17A8","LMX1B","SOX2","DNAJC12","SSTR2","CPNE4","CACNA2D11","ANKS1B","FAM213A1","CASD1"












#all genes positive word go
positive=unique(unlist(strsplit(temporal_go[temporal_go$Description %like% "positive", ]$geneID, "/")))
#all genes negative word go
negative=unique(unlist(strsplit(temporal_go[temporal_go$Description %like% "negative", ]$geneID, "/")))

only_positive = positive[!positive %in% negative]

only_negative = negative[!negative %in% positive]


temporal_go <- read.csv(file = 'only_C_lvl_2_all_genes_GO.csv')

#early
early = temporal_go[1:18,]
early=unique(unlist(strsplit(early$geneID, "/")))

#old
late = temporal_go[19:25,]
late=unique(unlist(strsplit(late$geneID, "/")))

only_early = early[!early %in% late]

only_late = late[!late %in% early]


#######
#MONOCLE3
#######
library("monocle3")


cds_short_names_df = data.frame(gene_short_name = rownames(only_control_cells_subset))
rownames(cds_short_names_df) = rownames(only_control_cells_subset)

cds <- new_cell_data_set(only_control_cells_subset@assays$RNA@data,
                         cell_metadata = only_control_cells_subset@meta.data,
                         gene_metadata = cds_short_names_df)
cds@int_colData$reducedDims@listData[["UMAP"]] = only_control_cells_subset@reductions[["umap"]]@cell.embeddings

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "vector_days")
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)



# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="4"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)



gene_fits <- fit_models(cds, model_formula_str = "~vector_days")
fit_coefs <- coefficient_table(gene_fits)

emb_time_terms <- fit_coefs %>% filter(term == "vector_days")

emb_time_terms %>% filter (q_value < 0.05) %>%
         select(gene_short_name, term, q_value, estimate)

plot_genes_violin(cds, group_cells_by="seurat_clusters", ncol=2) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
































































































