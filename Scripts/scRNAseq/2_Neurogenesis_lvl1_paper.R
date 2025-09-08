#####################################################################################################################
#
#####################################################################################################################


#Only all cells except ugly ones
#Idents(combined_seurat) <- "integrated_snn_res.0.6"
#only_control_cells_subset = subset(x = combined_seurat, idents=c("0","1","2","3","4","6","7","8","13","15","12","5"))


# STEP 1 - NEURONS DEVELOPMENT SUBSET CONTROL CELLS - CONTRASTS
#Generate new metadata grouping for the overall 5 groups of trajectory
only_control_cells_subset@meta.data$group_A = only_control_cells_subset@meta.data$integrated_snn_res.0.6
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "4"] <- "4"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "0"] <- "0"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "7"] <- "7"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "3"] <- "3"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "13"] <- "13"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "15"] <- "15"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "1"] <- "1"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "8"] <- "8"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "2"] <- "2"
only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "6"] <- "6"
#only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "12"] <- "12"
#only_control_cells_subset@meta.data$group_A[only_control_cells_subset@meta.data$group_A == "5"] <- "5"

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

list_of_contrasts = list(markers_BA,markers_CB,markers_DC,markers_ED)

names(list_of_contrasts) = c("markers_BA","markers_CB","markers_DC","markers_ED")

######################################################################################################################################################
#CLUSTERPROFILER
#GSEA for velocity_levels information ranking base on control vs BA contrast
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/lvl_1_neurons/GO_control/"

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
  #Fix posible problems with too much 0
  list_of_contrasts[[n]] <- list_of_contrasts[[n]][!(list_of_contrasts[[n]]$pct.1==0 & list_of_contrasts[[n]]$pct.2==0 & list_of_contrasts[[n]]$avg_log2FC==0) ,]
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
            minGSSize    = 10,
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
    #
    gseaplot_n=gseaplot2(ego, geneSetID = 1:3)
    pdf(file=paste0(OutPath,"gsea_gseaPlot_","BP","_contrast_",names(list_of_contrasts[n]),".pdf"), width=12, height=12)  
    print(gseaplot_n)
    dev.off()    
  }
  print(n)
}

######################################################################################################################################################
#TREEPLOT
#
geneList_list = list()
for (n in 1:length(list_of_contrasts)){
  print(n)
  #Fix posible problems with too much 0
  list_of_contrasts[[n]] <- list_of_contrasts[[n]][!(list_of_contrasts[[n]]$pct.1==0 & list_of_contrasts[[n]]$pct.2==0 & list_of_contrasts[[n]]$avg_log2FC==0) ,]
  #Create the geneList
  geneList <- list_of_contrasts[[n]]$avg_log2FC
  names(geneList) <- as.character(rownames(list_of_contrasts[[n]]))
  geneList = sort(geneList, decreasing = TRUE)
   
  annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', names(geneList))
  entrezgene = annotation$entrezgene_id[match(names(geneList), annotation$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]
  geneList_list = append(geneList_list, list(geneList))
  print(length(geneList_list))
}
names(geneList_list) = c("BA","CB","DC","ED")

geneList_list_tree = geneList_list
geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseGO", OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff  = 0.01)
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseKEGG", pvalueCutoff  = 0.01)
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseDO", pvalueCutoff  = 0.01)

geneList_list_tree <- pairwise_termsim(geneList_list_tree)
p1 <- treeplot(geneList_list_tree, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.0),
                    color = "NES", geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(2.5))

pdf(file= paste0("gsea_treeplot_C_kegg_BIEN.pdf"), width = 20, height =30)
print(p1)
dev.off()

######################################################################################################################################################
#####DO limma (Ta-Tb) - (Ca-Cb)
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


g <- factor(combined_seurat@meta.data$contrasts_groups)
b <- factor(combined_seurat@meta.data$vector_days)
design <- model.matrix(~0 + g + b)
design <- design[,setdiff(colnames(design), "b4")] # get to full rank

exp_matrix = combined_seurat@assays$RNA@data
#y <- voom(exp_matrix, design, plot = F)
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

######################################################################################################################################################
#PREPARING gsea input
list_of_contrasts_step2 = list(diff_expression_results_BA,diff_expression_results_CB,diff_expression_results_DC,diff_expression_results_ED)
names(list_of_contrasts_step2) = c("diff_expression_results_BA","diff_expression_results_CB","diff_expression_results_DC","diff_expression_results_ED")

OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/lvl_1_neurons/GO_T/"

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
            minGSSize    = 10,
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
    gseaplot_n=gseaplot2(ego, geneSetID = 1:3)
    pdf(file=paste0(OutPath,"gsea_gseaPlot_","BP","_contrast_",names(list_of_contrasts_step2[n]),".pdf"), width=12, height=12)  
    print(gseaplot_n)
    dev.off()      
  }
}

######################################################################################################################################################
#TREEPLOT
#
geneList_list = list()
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
  geneList_list = append(geneList_list, list(geneList))
}
names(geneList_list) = c("BA","CB","DC","ED")

geneList_list_tree = geneList_list
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseGO", OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff  = 0.01)
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseKEGG", pvalueCutoff  = 0.01)
geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseDO", pvalueCutoff  = 0.05)

geneList_list_tree <- pairwise_termsim(geneList_list_tree)
p1 <- treeplot(geneList_list_tree, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.5),
                    color = "NES", geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(4.0))

pdf(file= paste0("gsea_treeplot_T_do.pdf"), width = 20, height =30)
print(p1)
dev.off()




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



###########################################
#Auto generate all tables for each T contrast vs all C 
aa = list(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED)
aa_2 = list(potential_genes_BA,potential_genes_CB,potential_genes_DC,potential_genes_ED)

bb = list(names_sig_C_markers_BA,names_sig_C_markers_CB,names_sig_C_markers_DC,names_sig_C_markers_ED)
bb_2 = list(potential_markers_BA,potential_markers_CB,potential_markers_DC,potential_markers_ED)

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
  lista_gene_names = aa[[n]][!(aa[[n]] %in% c(names_sig_C_markers_BA,names_sig_C_markers_CB,names_sig_C_markers_DC,names_sig_C_markers_ED))]
  only_T = append(only_T, lista_gene_names)
}

only_T_conditions_all = unique(only_T)

#Only C
only_C_conditions_all_BA = names_sig_C_markers_BA[!(names_sig_C_markers_BA %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_conditions_all_CB = names_sig_C_markers_CB[!(names_sig_C_markers_CB %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_conditions_all_DC = names_sig_C_markers_DC[!(names_sig_C_markers_DC %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_conditions_all_ED = only_C_conditions_all_ED[!(only_C_conditions_all_ED %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC,names_sig_T_potential_ED))]
only_C_condition_all = unique(c(only_C_conditions_all_BA,only_C_conditions_all_CB,only_C_conditions_all_DC,only_C_conditions_all_ED))


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

#########
#TF
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

#How many TF do we have in out 1728? 
sum(rownames(out_table_plot) %in% sub(".*_", "", names(pfm))) #67 USE JASPAR 2022
tf_in_our_contrasts= rownames(out_table_plot)[rownames(out_table_plot) %in% sub(".*_", "", names(pfm))]

write.table(out_table[tf_in_our_contrasts,], file="TF_binary_of_interest_per_contrast.csv", sep=";", row.names = TRUE)
write.table(out_table_plot[tf_in_our_contrasts,], file="TF_of_interest_per_contrast.csv", sep=";", row.names = TRUE)

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
my_data[my_data < 0] <- -1
my_data[my_data > 0] <- 1
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "significant_lvl_1_genelist_both_and_T_only.pdf",
    width = 10,
    height = 10)
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
my_data[my_data < 0] <- -1
my_data[my_data > 0] <- 1
my_data_max=max(my_data)
my_data_min=min(my_data)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

pdf(file = "significant_lvl_1_genelist_C_only_LAST.pdf",
    width = 10,
    height = 10)
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

##################################
##################################
#VERSION FOR ALL VALUES PROVIDED regulat heatmaps scaled by rows
out_table_plot = out_table
out_table_plot$Events = NULL

out_table_plot[,1] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_A]))
out_table_plot[,2] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_B]))
out_table_plot[,3] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_C]))
out_table_plot[,4] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_D]))
out_table_plot[,5] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_E]))

out_table_plot[,6] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_A]))
out_table_plot[,7] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_B]))
out_table_plot[,8] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_C]))
out_table_plot[,9] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_D]))
out_table_plot[,10] = rowMeans(as.matrix(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_E]))

out_table_plot = t(scale(t(out_table_plot)))
##################################################################################################################################################
##################################################################################################################################################
#HEATMAPS SUMMARY 1-c+t+only_T AND 2-only_C
genelist_both_and_T_only = c(both_conditions_all,only_T_conditions_all)

out_table_plot_genelist_both_and_T_only= out_table_plot[genelist_both_and_T_only,]
# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 5)
tiempos=append(tiempos,rep("Treatment", 5))

#replace colnames by nice ones
table_plot_names =c("C_A","C_B","C_C","C_D","C_E","T_A","T_B","T_C","T_D","T_E")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("A","B","C","D","E","A","B","C","D","E")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("A" = "#d01c8b",
                                           "B" = "#f1b6da",
                                           "C" = "#f7f7f7",
                                           "D" = "#b8e186",
                                           "E" = "#4dac26")
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

#heatmap plot
pdf(file = "Expression_cascade_of_activation_transition_lvl_1_genelist_both_and_T_only_extra2columns.pdf",
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
tiempos=rep("Control", 5)
tiempos=append(tiempos,rep("Treatment", 5))

#replace colnames by nice ones
table_plot_names =c("C_A","C_B","C_C","C_D","C_E","T_A","T_B","T_C","T_D","T_E")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("A","B","C","D","E","A","B","C","D","E")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("A" = "#d01c8b",
                                           "B" = "#f1b6da",
                                           "C" = "#f7f7f7",
                                           "D" = "#b8e186",
                                           "E" = "#4dac26")
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

pdf(file = "Expression_cascade_of_activation_transition_lvl_1_genelist_C_only_LAST.pdf",
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

######################################################################################################################################################
######################################################################################################################################################
#TREEPLOT
# OUTCOMES GENE LIST FOR BOTH AND T together
#Create the geneList

#Create the GO
geneList <- genelist_both_and_T_only                      
ego <- enrichGO(gene       = geneList,
           universe      = rownames(combined_seurat),
           OrgDb         = org.Hs.eg.db,
           keyType       = 'SYMBOL',
           ont           = 'BP',
           pAdjustMethod = "BH",
           minGSSize     = 10,
           maxGSSize     = 500,
           pvalueCutoff  = 0.05)
ego_table <- as.data.frame(ego) 
write.csv(ego_table,"gsea_treeplot_outcomes_go_bp.csv")


ego2 <- pairwise_termsim(ego)
p1 <- treeplot(ego2, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.0),
                    geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(5.0))
                    
pdf(file= paste0("gsea_treeplot_outcomes_go_bp.pdf"), width = 20, height =30)
print(p1)
dev.off()

#Create the KEGG
geneList <- genelist_both_and_T_only
annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', geneList)
entrezgene = annotation$entrezgene_id[match(geneList, annotation$hgnc_symbol)]
geneList = as.character(entrezgene)
geneList = geneList[!is.na(geneList)]

geneList <- geneList
ego <- enrichKEGG(gene      = geneList,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
ego_table <- as.data.frame(ego)                  
write.csv(ego_table,"gsea_treeplot_outcomes_kegg.csv")


ego2 <- pairwise_termsim(ego)
p1 <- treeplot(ego2, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.0),
                    geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(5.0))
                    
pdf(file= paste0("gsea_treeplot_outcomes_kegg.pdf"), width = 20, height =30)
print(p1)
dev.off()

#Create the DO
library(DOSE)
geneList <- genelist_both_and_T_only
annotation = name_changer('ensembl', 'hsapiens_gene_ensembl', geneList)
entrezgene = annotation$entrezgene_id[match(geneList, annotation$hgnc_symbol)]
geneList = as.character(entrezgene)
geneList = geneList[!is.na(geneList)]

geneList <- geneList
ego <- enrichDO(gene          = geneList,
              ont           = "DO",
              pAdjustMethod = "BH",
              minGSSize     = 10,
              maxGSSize     = 500,
              pvalueCutoff  = 0.05,
              readable      = FALSE)
ego_table <- as.data.frame(ego) 
write.csv(ego_table,"gsea_treeplot_outcomes_do.csv")


ego2 <- pairwise_termsim(ego)
p1 <- treeplot(ego2, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.0),
                    geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(5.0))
                    
pdf(file= paste0("gsea_treeplot_outcomes_do.pdf"), width = 20, height =30)
print(p1)
dev.off()


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

#
identities=c("0","1","2","3","4","5","6","7","8","12","13","15")
my_color_palette <- hue_pal()(length(identities))
cds$group_A=factor(cds$group_A)
#
plot_cells(cds,
           color_cells_by = "group_A",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5) + scale_color_manual(values = my_color_palette)



# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="4"){
  cell_ids <- which(colData(cds)[, "group_A"] == time_bin)
  
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

#######
#CHECK GENES TO TRY YO FIT TRAYECTORY ON NEUROGENESIS LVL_1
#######
Idents(only_control_cells_subset) = only_control_cells_subset@meta.data$group_A
Idents(only_control_cells_subset) <- factor(Idents(only_control_cells_subset), levels= c("A","B","C","D","E"))

#level 1: Hes5/PAX6/DCX/CALB1
#level 2: NES/Sox11/MAP2/DCX/STMN1/SYP/CALB1

pdf(file="Characterization_neurogenesis.pdf", width=10, height=6)
features = unique(c("HES5","PAX6","DCX","SYP","CALB1"))
DotPlot(only_control_cells_subset, features = features, cols = c("blue", "yellow")) + RotatedAxis()  #only_control_cells_subset_neurons combined_seurat
dev.off()

#Plot Genes from Maite OPTION _1 AND _2
pdf(file="Characterization_neurogenesis_2_1_umap.pdf", width=12, height=12)
p1 = DimPlot(only_control_cells_subset, reduction = 'umap', group.by = "group_A",label = TRUE)

CombinePlots(list(p1))
dev.off()

pdf(file="Characterization_neurogenesis_2_2.pdf", width=12, height=12)
p2 = FeaturePlot(
      object = only_control_cells_subset,
      features = c("HES5"),
      ncol = 1
      )
CombinePlots(list(p2))
dev.off()

pdf(file="Characterization_neurogenesis_2_3.pdf", width=12, height=12)
p3 = FeaturePlot(
      object = only_control_cells_subset,
      features = c("PAX6"),
      ncol = 1
      )
CombinePlots(list(p3))
dev.off()

pdf(file="Characterization_neurogenesis_2_4.pdf", width=12, height=12)
p4 = FeaturePlot(
      object = only_control_cells_subset,
      features = c("DCX"),
      ncol = 1
      )
CombinePlots(list(p4))
dev.off()

pdf(file="Characterization_neurogenesis_2_5.pdf", width=12, height=12)
p5 = FeaturePlot(
      object = only_control_cells_subset,
      features = c("SYP"),
      ncol = 1
)
CombinePlots(list(p5))
dev.off()

pdf(file="Characterization_neurogenesis_2_6.pdf", width=12, height=12)
p5 = FeaturePlot(
      object = only_control_cells_subset,
      features = c("CALB1"),
      ncol = 1
)
CombinePlots(list(p5))
dev.off()

#######
#MOUNT SIGNIFICANT THINGS TABLES 
#######
lista_all[[1]]
lista_all[[2]]
lista_all[[3]]
lista_all[[4]]

#GENES ONLY IN CONTROL
lista_gene_names = names_sig_C_markers_DC[!(names_sig_C_markers_DC %in% c(names_sig_T_potential_BA,names_sig_T_potential_CB,names_sig_T_potential_DC))]
potential_markers_DC[lista_gene_names,]

#GENES ONLY IN TREATMENT
lista_gene_names = names_sig_T_potential_DC[!(names_sig_T_potential_DC %in% c(names_sig_C_markers_BA,names_sig_C_markers_CB,names_sig_C_markers_DC))]
potential_genes_DC[lista_gene_names,]

#######
#PSEUDOTIME HEATMAP
#######
library(RColorBrewer)

genes = c(both_conditions_all,only_T_conditions_all)

pt.matrix = as.matrix(cds@assays@data$counts[match(genes,rowData(cds)[,1]),order(pseudotime(cds))])
pt.matrix = t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix = t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) = genes

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
my_data <- out_table_plot_genelist_both_and_T_only
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

ht = Heatmap(pt.matrix,
             name                         = "z-score",
             col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
             show_row_names               = FALSE,
             show_column_names            = FALSE,
             row_names_gp                 = gpar(fontsize = 6),
             km = 6,
             row_title_rot                = 0,
             cluster_rows                 = TRUE,
             cluster_row_slices           = FALSE,
             cluster_columns              = FALSE,
             right_annotation             = rowAnno)

pdf(file = "Pseudotime_both_C_T_heatmap.pdf",width = 10,height = 10)
print(ht)
dev.off() 

#######
#PATHVIEW WITH MULTIPLE EXPRESSION VSALUES IN SAME PLOT
#######
library("pathview")

data.frame(ego)$ID
#"hsa03010" "hsa05012" "hsa04723" "hsa00190" "hsa04714" "hsa05010" "hsa04360" "hsa04550"

#Create the KEGG
genes_present_in_entrez_from_both_and_only_T= rownames(my_data)[rownames(my_data) %in% annotation$hgnc_symbol]
entrezgene = annotation$entrezgene_id[match(genes_present_in_entrez_from_both_and_only_T, annotation$hgnc_symbol)]

my_data_corrected = my_data[genes_present_in_entrez_from_both_and_only_T,]
rownames(my_data_corrected) = entrezgene

colnames(my_data_corrected) = colnames(my_data)
my_data_corrected$GeneID = rownames(my_data_corrected)

setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/Neurogenesis_lvl_1/pathways/")

for (n in 1:length(data.frame(ego)$ID)){
  pathview(gene.data  = as.matrix(my_data_corrected[1:8]),
                       pathway.id = data.frame(ego)$ID[n],
                       species    = "hsa",
                       limit      = c(-0.10, 0.10))
}


#######
#PSEUDOTIME HEATMAP ONLY_CONTROL
#######
library(RColorBrewer)
library("ComplexHeatmap")
library("circlize")

genes = c(only_C_condition_all)
genes_temp = c("TPM1","HES4","METRN","RPL22L1","GNG12","LAMA4","FRAS1","LIX1","IFITM3","IDH1","FAM181A","CYR61","CXCR4","FABP7","S100A13","MYO10","SMOC1","TFPI2","ZFP36L1","PTPRZ1","SPECC1","PON2","SOX3","IGFBP2","PAX6","VIM","ID4","GPM6B","NTRK2","AL139246.5","FGFBP3","PLP1","ENPP2","SFRP2","TTYH1",
"SLC1A3","FBXO32","ITGB8","SOX9","QKI","ATP6V0E1","ANXA5","EZR","LITAF","RGMB","SCRT2","ARL4D","ST18","INSM1","C1QL1","EPB41","DCC","ELAVL2","MAP6","DCX","TCF4","TFDP2","JAG1","BTG2","PCDH18","NEUROG2","MAGI1","SYT2","GADD45G","LRRN1","PHLDA1","ENC1","TAGLN3","PCBP4","MIAT","CBFA2T2","KCNQ1OT1",
"CDKN1C","PLXNA2","SINHCAF","TCF12","NHLH1","NEUROD4","NKAIN4","MEIS2","SCG3","GDAP1L1","GNG3","INA","STMN4","ZC2HC1A","CADM1","RUNX1T1","PGM2L1","CSRNP3","TMX4","NTM","PCDH9","GRIA4","PCLO","ZFHX4","SCN2A","MAPT","REEP1","CELF4","SH3BP5","SYT4","PRKAR2B","SPTBN1","NRXN1","TTC9B","ANK3","ANK2",
"GAP43","PCSK1N","TUBB2A","PPP2R2B")

genes_temp_lvl_2 = c("SPX","CBLN2","CACNA2D1","NTM","SNCA","FXYD6","CENPV","TFPI2","RFX4","PGAP1","SLC17A8","CRYBA2","FEV","GATA3-AS1","GATA3","ID4","NOVA1","ANK3","NRXN1","LSAMP","REEP1","TMX4","MAPT","GAP43","UCHL1","SPTBN1","SPOCK1","FGF13","LUZP2","SCN2A","TMEFF2","CDH7","NAV3","LIMCH1","RPRM","GRIA4","RALYL","TMSB4X","ENC1","DOK5","TUBB2A","PHLDA1","TMSB15A","OLFM1","STMN2","SLC17A6","MAB21L2","ZNF503","HOXB3","TUBB2B","TUBB","COTL1","ACTG1","MTRNR2L8","MTRNR2L10","SLC25A6","CLIC1","PDLIM7","CKB","PCBP4","NKAIN4","CDKN1C","NEAT1","MDK","CRABP1","IGDCC3","CNTN2","ARL4D","KLHL35","NHLH1","CRYBA1","GNG5","SINHCAF","ST18","RPS11","TCF12","RPS20","VCAN","ERBB4","RUNX1T1","TENM2","ZNF536","NR2F2","FOXP2","FIGN","LINC00461","NRP2","RPS27L","TAGLN3","CCND1","IRX5","BTG1","SOX11","SPATS2L","MEIS2","LHX1","NR2F1","RBMS1","LHX5-AS1","PCSK2","VIM")


pt.matrix = as.matrix(cds@assays@data$counts[match(genes,rowData(cds)[,1]),order(pseudotime(cds))])
pt.matrix = t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix = t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) = genes

#ADD colors and sizes to rownames base on groups 
# Rows to highlight
myRows <- only_C_condition_all
myRows_t <- only_C_condition_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(my_data) %in% myRows) #BOTH
row_idx_t <- which(rownames(my_data) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(my_data))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(my_data))
fontcolors[row_idx] <- 'd8b365'
fontcolors[row_idx_t] <- '#d8b365'

fontfaces <- rep('plain',nrow(my_data))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
my_data <- out_table_plot_genelist_C_only
rowAnno <- rowAnnotation(rows = anno_text(rownames(my_data), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

pt.matrix = pt.matrix[rownames(my_data),]

ht = Heatmap(pt.matrix,
             name                         = "z-score",
             col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
             show_row_names               = FALSE,
             show_column_names            = FALSE,
             row_names_gp                 = gpar(fontsize = 6),
             km = 6,
             row_title_rot                = 0,
             cluster_rows                 = TRUE,
             cluster_row_slices           = FALSE,
             cluster_columns              = FALSE,
             right_annotation             = rowAnno)

pdf(file = "Pseudotime_C_heatmap.pdf",width = 10,height = 10)
print(ht)
dev.off() 


#ALL green genes from lvl1 and lvl2
genelist_both_and_T_only = c("NOL7","CLUAP1","SEC61G","SRP14","SLIRP","SNRPF","TMEM59","ROMO1","UBE2E3","STRAP","ENAH",
"FRMD4A","RAMMET","YPEL5","RUNDC3A","KALRN","EIF3I","STRAP","ITGAE","RPS8","H3F3A","FAM216A","OAZ1","TKT",
"RPSA","CASP8AP2","PCYOX1","SFRP2","KRAS","ST13","IFT81","PRDX3","PFDN4","TPRKB","C11orf1","RPL21","HIP1R","ANKLE2","UBN1","ZNF618","LMO3",
"SMC6","ACVR2B","NUP54","ITFG1","TMEM169","INSR","DLG5")














#DO DE between T-C in cluster E 

Idents(combined_seurat) = "protocol_2"

de_genes_between_T_C_in_cluster_E= FindMarkers(object = combined_seurat, ident.1 = "basal", ident.2 = "control", only.pos = FALSE)
de_genes_between_T_C_in_cluster_E = de_genes_between_T_C_in_cluster_E[abs(de_genes_between_T_C_in_cluster_E$p_val_adj) < 0.05, ]



write.csv(top.table_dani,"top.table_dani.csv")

top.table_dani = read.csv("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/validation_tfg_alejandro_data/top.table_dani.csv")

rownames(de_genes_between_T_C_in_cluster_E)
top.table_dani$hgnc_symbol


overlapping_genes_significant_in_our_data = top.table_dani$hgnc_symbol[top.table_dani$hgnc_symbol %in% rownames(de_genes_between_T_C_in_cluster_E)]

temp = top.table_dani[top.table_dani$hgnc_symbol %in% overlapping_genes_significant_in_our_data,]
temp = temp[temp$AD_vs_Control.sig_genes !=0,]

rownames(temp)=unlist(temp['hgnc_symbol'])
temp[, c(2,22,26,28)]


de_genes_between_T_C_in_cluster_E[rownames(temp),]$avg_log2FC

temp_final=temp[, c(2,22,26,28)]

#PUTTING TOGETHER OVERLAPPING SIGNIFICANT GENES FROM BOTH ANALYSIS 
temp_final['mio_log2FC'] = de_genes_between_T_C_in_cluster_E[rownames(temp),]$avg_log2FC
temp_final['mio_p_val_adj'] = de_genes_between_T_C_in_cluster_E[rownames(temp),]$p_val_adj

#TAKE A LOOK OF SIGNS AND CHECK NUMBER OF COINCIDENDES :) ANY GENES ARE RELEVANT FOR US IN OUR ANALYSIS????
temp_final$same_direction = ifelse((temp_final$AD_vs_Control.logFC > 0 & temp_final$mio_log2FC > 0), 1, ifelse((temp_final$AD_vs_Control.logFC < 0 & temp_final$mio_log2FC < 0), -1, 0))
table(temp_final$same_direction)

relevant_genes = c('GUK1','ZNF428','ARL6IP5','APC','ENAH','DCTN3','SRP14','FAM57B','POU2F2','TUBB4A','STMN2','THSD7A','NCAM1','NSG1','NEFL','RBP1','NEFM','NUDT14','SNCG','ONECUT2','SCG2','UCHL1','TUBB2B','DDAH2','TUBB','EBF1','KIF5C','TUBA1A','PTMA','ROMO1','PTN','SEC61G','NDUFB1','SEC62','SLIRP','TMEM59','IER2','COTL1','PRPF40A','STRAP','IDH2','UBE2E3','LMO1','RND3','SCRG1','SOX11','HOXB2','DRAXIN','KLHL35','GLRX','CRABP1','CBLB','CRYBA1','SH3BGRL3','PPP1R14C','DLL1','CDC25B','DLL3','SPSB4','PRMT8','HES6','NIN','CCND1','ASCL1','SNRPF','C1orf61','RPL37','S100B','BCAN','TMEM161B-AS1','CLU','COPZ1','CD99','TIMP1','SPARC','B2M','MYL12A','SFRP1','OPRK1','FOS','PPM1K','MIR99AHG','IGFBP5','SLC2A1','TCF7L1','HES5','NES','RGMA','EIF4H','CDC34','POLR1D','GLUL','PRSS23','ZBTB16','DIO3','IGDCC3','SIPA1L2','NCALD','GNG5','RPS27L','MAGED2','PCDH17','MYL6','SKA2','ANXA2','ID1','ID3','CTGF','LIMCH1','CLUAP1','NOL7')

relevant_genes[relevant_genes %in% rownames(temp_final)]
# "STMN2" "NEFM"  "ASCL1"

write.csv(temp_final,"overlapping genes_between AD_C_in_both_our_Data_and_online_data.csv")





#Let's calculate also de rank correlation for overlapping total genes in the contrasts :)
de_genes_between_T_C_in_cluster_E = FindMarkers(object = combined_seurat, ident.1 = "ba", ident.2 = "control", only.pos = FALSE, min.pct = 0, logfc.threshold = 0, features = rownames(combined_seurat), return.thresh = 1.01)
top.table_dani

de_genes_between_T_C_in_cluster_E = de_genes_between_T_C_in_cluster_E[order(de_genes_between_T_C_in_cluster_E$avg_log2FC), ]
top.table_dani = top.table_dani[order(top.table_dani$AD_vs_Control.logFC), ]

overlapping_genes_significant_in_our_data = top.table_dani$hgnc_symbol[top.table_dani$hgnc_symbol %in% rownames(de_genes_between_T_C_in_cluster_E)]

de_genes_between_T_C_in_cluster_E = de_genes_between_T_C_in_cluster_E[overlapping_genes_significant_in_our_data,]
top.table_dani = top.table_dani[top.table_dani$hgnc_symbol %in% overlapping_genes_significant_in_our_data,]

res <- cor.test(de_genes_between_T_C_in_cluster_E$avg_log2FC, top.table_dani$AD_vs_Control.logFC, method = "pearson")






######################################################################
######################################################################

#New Heatmap with 2 extra columns mean of ba and mean of control
setwd("/ibex/scratch/projects/c2169/Xabier/Navarra/Neuro")
load("/ibex/scratch/projects/c2169/Xabier/Navarra/Neuro/21_09_2022_neurogenesis_lvl1.RData")

top.table_dani = read.csv("/ibex/project/c2169/Xabier/Navarra/Neuro/top.table_dani.csv")
temp_final = read.csv("/ibex/project/c2169/Xabier/Navarra/Neuro/overlapping genes_between AD_C_in_both_our_Data_and_online_data.csv")


library("Seurat")
library("ComplexHeatmap")
library("circlize")

#Plot HEATMAP with paper data for our list of interesting lvl_1 genes and see if ad and control groups are separated thx to them
rownames(top.table_dani) = make.unique(top.table_dani$hgnc_symbol, sep = ".")
paper_data_matrix_1 = as.matrix(top.table_dani[,3:20])

#Now for mature version 2
load("/ibex/scratch/projects/c2169/Xabier/Navarra/Neuro/28_09_2022_neurogenesis_lvl2.RData")
out_table_plot_genelist_both_and_T_only
both_conditions_all
only_T_conditions_all

#all genes for heatmap
#relevant_genes = c('GUK1','ZNF428','ARL6IP5','APC','ENAH','DCTN3','SRP14','FAM57B','POU2F2','TUBB4A','STMN2','THSD7A','NCAM1','NSG1','NEFL','RBP1','NEFM','NUDT14','SNCG','ONECUT2','SCG2','UCHL1','TUBB2B','DDAH2','TUBB','EBF1','KIF5C','TUBA1A','PTMA','ROMO1','PTN','SEC61G','NDUFB1','SEC62','SLIRP','TMEM59','IER2','COTL1','PRPF40A','STRAP','IDH2','UBE2E3','LMO1','RND3','SCRG1','SOX11','HOXB2','DRAXIN','KLHL35','GLRX','CRABP1','CBLB','CRYBA1','SH3BGRL3','PPP1R14C','DLL1','CDC25B','DLL3','SPSB4','PRMT8','HES6','NIN','CCND1','ASCL1','SNRPF','C1orf61','RPL37','S100B','BCAN','TMEM161B-AS1','CLU','COPZ1','CD99','TIMP1','SPARC','B2M','MYL12A','SFRP1','OPRK1','FOS','PPM1K','MIR99AHG','IGFBP5','SLC2A1','TCF7L1','HES5','NES','RGMA','EIF4H','CDC34','POLR1D','GLUL','PRSS23','ZBTB16','DIO3','IGDCC3','SIPA1L2','NCALD','GNG5','RPS27L','MAGED2','PCDH17','MYL6','SKA2','ANXA2','ID1','ID3','CTGF','LIMCH1','CLUAP1','NOL7')
relevant_genes = rownames(out_table_plot_genelist_both_and_T_only)

paper_data_matrix_1 = paper_data_matrix_1[rownames(paper_data_matrix_1) %in% relevant_genes,]
paper_data_matrix_1_scales = t(scale(t(paper_data_matrix_1)))

#both_conditions_all = c('GUK1','ZNF428','ARL6IP5','APC','DCTN3','FAM57B','POU2F2','TUBB4A','STMN2','THSD7A','NCAM1','NSG1','NEFL','RBP1','NEFM','NUDT14','SNCG','ONECUT2','SCG2','UCHL1','TUBB2B','DDAH2','TUBB','EBF1','KIF5C','TUBA1A','PTMA','PTN','NDUFB1','SEC62','IER2','COTL1','PRPF40A','IDH2','LMO1','RND3','SCRG1','SOX11','HOXB2','DRAXIN','KLHL35','GLRX','CRABP1','CBLB','CRYBA1','SH3BGRL3','PPP1R14C','DLL1','CDC25B','DLL3','SPSB4','PRMT8','HES6','NIN','CCND1','ASCL1','C1orf61','RPL37','S100B','BCAN','TMEM161B-AS1','CLU','COPZ1','CD99','TIMP1','SPARC','B2M','MYL12A','SFRP1','OPRK1','FOS','PPM1K','MIR99AHG','IGFBP5','SLC2A1','TCF7L1','HES5','NES','RGMA','EIF4H','CDC34','POLR1D','GLUL','PRSS23','ZBTB16','DIO3','IGDCC3','SIPA1L2','NCALD','GNG5','RPS27L','MAGED2','PCDH17','MYL6','SKA2','ANXA2','ID1','ID3','CTGF','LIMCH1')
both_conditions_all = both_conditions_all[both_conditions_all %in% rownames(paper_data_matrix_1)]

#only_T_conditions_all = c("ENAH", "SRP14", "ROMO1", "SEC61G", "SLIRP", "TMEM59", "STRAP", "UBE2E3", "SNRPF", "CLUAP1", "NOL7")
only_T_conditions_all = only_T_conditions_all[only_T_conditions_all %in% rownames(paper_data_matrix_1)]

#GENERATE HEATMAP PLOT
# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("AD", 8)
tiempos=append(tiempos,rep("Control", 10))

#replace colnames by nice ones
table_plot_names = colnames(paper_data_matrix_1)

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "AD" = "#5ab4ac"))

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

#ADD colors and sizes to rownames base on groups
# Rows to highlight
myRows <- both_conditions_all
myRows_t <- only_T_conditions_all
# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(paper_data_matrix_1) %in% myRows) #BOTH
row_idx_t <- which(rownames(paper_data_matrix_1) %in% myRows_t)  #ONLY in Treatment

fontsizes <- rep(5, nrow(paper_data_matrix_1))
fontsizes[row_idx] <- 5
fontsizes[row_idx_t] <- 5

fontcolors <- rep('#d8b365', nrow(paper_data_matrix_1))
fontcolors[row_idx] <- 'black'
fontcolors[row_idx_t] <- '#5ab4ac'

fontfaces <- rep('plain',nrow(paper_data_matrix_1))
fontfaces[row_idx] <- 'bold'
fontfaces[row_idx_t] <- 'bold'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(paper_data_matrix_1), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

#Gradient for EXPESSION data
my_data_max=max(paper_data_matrix_1_scales)
my_data_min=min(paper_data_matrix_1_scales)
col_fun = colorRamp2(c(my_data_min, 0, my_data_max), c("blue", "white", "red"))

#Compute average of ba and average of control and add them to the matrix
ad_average = rowMeans(paper_data_matrix_1_scales[,1:8])
control_average = rowMeans(paper_data_matrix_1_scales[,9:18])

#Add the additional 2 columns information to generate the new plot
paper_data_matrix_1_scales = data.frame(paper_data_matrix_1_scales)

paper_data_matrix_1_scales["AD_Average"] = ad_average
paper_data_matrix_1_scales["Control_Average"] = control_average

library(dplyr)
paper_data_matrix_1_scales <- paper_data_matrix_1_scales %>%
  arrange(AD_Average - Control_Average)

#Fix annotation data to acomodate the 2 new columns
new_rows <- data.frame(
  Condition = c("AD_Average","Control_Average")
)
meta_data_plot <- rbind(meta_data_plot, new_rows)

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                          "AD" = "#5ab4ac",
                                          "Control_Average" = "#b3844e",
                                          "AD_Average" = "#47928a"
                                          )
                          )

ha = HeatmapAnnotation(df = meta_data_plot,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

#New plot version with the two columns
pdf(file = "Online_dataset_1_Leonardo_with_lvl_1_interesting_genes_2222_MATURE.pdf",
    width = 10,
    height = 10)
Heatmap(
  paper_data_matrix_1_scales,
  heatmap_legend_param = list(title = "Neurogenesis in Control and Treatment"),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = colnames(paper_data_matrix_1_scales),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  bottom_annotation = NULL,
  top_annotation = ha,
  right_annotation = rowAnno
)
dev.off()

#Plot PCA on online dataset using only the genes interesting from our data
#lets load the normalized data to subset genes and compute PCA (lets see if samples are separated for AD and Control groups in plot)

load("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/validation_tfg_alejandro_data/session.RData")
norm_counts

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#For random PCA
#random_relevant_genes_selection = sample(rownames(norm_counts), 111)   ##list_sig_genes from validation data

anno = getBM(
  values = rownames(norm_counts), #El env�o
  filters = c("ensembl_gene_id"), #Categor�a de env�o
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), #Categor�a de retorno
  mart = human)
our_interesting_genes = anno[(anno$hgnc_symbol %in% relevant_genes),]
interesting_genes_subseted_nom_counts = norm_counts[our_interesting_genes$ensembl_gene_id,]

potential_genes_1_13
potential_genes_6_1
potential_genes_8_6
potential_genes_2_13
our_interesting_genes = anno[(anno$hgnc_symbol %in% rownames(potential_genes_6_1)),]
interesting_genes_subseted_nom_counts = norm_counts[our_interesting_genes$ensembl_gene_id,]

# PCA
pca <- prcomp(t(interesting_genes_subseted_nom_counts))
pr <- summary(pca)$importance[,1:5]
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], groups=design_matrix$combine)
pdf(file = "PCA_LVL2_ALL_RELEVANT_GENES_only_6_1.pdf", width = 10, height = 10)
pca1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = groups, label=rownames(mplot))) +
  geom_point(size=4) +
  scale_color_manual(values=color.palette) +
  labs(title="PCA - Normalized counts",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))
pca1 <- pca1 + geom_text(aes(label=rownames(design_matrix)),hjust=-0.1, vjust=-0.1)
pca1 <- pca1 + theme_classic()
pca1
dev.off()

#Generate a gene_set with the different gene subsets to try and compute a GSEA against the ranked genes logFC
library("AUCell")

geneSets <- list(geneSet1=c('GUK1','ZNF428','ARL6IP5','APC','ENAH','DCTN3','SRP14','FAM57B','POU2F2','TUBB4A','STMN2','THSD7A','NCAM1','NSG1','NEFL','RBP1','NEFM','NUDT14','SNCG','ONECUT2','SCG2','UCHL1','TUBB2B','DDAH2','TUBB','EBF1','KIF5C','TUBA1A','PTMA','ROMO1','PTN','SEC61G','NDUFB1','SEC62','SLIRP','TMEM59','IER2','COTL1','PRPF40A','STRAP','IDH2','UBE2E3','LMO1','RND3','SCRG1','SOX11','HOXB2','DRAXIN','KLHL35','GLRX','CRABP1','CBLB','CRYBA1','SH3BGRL3','PPP1R14C','DLL1','CDC25B','DLL3','SPSB4','PRMT8','HES6','NIN','CCND1','ASCL1','SNRPF','C1orf61','RPL37','S100B','BCAN','TMEM161B-AS1','CLU','COPZ1','CD99','TIMP1','SPARC','B2M','MYL12A','SFRP1','OPRK1','FOS','PPM1K','MIR99AHG','IGFBP5','SLC2A1','TCF7L1','HES5','NES','RGMA','EIF4H','CDC34','POLR1D','GLUL','PRSS23','ZBTB16','DIO3','IGDCC3','SIPA1L2','NCALD','GNG5','RPS27L','MAGED2','PCDH17','MYL6','SKA2','ANXA2','ID1','ID3','CTGF','LIMCH1','CLUAP1','NOL7'))
geneSets = our_interesting_genes$ensembl_gene_id[our_interesting_genes$hgnc_symbol %in% geneSets$geneSet1]

# Calculate enrichment scores
cells_rankings <- AUCell_buildRankings(counts_matrix1, plotStats=FALSE)
cells_AUC <- AUCell_run(counts_matrix1, geneSets, aucMaxRank=nrow(cells_rankings)*0.05)


######################################################################################################################################################
#CLUSTERPROFILER
#GSEA
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/validation_tfg_alejandro_data/"

library("enrichplot")
library("clusterProfiler")
library("org.Hs.eg.db")

list_of_contrasts = list(top.table_dani)
list_of_gsea = list()

for (n in 1:length(list_of_contrasts)){
  print(n)
  #Fix posible problems with too much 0
  list_of_contrasts[[n]] <- list_of_contrasts[[n]][!(list_of_contrasts[[n]]$AD_vs_Control.logFC==0) ,]
  #Create the geneList
  geneList <- list_of_contrasts[[n]]$AD_vs_Control.logFC
  names(geneList) <- as.character(list_of_contrasts[[n]]$hgnc_symbol)
  geneList = sort(geneList, decreasing = TRUE)

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  anno = getBM(
  values = names(geneList), #El env�o
  filters = c("hgnc_symbol"), #Categor�a de env�o
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), #Categor�a de retorno
  mart = human)

  entrezgene = anno$entrezgene_id[match(names(geneList), anno$hgnc_symbol)]
  names(geneList) = as.character(entrezgene)
  geneList = geneList[!is.na(names(geneList))]
  
  #GSEA--------------------------------------------------------------------------------------------------------------------------si no se peude con esto con el de abajo
  ego <- gseGO(geneList    = geneList,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 10,
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
    #
    gseaplot_n=gseaplot2(ego, geneSetID = 1:3)
    pdf(file=paste0(OutPath,"gsea_gseaPlot_","BP","_contrast_",names(list_of_contrasts[n]),".pdf"), width=12, height=12)  
    print(gseaplot_n)
    dev.off()    
  }
  print(n)
}








#Generate the ranked genelist from human online data
list_of_contrasts = top.table_dani
list_of_contrasts = list_of_contrasts[!(list_of_contrasts['AD_vs_Control.logFC']==0) ,]
#Create the geneList
geneList = list_of_contrasts['AD_vs_Control.logFC'] 
geneList = unlist(geneList)
names(geneList) = as.character(list_of_contrasts$entrezgene_id)
geneList = sort(geneList, decreasing = TRUE)
geneList = geneList[!is.na(names(geneList))]

#Get entrezID from my list of interesting genes
geneSets <- list(geneSet1=c('GUK1','ZNF428','ARL6IP5','APC','ENAH','DCTN3','SRP14','FAM57B','POU2F2','TUBB4A','STMN2','THSD7A','NCAM1','NSG1','NEFL','RBP1','NEFM','NUDT14','SNCG','ONECUT2','SCG2','UCHL1','TUBB2B','DDAH2','TUBB','EBF1','KIF5C','TUBA1A','PTMA','ROMO1','PTN','SEC61G','NDUFB1','SEC62','SLIRP','TMEM59','IER2','COTL1','PRPF40A','STRAP','IDH2','UBE2E3','LMO1','RND3','SCRG1','SOX11','HOXB2','DRAXIN','KLHL35','GLRX','CRABP1','CBLB','CRYBA1','SH3BGRL3','PPP1R14C','DLL1','CDC25B','DLL3','SPSB4','PRMT8','HES6','NIN','CCND1','ASCL1','SNRPF','C1orf61','RPL37','S100B','BCAN','TMEM161B-AS1','CLU','COPZ1','CD99','TIMP1','SPARC','B2M','MYL12A','SFRP1','OPRK1','FOS','PPM1K','MIR99AHG','IGFBP5','SLC2A1','TCF7L1','HES5','NES','RGMA','EIF4H','CDC34','POLR1D','GLUL','PRSS23','ZBTB16','DIO3','IGDCC3','SIPA1L2','NCALD','GNG5','RPS27L','MAGED2','PCDH17','MYL6','SKA2','ANXA2','ID1','ID3','CTGF','LIMCH1','CLUAP1','NOL7'))

#FOR LVL2
geneSets <- list(geneSet1=c('RPL21','AKAP12','HIP1R','THSD7A','YPEL5','SRSF10','INSR','ZNF618','KALRN','TMEM169','KMT2A','LMO3','RPSA','ACVR2B','TKT','RPS2','RBP1','RPL10','RPL18','RPL29','RPL18A','MTRNR2L12','POU2F2','VSTM2L','TUBA1A','C1QL1','PMEL','RGMB','CBLB','GLRX','NKX6-1','EPB41','CYP26A1','SFRP1','SOX3','DCC','TPRKB','SFRP2','NUP54','SMC6','SRRM4','PTPRG','ZNF703','ACTB','SNCG','TUBA1B','MAP4','MIR99AHG','ONECUT2','FJX1','IGLON5','CNR1','CPE','SYT4','EEF1B2','LINC01828','PCDH9','NEFL','NEFM','EBF1','SLIT2','COL19A1','MEIS1','PCBP1','PCYOX1','MRPL34','SLC25A4','GNG2','DTD1','CCNI','ONECUT1','CIB2','TXN','SHOX2','SMS','RABAC1','ZNF706','EEF1A2','RUNDC3A','CXXC5','GAPDH','OAZ1','PRDX3','RAMMET','BLOC1S2','ATP6V1A','ANKLE2','PFDN4','MGAT4C','ZFHX4','RTN1','NME1','TSPAN13','PLEKHA5','MAPRE2','GNAO1','MYO5A','NAP1L3','YBX1','TMOD2','WRB','KRAS','FAM216A','YWHAH','C19orf70','RAC3','CADM2','TMEM14A','SV2A','VGF','CALY','ZMYND11','CHMP5','PPIP5K2','HSPH1','CELF4','PCLO','SCG2','PREPL','TPM1','NSG2','TMSB10','SYT1','PEG10','CHGB','GCHFR','DNAJC12','CADM1','PAM','ZSCAN18','MEST','ANKS1B','CASD1','CALB1','GATA2','CPNE4','SSTR2','ISG15','DLL3','GNAI2','SF3B5','RND3','IRX3','CRNDE','PDZRN4','SOX2','ASCL1','FOS','EGR1','UBN1','ST13','ITGAE','FRMD4A','NXPH1','EIF3I','STRAP','DLG5','ITFG1','SCRG1','MALAT1','PAX6','EPHB1','EPHA4','MID1','ANKRD10','PBX3','SNRPF','C1orf61','DUSP4','CDH2','UQCRQ','SEC62','GPM6A','SYNRG','ZFHX3','DSEL','DENR','C11orf1','IFT81','CASP8AP2','DICER1','TOP2B','C6orf62','TSHZ2','CDH8','CHD9','ARHGEF12','FAT3','ZMIZ1','MLLT3','SCD5','FAM200B','NDUFB3','CNTNAP2','SLIRP','MT-ND3','NDUFA3','ATP5F1E','NDUFB1','NDUFC1','ATP5MD','MRPL33','RPS29','RPL38','RPS8','RPL36','RPL39','RPL37','RPS27','H3F3A'))

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
anno = getBM(
values = names(geneList), #El env�o
filters = c("hgnc_symbol"), #Categor�a de env�o
attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), #Categor�a de retorno
mart = human)

entrezgene = anno$entrezgene_id[match(geneSets$geneSet1, anno$hgnc_symbol)]
my_geneSet = as.character(entrezgene)
my_geneSet = my_geneSet[!is.na(my_geneSet)]
  
#MIRAR CON QUE LIBRERIA SE PUEDE HACER CUSTOM GENE SETS PARA CORRER GSEA
library(fgsea)
examplePathways
exampleRanks
my_geneSet_custom = examplePathways[1]
names(my_geneSet_custom) = c("custom_geneSet_sig_Genes")
my_geneSet_custom$custom_geneSet_sig_Genes = my_geneSet



fgseaRes <- fgsea(pathways = my_geneSet_custom_test, 
                  stats    = geneList,
                  minSize  = 10,
                  maxSize  = 500)


#Lets try to take top 50 genes base on PC1 contribution and recalculate the fgsea
pca_matrix = data.frame(pca$rotation)
pca_matrix["PC1"]

values = data.frame(abs(pca_matrix[,"PC1"]))
rownames(values) = rownames(pca_matrix)
values = as.matrix(values)
values = data.frame(sort(values[,1], decreasing=TRUE))


#Blucle searcher
list_of_partial_pval = list()

for (n in 1:dim(values)[1]){
  genes_topN = rownames(values)[1:n]
  
  entrezgene = anno$entrezgene_id[match(genes_topN, anno$ensembl_gene_id)]
  my_geneSet = as.character(entrezgene)
  my_geneSet = my_geneSet[!is.na(my_geneSet)]
   
  my_geneSet_custom = examplePathways[1]
  names(my_geneSet_custom) = c("custom_geneSet_sig_Genes")
  my_geneSet_custom$custom_geneSet_sig_Genes = my_geneSet

  my_geneSet_custom
  geneList
  fgseaRes <- fgsea(pathways = my_geneSet_custom, 
                    stats    = geneList,
                    minSize  = 10,
                    maxSize  = 500)
  partial_pval = fgseaRes$pval 
      
  paste0("For",n,"Top genes")
  print(my_geneSet_custom)
               
  list_of_partial_pval = append(list_of_partial_pval, partial_pval)
}

genes_topN_significant = genes_topN[1:13]

hgnc_symbol = anno$hgnc_symbol[match(genes_topN_significant, anno$ensembl_gene_id)]
hgnc_symbol = as.character(hgnc_symbol)
my_geneSet = hgnc_symbol[!is.na(hgnc_symbol)]


#Do ramdom selection 200 times with length 13 and save pvalues
list_of_partial_random_pval = list()
for (n in 1:200){
  random_n_length_gene_selection= names(sample(geneList, 13))

  fgseaRes <- fgsea(pathways = my_geneSet_custom, 
                      stats    = geneList,
                      minSize  = 10,
                      maxSize  = 500)
  print(n)
  partial_pval = fgseaRes$pval 
  list_of_partial_random_pval = append(list_of_partial_random_pval, partial_pval)
}
as.numeric(list_of_partial_random_pval)


#PERFORM SIGNIFICANDE PLOT FOR GENESETS OF EACH TRANSITION AND GENERATE A XY PLOT WITH RESULTS FOR EACH SET
#Bucle searcher
potential_genes_BA
potential_genes_CB
potential_genes_DC
potential_genes_ED

#For lvl2 
potential_genes_1_13
potential_genes_6_1
potential_genes_8_6
potential_genes_2_13

list_of_partial_pval = list()

for (n in 1:dim(potential_genes_1_13)[1]){
  genes_topN = rownames(potential_genes_1_13)[1:n]
  
  entrezgene = anno$entrezgene_id[match(genes_topN, anno$hgnc_symbol)]
  my_geneSet = as.character(entrezgene)
  my_geneSet = my_geneSet[!is.na(my_geneSet)]
   
  my_geneSet_custom = examplePathways[1]
  names(my_geneSet_custom) = c("custom_geneSet_sig_Genes")
  my_geneSet_custom$custom_geneSet_sig_Genes = my_geneSet

  my_geneSet_custom
  geneList
  fgseaRes <- fgsea(pathways = my_geneSet_custom, 
                    stats    = geneList,
                    minSize  = 10,
                    maxSize  = 500)
  partial_pval = fgseaRes$pval 
      
  paste0("For",n,"Top genes")
  print(my_geneSet_custom)
               
  list_of_partial_pval = append(list_of_partial_pval, partial_pval)
}

#RESULTS MATCHING AGAINS OUR TRANSITIONS (YES CHANGES APPEARS HERE)
#SET AB (4)             0.1303462

#SET BC  (21)           0.004547561

#SET CD  (14)           0.1007905

#SET DE  (14)           0.1428571

#FOR LVL2
#SET 1_13 (41)          0.009806458     

#SET 6_1  (96)          0.0001682029

#SET 8_6  (35)          0.0001608246

#SET 2_13  (75)         0.01051655 

#Bucle searcher FOR EACH STATIC POINT A B C D E
Idents(combined_seurat) = "contrasts_groups"

gene_CvsT_changes_A= FindMarkers(object = combined_seurat, ident.1 = "T_group_A", ident.2 = "C_group_A", only.pos = FALSE, min.pct = 0, logfc.threshold = 0.2)     #lvl1 with 0.2, lvl2 with 0.15
gene_CvsT_changes_A = gene_CvsT_changes_A[abs(gene_CvsT_changes_A$p_val_adj) < 0.05, ]
dim(gene_CvsT_changes_A)

gene_CvsT_changes_B= FindMarkers(object = combined_seurat, ident.1 = "T_group_B", ident.2 = "C_group_B", only.pos = FALSE, min.pct = 0, logfc.threshold = 0.2)
gene_CvsT_changes_B = gene_CvsT_changes_B[abs(gene_CvsT_changes_B$p_val_adj) < 0.05, ]
dim(gene_CvsT_changes_B)

gene_CvsT_changes_C= FindMarkers(object = combined_seurat, ident.1 = "T_group_C", ident.2 = "C_group_C", only.pos = FALSE, min.pct = 0, logfc.threshold = 0.2)
gene_CvsT_changes_C = gene_CvsT_changes_C[abs(gene_CvsT_changes_C$p_val_adj) < 0.05, ]
dim(gene_CvsT_changes_C)

gene_CvsT_changes_D= FindMarkers(object = combined_seurat, ident.1 = "T_group_D", ident.2 = "C_group_D", only.pos = FALSE, min.pct = 0, logfc.threshold = 0.2)
gene_CvsT_changes_D = gene_CvsT_changes_D[abs(gene_CvsT_changes_D$p_val_adj) < 0.05, ]
dim(gene_CvsT_changes_D)

gene_CvsT_changes_E= FindMarkers(object = combined_seurat, ident.1 = "T_group_E", ident.2 = "C_group_E", only.pos = FALSE, min.pct = 0, logfc.threshold = 0.2)
gene_CvsT_changes_E = gene_CvsT_changes_E[abs(gene_CvsT_changes_E$p_val_adj) < 0.05, ]
dim(gene_CvsT_changes_E)


list_of_partial_pval = list()

for (n in dim(gene_CvsT_changes_E)[1]:dim(gene_CvsT_changes_E)[1]){
  genes_topN = rownames(gene_CvsT_changes_E)[1:n]
  
  entrezgene = anno$entrezgene_id[match(genes_topN, anno$hgnc_symbol)]
  my_geneSet = as.character(entrezgene)
  my_geneSet = my_geneSet[!is.na(my_geneSet)]
   
  my_geneSet_custom = examplePathways[1]
  names(my_geneSet_custom) = c("custom_geneSet_sig_Genes")
  my_geneSet_custom$custom_geneSet_sig_Genes = my_geneSet

  my_geneSet_custom
  geneList
  fgseaRes <- fgsea(pathways = my_geneSet_custom, 
                    stats    = geneList,
                    minSize  = 10,
                    maxSize  = 500)
  partial_pval = fgseaRes$pval 
      
  paste0("For",n,"Top genes")
  print(my_geneSet_custom)
               
  list_of_partial_pval = append(list_of_partial_pval, partial_pval)
}

#RESULTS MATCHING AGAINS OUR static points (YES the trend is waht it should be) logChange 0.2
#A(21) 0.958498    

#B(8) 0.26         

#C(19) 0.33        

#D(10) 0.086629    

#E(10) 0.0057666

#FOR LVL2
#SET 13 (41)          0.9623016   

#SET 1  (96)          0.4074844

#SET 6  (35)          0.2743191

#SET 8 (75)           0.544592

#SET 2 (75)          0.9941176

#MAKE plot xy with both bar information and density line
values = c(0.9623016, 0.9941176, 0.4074844, 0.2743191 , 0.544592)
group = c('13', '2', '1', '6', '8')
pVAL = -log10(values)

data_summary = data.frame (Clusters =  c('13', '2', '1', '6', '8'),
                           pVAL = -log10(values)
                           )

#output_path = "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/validation_tfg_alejandro_data/"                          
output_path = "/ibex/scratch/projects/c2169/Xabier/Navarra/Neuro"
pdf(file=paste0(output_path,"GeneSets_significance_trend_by_static_clusters_CvsT_LVL2.pdf"), width=10, height=10)
ggplot(data_summary, aes(x = factor(Clusters, levels=c('13', '1', '6', '8', '2')), y = pVAL, fill = Clusters, colour = Clusters)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  labs(x="Developmental stages", y="-log10(P-Value)") +
  theme_bw() + 
  theme(text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.1, 0.75)) #+ 
  #geom_line(size = 1, color="blue", group = 1)
dev.off()

#SET 1_13 (41)          0.009806458     

#SET 6_1  (96)          0.0001682029

#SET 8_6  (35)          0.0001608246

#SET 2_13  (75)         0.01051655 
######
values = c(0.009806458, 0.01051655, 0.0001682029, 0.0001608246)
group = c('13-1', '13-2', '1-6', '6-8')
pVAL = -log10(values)

data_summary = data.frame (Transitions = c('13-1', '13-2', '1-6', '6-8'),
                           pVAL = -log10(values)
                           )

output_path = "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/validation_tfg_alejandro_data/"                          
pdf(file=paste0(output_path,"GeneSets_significance_trend_by_transitions_CvsT_LVL2.pdf"), width=10, height=10)
ggplot(data_summary, aes(x = factor(Transitions, levels = c('13-1', '13-2', '1-6', '6-8')), y = pVAL, fill = Transitions, colour = Transitions)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  labs(x="Developmental transitions", y="-log10(P-Value)") +
  theme_bw() + 
  theme(text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.1, 0.75)) + 
  geom_line(size = 1, color="blue", group = 1)
dev.off()



##################
#Compute a classifier model to allocate samples base on -our genes of interest/their genes of interest/random geneset of same length- as AorB
#Then do aucROC to check if we get more significance or less than a random gene selection

set.seed(42)

#Subset only the raiables of interest base on our DE gene list
table_dani = data.frame(top.table_dani)
rownames(table_dani) = make.unique(as.character(unlist(table_dani['hgnc_symbol'])))
overlap_geneSets = geneSets[[1]][geneSets[[1]] %in% rownames(table_dani)]

####
#for our list of genes relevants
#table_dani_base_our_genes = table_dani[overlap_geneSets,]
#FOR RANSOM SELECTION OF SAME LENTH N
#random_n_length_gene_selection= sample(rownames(table_dani), 802)
#table_dani_base_our_genes = table_dani[random_n_length_gene_selection,]
#For their genes of interest paper
their_genes = anno[(anno$ensembl_gene_id %in% list_sig_genes),]
their_genes = lapply(list(their_genes$hgnc_symbol), function(z){ z[!is.na(z) & z != ""]})
table_dani_base_our_genes = table_dani[unlist(their_genes),]
table_dani_base_our_genes = table_dani_base_our_genes[!grepl("NA.", rownames(table_dani_base_our_genes)),]
###

#keep only samples columns and Transpose it to generate the required data format
table_dani_base_our_genes = table_dani_base_our_genes[,3:20]
#This line only for check same size with their genes
#table_dani_base_our_genes = table_dani_base_our_genes[1:103,]
table_dani_base_our_genes = data.frame(t(table_dani_base_our_genes))
table_dani_base_our_genes['NA.'] = NULL #802 genes


#add a column with category
table_dani_base_our_genes['category'] = c(rep(0,8),rep(1,10))

#Lest generate the train and test datasets (they will be the same in this scenario)
#default_idx = sample(nrow(Default), 5000)
#default_trn = Default[default_idx, ]
#default_tst = Default[-default_idx, ]
default_trn = table_dani_base_our_genes[c(1,2,3,4,9,10,11,12,13),]
default_tst = table_dani_base_our_genes[c(5,6,7,8,14,15,16,17,18),]

#Train a logistic regresion model
model_1 = glm(category ~ ., data = default_trn, family = "binomial")

#Function to check different probability cut-offs
get_logistic_pred = function(mod, data, res = "y", pos = 1, neg = 0, cut = 0.5) {
  probs = predict(mod, newdata = data, type = "response")
  ifelse(probs > cut, pos, neg)
}

#Lets check low-med-high prob
test_pred_10 = get_logistic_pred(model_1, data = default_tst, res = "category", pos = "1", neg = "0", cut = 0.1)
test_pred_50 = get_logistic_pred(model_1, data = default_tst, res = "category", pos = "1", neg = "0", cut = 0.5)
test_pred_90 = get_logistic_pred(model_1, data = default_tst, res = "category", pos = "1", neg = "0", cut = 0.9)

#Evaluate accuracy, sensitivity, and specificity for the clussifiers
test_tab_10 = table(predicted = test_pred_10, actual = default_tst$category)
test_tab_50 = table(predicted = test_pred_50, actual = default_tst$category)
test_tab_90 = table(predicted = test_pred_90, actual = default_tst$category)

library(caret)
test_con_mat_10 = confusionMatrix(test_tab_10, positive = "1")
test_con_mat_50 = confusionMatrix(test_tab_50, positive = "1")
test_con_mat_90 = confusionMatrix(test_tab_90, positive = "1")

#prepare and print oputputs
metrics = rbind(
  
  c(test_con_mat_10$overall["Accuracy"], 
    test_con_mat_10$byClass["Sensitivity"], 
    test_con_mat_10$byClass["Specificity"]),
  
  c(test_con_mat_50$overall["Accuracy"], 
    test_con_mat_50$byClass["Sensitivity"], 
    test_con_mat_50$byClass["Specificity"]),
  
  c(test_con_mat_90$overall["Accuracy"], 
    test_con_mat_90$byClass["Sensitivity"], 
    test_con_mat_90$byClass["Specificity"])

)

rownames(metrics) = c("c = 0.10", "c = 0.50", "c = 0.90")
metrics

#our genes
#          Accuracy Sensitivity Specificity
#c = 0.10 0.6666667         0.8         0.5
#c = 0.50 0.6666667         0.8         0.5
#c = 0.90 0.6666667         0.8         0.5

#Random 802 selected genes
#          Accuracy Sensitivity Specificity
#c = 0.10 0.3333333         0.2         0.5
#c = 0.50 0.3333333         0.2         0.5
#c = 0.90 0.3333333         0.2         0.5

#Their genes
#          Accuracy Sensitivity Specificity
#c = 0.10 0.4444444         0.4         0.5
#c = 0.50 0.4444444         0.4         0.5
#c = 0.90 0.4444444         0.4         0.5


library(pROC)
test_prob = predict(model_1, newdata = default_tst, type = "response")

#RANDOM
test_roc = plot(roc(default_tst$category ~ test_prob), print.auc = TRUE, col = "black")
#OUR
test_roc = plot(roc(default_tst$category ~ test_prob), print.auc = TRUE, col = "red", print.auc.y = .2, add = TRUE)
#THEIR
test_roc = plot(roc(default_tst$category ~ test_prob), print.auc = TRUE, col = "orange", print.auc.y = .8, add = TRUE)
dev.off()

as.numeric(test_roc$auc)

#OUR SET OF GENES = 0.65 
#THEIR SET OF GENES = 0.425


############
#PROPORTION ANALYSIS ALONG NEUROGENESIS DEVELOPMENTAL BETWEEN CONTROL AND BETA CONDITIONS

############
# 1.Subset only relevant cells A-B-C-D-E
combined_seurat <- subset(combined_seurat, subset = contrasts_groups != "Rest")


############
# 2.Extract metadata  sample=protocol/condition=protocol_2/clusters=clusters
metadata <- combined_seurat@meta.data


############
# 3.Extract relevant information from cluster and condition into 2 columns
# Split the column using sub and gsub
metadata$clusters <- sub("^.*_", "", metadata$contrasts_groups)
metadata$barcode = rownames(metadata)

############
# 4.Find the minimum number of cells across all samples
library(dplyr)
min_cells <- metadata %>%
  group_by(protocol) %>%
  summarize(CellCount = n()) %>%
  pull(CellCount) %>%
  min()

# Subsample the same number of cells from each sample
subsampled_metadata <- metadata %>%
  group_by(protocol) %>%
  sample_n(min_cells) %>%
  ungroup()

# Check proportion by sample
table(subsampled_metadata$protocol)

############
# Calculate proportions for subsampled data for each cluster and condition
proportions_df <- subsampled_metadata %>%
  group_by(protocol, protocol_2, clusters) %>%
  summarize(CellCount = n(), .groups = "drop") %>%
  group_by(protocol, protocol_2) %>%
  mutate(Proportion = CellCount / sum(CellCount)) %>%
  ungroup()

head(proportions_df)


############
#5.Aggregate Across Samples
# Aggregate proportions
summary_df <- proportions_df %>%
  group_by(protocol_2, clusters) %>%
  summarize(
    MeanProportion = mean(Proportion),
    SD = sd(Proportion),
    .groups = "drop"
  )

summary_df


############
#6.Statistical Testing
# Create contingency table from subsampled data
library(tibble)
library(tidyr)
contingency_table <- subsampled_metadata %>%
  group_by(clusters, protocol_2) %>%
  summarize(CellCount = n(), .groups = "drop") %>%
  pivot_wider(names_from = protocol_2, values_from = CellCount, values_fill = 0) %>%
  column_to_rownames("clusters")

# Perform chi-square test
chisq_test <- chisq.test(as.matrix(contingency_table))
print(chisq_test)


############
# 6b. Ordinal Logistic Regression for Developmental Trajectory (A to E)
# Using the 'MASS' package for ordinal logistic regression (polr function)
library(MASS)

# Convert clusters into an ordered factor based on their developmental trajectory
subsampled_metadata$clusters <- factor(subsampled_metadata$clusters, 
                                       levels = c("A", "B", "C", "D", "E"), 
                                       ordered = TRUE)

# Fit the ordinal logistic regression model
model <- polr(clusters ~ protocol_2, data = subsampled_metadata, method = "logistic")

# Print model summary
summary(model)

# Check significance of the condition (protocol_2)
coef(model)
confint(model)
# You can check the p-values from the model summary, though note that the polr() function doesn't directly give p-values.
# You might need to use the Anova function from the car package if you need more detailed testing.


############
#7.Visualization
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/Neurogenesis_lvl_1/proportion_analysis/")
library(ggplot2)

summary_df$protocol_2 <- factor(summary_df$protocol_2, levels = c("control", "ba"))
pdf("cluster_proportions_by_condition.pdf", width = 8, height = 6)
ggplot(summary_df, aes(x = clusters, y = MeanProportion, fill = protocol_2)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = MeanProportion - SD, ymax = MeanProportion + SD),
    position = position_dodge(0.9), width = 0.2
  ) +
  theme_minimal() +
  labs(
    title = "Cluster Proportions by Condition",
    x = "Cluster",
    y = "Proportion"
  )
dev.off()

proportions_df$protocol_2 <- factor(proportions_df$protocol_2, levels = c("control", "ba"))
proportions_df$protocol <- factor(proportions_df$protocol, 
                                  levels = c(paste0("data_rna_control_", c("7", "13", "20"), "_days"),
                                             paste0("data_rna_ba_", c("7", "13", "20"), "_days")))
pdf("cluster_proportions_by_condition_2.pdf", width = 12, height = 6)  # Increase width
ggplot(proportions_df, aes(x = protocol, y = Proportion, fill = clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ protocol_2, scales = "free_x") +  # Allow x-axis to vary by facet
  theme_minimal() +
  labs(
    title = "Cluster Distributions Across Samples",
    x = "Sample",
    y = "Proportion"
  ) +
  # Rotate x-axis labels by 45 degrees
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10),   # Adjust facet label size
        axis.text = element_text(size = 8),      # Adjust axis text size
        axis.title = element_text(size = 10))    # Adjust axis title size
dev.off()





###############################
###############################
#
library(scProportionTest)

# Subset only unbiased selected cells from each sample to check also this analysis
combined_seurat@meta.data = metadata
combined_seurat <- subset(combined_seurat, cells = subsampled_metadata$barcode)


# Preparing the data using sc_utils on the combined_seurat object
prop_test <- scProportionTest::sc_utils(combined_seurat)

# Initialize an empty list to store the results
result_list <- list()

# Iterating over the conditions (control vs ba)
for (condition in c('control', 'ba')) {
  message("Running test for: ", condition)
  
  # Check if both control and the current condition are present in the data
  if (sum(combined_seurat@meta.data$protocol_2 == "control") > 0 & sum(combined_seurat@meta.data$protocol_2 == condition) > 0) {
    # Running the permutation test between control and the current condition
    result_list[[condition]] <- scProportionTest::permutation_test(
      prop_test, 
      cluster_identity = "clusters",  # Cluster identity is stored in clusters
      sample_1 = "control",           # Comparing control condition
      sample_2 = condition,           # Comparing ba condition for each iteration
      sample_identity = "protocol_2"  # Sample identity is based on protocol_2 (conditions)
    )
  } else {
    message("Skipping test for ", condition, " because there are no cells for this condition.")
  }
}


for (contrast in names(result_list)){
    result_list[[contrast]]@results$permutation$Group <- paste0(contrast, ' vs Control')
}

plotter <- as.data.frame(rbind(
    #reshape2::melt(result_list[['control']]@results$permutation[,c('obs_log2FD', 'FDR', 'clusters', 'Group')], id.vars = c('clusters', 'Group', 'obs_log2FD', 'FDR')),
    reshape2::melt(result_list[['ba']]@results$permutation[,c('obs_log2FD', 'FDR', 'clusters', 'Group')], id.vars = c('clusters', 'Group', 'obs_log2FD', 'FDR'))
	# ,
    # reshape2::melt(result_list[['X3']]@results$permutation[,c('obs_log2FD', 'FDR', 'clusters', 'Group')], id.vars = c('clusters', 'Group', 'obs_log2FD', 'FDR'))
))

#plotter$clusters <- factor(plotter$clusters, levels = group_ordered)
#plotter$clusters <- factor(stringr::str_wrap(gsub('_', ' ', plotter$clusters), width = 20), levels = stringr::str_wrap(gsub('_', ' ', group_ordered), width = 20))
plotter$clusters <- factor(plotter$clusters)
plotter$clusters <- factor(stringr::str_wrap(gsub('_', ' ', plotter$clusters), width = 20))

KO_colors <- c('ba vs Control' = "#FF6347") 

pdf("scProportionTest_package_plot.pdf", width = 8, height = 6)
ggplot(plotter, aes(x=obs_log2FD, y=clusters, fill=Group)) +
    geom_point(aes(size=FDR), alpha =0.8, color='black', shape=21)  +
    guides(size = guide_legend(reverse=F), fill = guide_legend(override.aes = list(size = 4, alpha=1)))+
    scale_size("P.value (FDR)",range=c(4,0.01), breaks=c(0.0015, 0.01, 0.05, 0.1, 0.5)) +
    scale_fill_manual(values=c("ba vs Control"  = paste(KO_colors['ba vs Control'])))+ 
    theme_bw() +
    geom_vline(xintercept = c(-0.2,0.2), linetype="dotted", color = "black", linewidth=0.5) + 
    theme(
		axis.title=element_text(family='Helvetica', size=7), 
		axis.text=element_text(family='Helvetica', size=5), 
		axis.title.y=element_blank(),
    legend.position="bottom",
    legend.spacing.x = unit(5, 'mm'),
    #legend.key.spacing.x = unit(0.5, 'mm'),
    #legend.key.spacing.y = unit(0.5, 'mm'),
		legend.text=element_text(size=6),
		legend.title=element_text(family='Helvetica', size=7, face='bold'),
		axis.ticks.y=element_blank()) + 
    xlab('log2Fd')
dev.off()




# Export raw counts for GEO upload
library(Matrix)

#Raw matrix
exp_matrix <- combined_seurat@assays$RNA@counts
writeMM(exp_matrix, file = "expression_matrix.mtx")
#Axis
write.csv(rownames(exp_matrix), "genes.csv", row.names = FALSE)
write.csv(colnames(exp_matrix), "barcodes.csv", row.names = FALSE)

# Export for shiny
saveRDS(combined_seurat, file = "/ibex/project/c2169/Xabier/Navarra/Neuro/data/combined_seurat.rds")

# Get how many cells for each time, condition and technology we have
# Summarize by time (vector_days), condition (protocol_2), and technology (protocol)
library(dplyr)

meta_df <- combined_seurat@meta.data

cell_summary <- meta_df %>%
  group_by(vector_days, protocol_2, protocol) %>%
  summarise(n_cells = n()) %>%
  arrange(vector_days, protocol_2)

print(cell_summary)








# Load scATACseq object to get how many cells for each time, condition and technology we have
library(GenomicRanges)
library(rtracklayer)   # For importing BED files

# List of per-sample peak files (edit this to match your actual paths)
peak_files <- list.files(path = "/ibex/scratch/projects/c2169/Xabier/Navarra/Neuro/peaks", pattern = "*.bed", full.names = TRUE)

# Import all peak files and store as GRanges objects
peak_list <- lapply(peak_files, import)

# Combine all GRanges into one
all_peaks <- do.call(c, peak_list)

# Merge overlapping peaks to generate consensus
combined.peaks <- reduce(all_peaks)

# Now you can continue with your existing code to write combined_peaks_neurogenesis.bed
df <- data.frame(
  seqnames = seqnames(combined.peaks),
  starts = start(combined.peaks) - 1,
  ends = end(combined.peaks),
  names = rep(".", length(combined.peaks)),
  scores = rep(".", length(combined.peaks)),
  strands = strand(combined.peaks)
)

write.table(df, file = "/ibex/scratch/projects/c2169/Xabier/Navarra/Neuro/peaks/combined_peaks_neurogenesis.bed", quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)






#################Raghad_H1
pdf(file="Raghad_H1_neurogenesis.pdf", width=10, height=6)
features = unique(c("H1F0"))
DotPlot(only_control_cells_subset, features = features, cols = c("blue", "yellow")) + RotatedAxis()  #only_control_cells_subset_neurons combined_seurat
dev.off()

library(Seurat)
library(dplyr)

# Store results here
de_results <- list()

# Create mapping vector
group_mapping <- c(
  "C_group_A" = "A", "T_group_A" = "A",
  "C_group_B" = "B", "T_group_B" = "B",
  "C_group_C" = "C", "T_group_C" = "C",
  "C_group_D" = "D", "T_group_D" = "D",
  "C_group_E" = "E", "T_group_E" = "E"
)
# Apply mapping to create new column
combined_seurat@meta.data$group_letter <- group_mapping[combined_seurat@meta.data$contrasts_groups]

combined_seurat@meta.data$group_letter <- factor(
  combined_seurat@meta.data$group_letter,
  levels = c("A", "B", "C", "D", "E")
)
clusters <- levels(combined_seurat@meta.data$group_letter)
Idents(combined_seurat) <- "group_letter"

de_results <- list()

for (clust in clusters) {
  # Subset to one cluster
  cluster_cells <- subset(combined_seurat, idents = clust)
  #
  # Only use cells with 'ba' or 'Control'
  cluster_cells <- subset(cluster_cells, subset = protocol_2 %in% c("ba", "control"))
  #
  # Set identity to protocol_2 for DE
  Idents(cluster_cells) <- "protocol_2"
  #
  # Run DE
  de <- FindMarkers(cluster_cells, ident.1 = "ba", ident.2 = "control", features = "H1F0", logfc.threshold = 0)
  #
  # Store result
  de_results[[clust]] <- de
}

library(Seurat)
library(dplyr)
library(ggplot2)

# --- Your DE results part ---
# Assume you have already run FindMarkers per cluster and stored results in `h1f0_de`

# Create significance labels function
get_signif_label <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

h1f0_de$signif_label <- sapply(h1f0_de$p_val_adj, get_signif_label)

# --- Prepare plot data ---
plot_data <- FetchData(combined_seurat, vars = c("H1F0", "protocol_2", "group_letter"))
plot_data <- plot_data[!is.na(plot_data$group_letter), ]

# Make sure protocol_2 is a factor with correct order
plot_data$protocol_2 <- factor(plot_data$protocol_2, levels = c("control", "ba"))

# Prepare label data for significance stars between "control" and "ba"
x_control <- which(levels(plot_data$protocol_2) == "control") # should be 1
x_ba <- which(levels(plot_data$protocol_2) == "ba")           # should be 2

label_data <- plot_data %>%
  group_by(group_letter) %>%
  summarise(ypos = max(H1F0, na.rm = TRUE) + 0.3, .groups = "drop") %>%
  left_join(h1f0_de %>% select(cluster, signif_label), by = c("group_letter" = "cluster")) %>%
  mutate(
    x_pos = (x_control + x_ba) / 2,
    x_control = x_control,
    x_ba = x_ba
  )

# --- Plot ---
ggplot(plot_data, aes(x = protocol_2, y = H1F0, fill = protocol_2)) +
  geom_violin(trim = TRUE) +
  facet_wrap(~ group_letter, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(strip.text = element_text(size = 14, face = "bold")) +

  # Add bracket lines
  geom_segment(data = label_data,
               aes(x = x_control, xend = x_ba, y = ypos - 0.05, yend = ypos - 0.05),
               inherit.aes = FALSE, size = 0.5) +
  geom_segment(data = label_data,
               aes(x = x_control, xend = x_control, y = ypos - 0.05, yend = ypos - 0.1),
               inherit.aes = FALSE, size = 0.5) +
  geom_segment(data = label_data,
               aes(x = x_ba, xend = x_ba, y = ypos - 0.05, yend = ypos - 0.1),
               inherit.aes = FALSE, size = 0.5) +

  # Add significance stars
  geom_text(data = label_data,
            aes(x = x_pos, y = ypos, label = signif_label),
            inherit.aes = FALSE,
            size = 6, fontface = "bold") +

  labs(title = "H1F0 expression by protocol and group_letter",
       y = "Expression Level", x = "Protocol")