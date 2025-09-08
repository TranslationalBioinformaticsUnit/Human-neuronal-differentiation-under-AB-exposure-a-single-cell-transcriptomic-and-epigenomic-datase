#####################################################################################################################
#
#####################################################################################################################
# STEP 1 - NEURONS LVL2 DEVELOPMENT SUBSET CONTROL CELLS - CONTRASTS
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/lvl_2_neurons/")

only_control_cells_subset_neurons = subset(x = only_control_cells_subset, idents=c("1","2","6","8","13","15"))
Idents(only_control_cells_subset_neurons) = "integrated_snn_res.0.6"

potential_markers_E_1_13 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "1", ident.2 = "13", only.pos = FALSE)
potential_markers_E_6_1  = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "6", ident.2 = "1", only.pos = FALSE)
potential_markers_E_8_6 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "8", ident.2 = "6", only.pos = FALSE)
potential_markers_E_2_13 = FindMarkers(object = only_control_cells_subset_neurons, ident.1 = "2", ident.2 = "13", only.pos = FALSE)

potential_markers_E_1_13 = potential_markers_E_1_13[abs(potential_markers_E_1_13$p_val_adj) < 0.05, ]
potential_markers_E_6_1 = potential_markers_E_6_1[abs(potential_markers_E_6_1$p_val_adj) < 0.05, ]
potential_markers_E_8_6 = potential_markers_E_8_6[abs(potential_markers_E_8_6$p_val_adj) < 0.05, ]
potential_markers_E_2_13 = potential_markers_E_2_13[abs(potential_markers_E_2_13$p_val_adj) < 0.05, ]

#Overall developmental 13-1-6-8  // 13-2
markers_E_1_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "1", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_6_1 = FindMarkers(object = only_control_cells_subset, ident.1 = "6", ident.2 = "1", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_8_6 = FindMarkers(object = only_control_cells_subset, ident.1 = "8", ident.2 = "6", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
markers_E_2_13 = FindMarkers(object = only_control_cells_subset, ident.1 = "2", ident.2 = "13", only.pos = FALSE, min.pct = 0, logfc.threshold = 0)

list_of_contrasts = list(markers_E_1_13,markers_E_6_1,markers_E_8_6,markers_E_2_13)
                      
names(list_of_contrasts) = c("markers_E_1_13","markers_E_6_1","markers_E_8_6","markers_E_2_13")

######################################################################################################################################################
#CLUSTERPROFILER
#GSEA for velocity_levels information ranking base on control vs BA contrast
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/lvl_2_neurons/GO_control/"
OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/Neurogenesis_lvl_2/"

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
names(geneList_list) = c("1-13","6-1","8-6","2-13")

geneList_list_tree = geneList_list
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseGO", OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff  = 0.01)
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseKEGG", pvalueCutoff  = 0.01)
geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseDO", pvalueCutoff  = 0.01)

geneList_list_tree <- pairwise_termsim(geneList_list_tree)
p1 <- treeplot(geneList_list_tree, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.0),
                    color = "NES", geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(2.5))

pdf(file= paste0("gsea_treeplot_C_do.pdf"), width = 20, height =30)
print(p1)
dev.off()

######################################################################################################################################################
#####DO limma (Ta-Tb) - (Ca-Cb)
library(limma)

#solo Tratamiento
Idents(combined_seurat) <- "protocol_2"
T_group_E_8 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(8) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_1 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(1) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_6 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(6) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_13 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(13) & combined_seurat@meta.data$protocol_2 == "ba"]
T_group_E_2 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(2) & combined_seurat@meta.data$protocol_2 == "ba"]


#solo Control
C_group_E_8 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(8) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_1 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(1) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_6 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(6) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_13 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(13) & combined_seurat@meta.data$protocol_2 == "control"]
C_group_E_2 = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$seurat_clusters == c(2) & combined_seurat@meta.data$protocol_2 == "control"]


combined_seurat@meta.data$contrasts_groups = "Rest"
combined_seurat@meta.data[T_group_E_8,]$contrasts_groups = "T_group_E_8"
combined_seurat@meta.data[T_group_E_1,]$contrasts_groups = "T_group_E_1"
combined_seurat@meta.data[T_group_E_6,]$contrasts_groups = "T_group_E_6"
combined_seurat@meta.data[T_group_E_13,]$contrasts_groups = "T_group_E_13"
combined_seurat@meta.data[T_group_E_2,]$contrasts_groups = "T_group_E_2"

combined_seurat@meta.data[C_group_E_8,]$contrasts_groups = "C_group_E_8"
combined_seurat@meta.data[C_group_E_1,]$contrasts_groups = "C_group_E_1"
combined_seurat@meta.data[C_group_E_6,]$contrasts_groups = "C_group_E_6"
combined_seurat@meta.data[C_group_E_13,]$contrasts_groups = "C_group_E_13"
combined_seurat@meta.data[C_group_E_2,]$contrasts_groups = "C_group_E_2"


g <- factor(combined_seurat@meta.data$contrasts_groups)
b <- factor(combined_seurat@meta.data$vector_days)

design <- model.matrix(~0 + g + b)
design <- design[,setdiff(colnames(design), "b4")] # get to full rank

exp_matrix = combined_seurat@assays$RNA@data
#y <- voom(exp_matrix, design, plot = F)
#corfit <- duplicateCorrelation(exp_matrix, design, block = combined_seurat@meta.data$vector_days)

fit <- lmFit(exp_matrix, design)#, block = combined_seurat@meta.data$vector_days, correlation = corfit$consensus)

contrast.matrix <- makeContrasts((gT_group_E_2 - gT_group_E_13) - (gC_group_E_2 - gC_group_E_13), levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)

summary(decideTests(fit2))
diff_expression_results <- topTreat(fit2, coef=1, n=Inf)

# STEP2 - Diff and significative changes of genes in contrast for Tx-Ty treatment strict derived changes
potential_genes_1_13 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_6_1 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_8_6 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
potential_genes_2_13 = diff_expression_results[abs(diff_expression_results$logFC) > 0.1 & diff_expression_results$adj.P.Val < 0.05, ]
write.csv(potential_genes_2_13, file = "Potential_markers_step_2_2_13.csv")

diff_expression_results_1_13 = diff_expression_results
diff_expression_results_6_1 = diff_expression_results
diff_expression_results_8_6 = diff_expression_results
diff_expression_results_2_13 = diff_expression_results


######################################################################################################################################################
#PREPARING gsea input
list_of_contrasts_step2 = list(diff_expression_results_1_13,diff_expression_results_6_1,diff_expression_results_8_6,diff_expression_results_2_13)
names(list_of_contrasts_step2) = c("diff_expression_results_1_13","diff_expression_results_6_1","diff_expression_results_8_6","diff_expression_results_2_13")

OutPath<- "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/lvl_2_neurons/GO_T/"

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
names(geneList_list) = c("1-13","6-1","8-6","2-13") #oldNames names(geneList_list) = c("BA","CB","DC","ED")

geneList_list_tree = geneList_list
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseGO", OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff  = 0.01)
#geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseKEGG", pvalueCutoff  = 0.01)
geneList_list_tree <- compareCluster(geneClusters=geneList_list_tree, fun = "gseDO", pvalueCutoff  = 0.05)

#eXPORT CSV WITH ENRICH UP/DOWN
write.csv(geneList_list_tree,paste0(OutPath,"gsea_treeplot_T","_kegg.csv"))

geneList_list_tree <- pairwise_termsim(geneList_list_tree)
p1 <- treeplot(geneList_list_tree, hclust_method = "average", showCategory=50, nCluster = 5, offset = rel(1.5),
                    color = "NES", geneClusterPanel = "heatMap",fontsize = 5,hexpand = .1, offset_tiplab = rel(4.0))

pdf(file= paste0("gsea_treeplot_T_do.pdf"), width = 20, height =30)
print(p1)
dev.off()




######################################################################################################################################################
# STEP 3 - create matrix with pertainti of genes on each category

names_sig_C_markers_E_1_13 = rownames(potential_markers_E_1_13)
names_sig_C_markers_E_6_1 = rownames(potential_markers_E_6_1)
names_sig_C_markers_E_8_6 = rownames(potential_markers_E_8_6)
names_sig_C_markers_E_2_13 = rownames(potential_markers_E_2_13)
write.table(potential_markers_E_2_13, file="potential_markers_E_2_13.csv", sep=";", row.names = TRUE)

names_sig_T_potential_E_1_13 = rownames(potential_genes_1_13)
names_sig_T_potential_E_6_1 = rownames(potential_genes_6_1)
names_sig_T_potential_E_8_6 = rownames(potential_genes_8_6)
names_sig_T_potential_E_2_13 = rownames(potential_genes_2_13)


out <- table(stack(mget(ls(pattern = "^names_sig_"))))
out <- out[, c("names_sig_C_markers_E_1_13", "names_sig_C_markers_E_6_1", "names_sig_C_markers_E_8_6", "names_sig_C_markers_E_2_13", "names_sig_T_potential_E_1_13", "names_sig_T_potential_E_6_1", "names_sig_T_potential_E_8_6", "names_sig_T_potential_E_2_13")]

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
aa = list(names_sig_T_potential_E_1_13,names_sig_T_potential_E_6_1,names_sig_T_potential_E_8_6,names_sig_T_potential_E_2_13)
aa_2 = list(potential_genes_1_13,potential_genes_6_1,potential_genes_8_6,potential_genes_2_13)

bb = list(names_sig_C_markers_E_1_13,names_sig_C_markers_E_6_1,names_sig_C_markers_E_8_6,names_sig_C_markers_E_2_13)
bb_2 = list(potential_markers_E_1_13,potential_markers_E_6_1,potential_markers_E_8_6,potential_markers_E_2_13)

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
  lista_gene_names = aa[[n]][!(aa[[n]] %in% c(names_sig_C_markers_E_1_13,names_sig_C_markers_E_6_1,names_sig_C_markers_E_8_6,names_sig_C_markers_E_2_13))]
  only_T = append(only_T, lista_gene_names)
}

only_T_conditions_all = unique(only_T)

#Only C
only_C_conditions_all_E_1_13 = names_sig_C_markers_E_1_13[!(names_sig_C_markers_E_1_13 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_6_1,names_sig_T_potential_E_8_6,names_sig_T_potential_E_2_13))]
only_C_conditions_all_E_6_1 = names_sig_C_markers_E_6_1[!(names_sig_C_markers_E_6_1 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_6_1,names_sig_T_potential_E_8_6,names_sig_T_potential_E_2_13))]
only_C_conditions_all_E_8_6 = names_sig_C_markers_E_8_6[!(names_sig_C_markers_E_8_6 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_6_1,names_sig_T_potential_E_8_6,names_sig_T_potential_E_2_13))]
only_C_conditions_all_E_2_13 = names_sig_C_markers_E_2_13[!(names_sig_C_markers_E_2_13 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_6_1,names_sig_T_potential_E_8_6,names_sig_T_potential_E_2_13))]
only_C_condition_all = unique(c(only_C_conditions_all_E_1_13,only_C_conditions_all_E_6_1,only_C_conditions_all_E_8_6,only_C_conditions_all_E_2_13))


######################################################################################################################################################
#REPEAT HEATMAP WITH EXPRESSION VALUES
#Table of binarized cascade activation + event(sum of cascade activations)
out_table_plot = out_table
out_table_plot$Events = NULL

#Fom binary to logFC matrix
out_table_plot[,1] = replace(out_table_plot[,1], out_table_plot[,1]==1, potential_markers_E_1_13[rownames(out_table_plot)[out_table_plot[,1]==1],]$avg_log2FC)
out_table_plot[,2] = replace(out_table_plot[,2], out_table_plot[,2]==1, potential_markers_E_6_1[rownames(out_table_plot)[out_table_plot[,2]==1],]$avg_log2FC)
out_table_plot[,3] = replace(out_table_plot[,3], out_table_plot[,3]==1, potential_markers_E_8_6[rownames(out_table_plot)[out_table_plot[,3]==1],]$avg_log2FC)
out_table_plot[,4] = replace(out_table_plot[,4], out_table_plot[,4]==1, potential_markers_E_2_13[rownames(out_table_plot)[out_table_plot[,4]==1],]$avg_log2FC)

out_table_plot[,5] = replace(out_table_plot[,5], out_table_plot[,5]==1, potential_genes_1_13[rownames(out_table_plot)[out_table_plot[,5]==1],]$logFC)
out_table_plot[,6] = replace(out_table_plot[,6], out_table_plot[,6]==1, potential_genes_6_1[rownames(out_table_plot)[out_table_plot[,6]==1],]$logFC)
out_table_plot[,7] = replace(out_table_plot[,7], out_table_plot[,7]==1, potential_genes_8_6[rownames(out_table_plot)[out_table_plot[,7]==1],]$logFC)
out_table_plot[,8] = replace(out_table_plot[,8], out_table_plot[,8]==1, potential_genes_2_13[rownames(out_table_plot)[out_table_plot[,8]==1],]$logFC)

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

#########
#HEATMAPS SUMMARY 1-c+t+only_T AND 2-only_C
genelist_both_and_T_only = c(both_conditions_all,only_T_conditions_all)
out_table_plot_genelist_both_and_T_only= out_table_plot[genelist_both_and_T_only,]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_1_13","C_6_1","C_8_6","C_2_13","T_1_13","T_6_1","T_8_6","T_2_13")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("1_13","6_1","8_6","2_13","1_13","6_1","8_6","2_13")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("1_13" = "#d01c8b",
                                           "6_1" = "#f1b6da",
                                           "8_6" = "#b8e186",
                                           "2_13" = "#4dac26")
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
    height = 15)
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
temp = potential_markers_E_1_13[order(potential_markers_E_1_13$avg_log2FC),]
temp_1_13 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_6_1[order(potential_markers_E_6_1$avg_log2FC),]
temp_6_1 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_8_6[order(potential_markers_E_8_6$avg_log2FC),]
temp_8_6 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_2_13[order(potential_markers_E_2_13$avg_log2FC),]
temp_2_13 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp_all = unique(c(temp_1_13,temp_6_1,temp_8_6,temp_2_13)) #116

only_C_condition_all_subset = only_C_condition_all[only_C_condition_all %in% temp_all]
out_table_plot_genelist_C_only= out_table_plot[only_C_condition_all_subset,]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 4)
tiempos=append(tiempos,rep("Treatment", 4))

#replace colnames by nice ones
table_plot_names =c("C_1_13","C_6_1","C_8_6","C_2_13","T_1_13","T_6_1","T_8_6","T_2_13")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("1_13","6_1","8_6","2_13","1_13","6_1","8_6","2_13")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("1_13" = "#d01c8b",
                                           "6_1" = "#f1b6da",
                                           "8_6" = "#b8e186",
                                           "2_13" = "#4dac26")
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

pdf(file = "significant_lvl_1_genelist_C_only.pdf",
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

out_table_plot[,1] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_E_13]))
out_table_plot[,2] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_E_1]))
out_table_plot[,3] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_E_6]))
out_table_plot[,4] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_E_8]))
out_table_plot[,5] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),C_group_E_2]))

out_table_plot[,6] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_E_13]))
out_table_plot[,7] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_E_1]))
out_table_plot[,8] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_E_6]))
out_table_plot[,9] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_E_8]))
out_table_plot[,10] = rowMeans(data.frame(combined_seurat@assays$RNA@data[rownames(out_table_plot),T_group_E_2]))

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
table_plot_names =c("C_13","C_1","C_6","C_8","C_2","T_13","T_1","T_6","T_8","T_2")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("13","1","6","8","2","13","1","6","8","2")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("13" = "#d01c8b",
                                           "1" = "#f1b6da",
                                           "6" = "#f7f7f7",
                                           "8" = "#b8e186",
                                           "2" = "#4dac26")
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
    height = 15)
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
temp_1_13 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_6_1[order(potential_markers_E_6_1$avg_log2FC),]
temp_6_1 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_8_6[order(potential_markers_E_8_6$avg_log2FC),]
temp_8_6 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp = potential_markers_E_2_13[order(potential_markers_E_2_13$avg_log2FC),]
temp_2_13 = c(rownames(head(temp,30)), rownames(tail(temp,30)))
temp_all = unique(c(temp_1_13,temp_6_1,temp_8_6,temp_2_13)) #116

only_C_condition_all_subset = only_C_condition_all[only_C_condition_all %in% temp_all]
out_table_plot_genelist_C_only= out_table_plot[only_C_condition_all_subset,]

# Annotations
#meta_data <- combined_seurat@meta.data
tiempos=rep("Control", 5)
tiempos=append(tiempos,rep("Treatment", 5))

#replace colnames by nice ones
table_plot_names =c("C_13","C_1","C_6","C_8","C_2","T_13","T_1","T_6","T_8","T_2")

#Create metadata
meta_data_plot = data.frame(Neurogenesis=table_plot_names, Condition=tiempos)
rownames(meta_data_plot)=meta_data_plot$Neurogenesis
meta_data_plot$Neurogenesis = NULL
#Second annotation
meta_data_plot$Transition = c("13","1","6","8","2","13","1","6","8","2")

annotation_colors <- list("Condition" = c("Control" = "#d8b365",
                                         "Treatment" = "#5ab4ac"),
                          "Transition" = c("13" = "#d01c8b",
                                           "1" = "#f1b6da",
                                           "6" = "#f7f7f7",
                                           "8" = "#b8e186",
                                           "2" = "#4dac26")
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

cds_short_names_df = data.frame(gene_short_name = rownames(only_control_cells_subset_neurons))
rownames(cds_short_names_df) = rownames(only_control_cells_subset_neurons)

cds <- new_cell_data_set(only_control_cells_subset_neurons@assays$RNA@data,
                         cell_metadata = only_control_cells_subset_neurons@meta.data,
                         gene_metadata = cds_short_names_df)
cds@int_colData$reducedDims@listData[["UMAP"]] = only_control_cells_subset_neurons@reductions[["umap"]]@cell.embeddings

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
           color_cells_by = "integrated_snn_res.0.6",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)



# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="13"){
  cell_ids <- which(colData(cds)[, "integrated_snn_res.0.6"] == time_bin)
  
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
Idents(only_control_cells_subset_neurons) = only_control_cells_subset_neurons@meta.data$integrated_snn_res.0.6
Idents(only_control_cells_subset_neurons) <- factor(Idents(only_control_cells_subset_neurons), levels= c("13","15","1","2","6","8"))

#level 1: Hes5/PAX6/DCX/CALB1
#level 2: NES/Sox11/MAP2/DCX/STMN1/SYP/CALB1

pdf(file="Characterization_neurogenesis.pdf", width=10, height=6)
features = unique(c("NES","SOX11","MAP2","DCX","STMN1","SYP","CALB1"))
DotPlot(only_control_cells_subset_neurons, features = features, cols = c("blue", "yellow")) + RotatedAxis()  #only_control_cells_subset_neurons combined_seurat
dev.off()

#Plot Genes from Maite OPTION _1 AND _2
pdf(file="Characterization_neurogenesis_2_1.pdf", width=12, height=12)
p1 = DimPlot(only_control_cells_subset_neurons, reduction = 'umap', group.by = "integrated_snn_res.0.6",label = TRUE)

CombinePlots(list(p1))
dev.off()

pdf(file="Characterization_neurogenesis_2_2.pdf", width=12, height=12)
p2 = FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("NES"),
      ncol = 1
      )
CombinePlots(list(p2))
dev.off()

pdf(file="Characterization_neurogenesis_2_3.pdf", width=12, height=12)
p3 = FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("SOX11"),
      ncol = 1
      )
CombinePlots(list(p3))
dev.off()

pdf(file="Characterization_neurogenesis_2_4.pdf", width=12, height=12)
p4 = FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("MAP2"),
      ncol = 1
      )
CombinePlots(list(p4))
dev.off()

pdf(file="Characterization_neurogenesis_2_5.pdf", width=12, height=12)
p5 = FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("DCX"),
      ncol = 1
)
CombinePlots(list(p5))
dev.off()

pdf(file="Characterization_neurogenesis_2_6.pdf", width=12, height=12)
p6 = FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("STMN1"),
      ncol = 1
)
CombinePlots(list(p6))
dev.off()

pdf(file="Characterization_neurogenesis_2_7.pdf", width=12, height=12)
p7 = FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("SYP"),
      ncol = 1
)
CombinePlots(list(p7))
dev.off()

pdf(file="Characterization_neurogenesis_2_8.pdf", width=12, height=12)
p8= FeaturePlot(
      object = only_control_cells_subset_neurons,
      features = c("CALB1"),
      ncol = 1
)
CombinePlots(list(p8))
dev.off()

#######
#MOUNT SIGNIFICANT THINGS TABLES 
#######
lista_all[[1]]
lista_all[[2]]
lista_all[[3]]
lista_all[[4]]

#GENES ONLY IN CONTROL
lista_gene_names = names_sig_C_markers_E_2_13[!(names_sig_C_markers_E_2_13 %in% c(names_sig_T_potential_E_1_13,names_sig_T_potential_E_6_1,names_sig_T_potential_E_8_6,names_sig_T_potential_E_2_13))]
potential_markers_E_2_13[lista_gene_names,]

#GENES ONLY IN TREATMENT
lista_gene_names = names_sig_T_potential_E_2_13[!(names_sig_T_potential_E_2_13 %in% c(names_sig_C_markers_E_1_13,names_sig_C_markers_E_6_1,names_sig_C_markers_E_8_6,names_sig_C_markers_E_2_13))]
potential_genes_2_13[lista_gene_names,]


names_sig_C_markers_E_1_13 = rownames(potential_markers_E_1_13)
names_sig_C_markers_E_6_1 = rownames(potential_markers_E_6_1)
names_sig_C_markers_E_8_6 = rownames(potential_markers_E_8_6)
names_sig_C_markers_E_2_13 = rownames(potential_markers_E_2_13)

names_sig_T_potential_E_1_13 = rownames(potential_genes_1_13)
names_sig_T_potential_E_6_1 = rownames(potential_genes_6_1)
names_sig_T_potential_E_8_6 = rownames(potential_genes_8_6)
names_sig_T_potential_E_2_13 = rownames(potential_genes_2_13)


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

setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/Neurogenesis_lvl_2/pathways/")

for (n in 1:length(data.frame(ego)$ID)){
  pathview(gene.data  = as.matrix(my_data_corrected[1:8]),
                       pathway.id = data.frame(ego)$ID[n],
                       species    = "hsa",
                       limit      = c(-0.10, 0.10))
}





















