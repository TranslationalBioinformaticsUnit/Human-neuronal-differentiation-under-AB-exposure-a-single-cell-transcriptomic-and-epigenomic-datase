############
# LIBRARIES #
############

library(readxl)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(biomaRt)
library(xlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(RColorBrewer)
library(enrichplot)
library(pathview)
library(tidyr)
library(BuenColors)

options(java.parameters = "-Xmx1024m")
set.seed(23) # For reproducibility

##########################################################################
# 1. LOAD RAW COUNT DATA
##########################################################################

setwd("your_path_here")

df <- read.table("raw_count_slide3.txt", header = TRUE, sep = "\t")
df <- unique(df)
rownames(df) <- df[,1]
df <- df[,-1]

colnames(df) <- c(paste("AD", 1:8, sep="_"), paste("Control", 1:10, sep="_"))
counts_matrix <- as.matrix(df)
summary(counts_matrix)

##########################################################################
# 2. DATA PROCESSING AND FILTERING
##########################################################################

# Keep genes with counts >=1 in at least 3 samples
counts_matrix1 <- counts_matrix[rowSums(counts_matrix >= 1) > 2, ]
head(counts_matrix1)
length(rownames(counts_matrix1))

# Quality control: raw counts boxplot
setwd("your_path_here/Output")

tiff("Raw_counts_boxplot.tiff", units="cm", width=15, height=15, res=300)
color.palette <- paletteer_d("ggsci::nrc_npg")
boxplot(log2(counts_matrix1),
        las  = 2,
        col  = color.palette[2],
        ylab = 'Raw Counts',
        main = 'Raw counts boxplot')
dev.off()

##########################################################################
# 3. NORMALIZATION (TMM)
##########################################################################

# Design matrix
design_matrix <- data.frame(
  combine = c(rep("AD",8), rep("Control",10)),
  gender  = c("Female","Female","Female","Female","Male","Male","Male","Female",
              "Female","Female","Female","Female","Male","Male","Female",
              "Female","Male","Male","Male"),
  age     = c("88","95","95","100","99","83","90","84","87","80",
              "84","77","55","72","78","83","80","74")
)
rownames(design_matrix) <- colnames(counts_matrix1)
comb <- as.factor(design_matrix$combine)
gender <- as.factor(design_matrix$gender)
age <- factor(design_matrix$age, levels = c(55, 72, 74, 77, 78, 80, 83, 84, 87, 88, 90, 95, 99, 100))

# Create DGEList
d0 <- DGEList(counts = counts_matrix1, group = comb)
d0$samples$comb <- comb
d0$samples$gender <- gender
d0$samples$age <- age

design <- model.matrix(~0 + comb + gender + age)
head(design, 20)

# Filter lowly expressed genes
keep.exprs <- filterByExpr(d0, group = gender)
d0 <- d0[keep.exprs, , keep.lib.sizes = FALSE]

# TMM normalization and voom transformation
d0 <- calcNormFactors(d0, method = "TMM")
v <- voom(d0, design, plot = TRUE)
norm_counts <- v$E

##########################################################################
# 4. QUALITY CONTROL: BOXPLOT, DENSITY, PCA, MDS, HEATMAP
##########################################################################

# Boxplot after normalization
tiff("Norm_boxplot.tiff", units="cm", width=15, height=15, res=300)
boxplot(norm_counts,
        las = 3,
        col = color.palette[2],
        ylab = 'Normalized expression levels',
        main = 'Normalized counts boxplot')
dev.off()

# Density plot
tiff("Density.tiff", units="in", width=10, height=10, res=300)
plot(density(norm_counts[,1]), lwd = 5, col = sample(color.palette,1),
     xlab = "Expression values", ylab = "Density",
     main = "Distribution of transformed data", ylim=c(0,0.3), xlim=c(-3,17))
for (i in 2:ncol(norm_counts)) {
  lines(density(norm_counts[,i]), lwd = 5, col = sample(color.palette,1))
}
dev.off()

# PCA
pca <- prcomp(t(norm_counts))
pr <- summary(pca)$importance[,1:5]
mplot <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], groups = design_matrix$combine)
tiff("PCA.tiff", units="in", width=7, height=5, res=300)
pca1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = groups, label = rownames(mplot))) +
  geom_point(size=4) +
  scale_color_manual(values=color.palette) +
  labs(title="PCA - Normalized counts",
       x=paste("PC1 (", round(pr[2,1]*100,0),"%)", sep=""),
       y=paste("PC2 (", round(pr[2,2]*100,0),"%)", sep="")) +
  geom_text(aes(label=rownames(design_matrix)), hjust=-0.1, vjust=-0.1) +
  theme_classic()
pca1
dev.off()

# MDS by condition
d <- dist(t(norm_counts))
fit <- cmdscale(d, eig = TRUE, k = 2)
x <- fit$points[,1]
y <- fit$points[,2]
mplot <- data.frame(Coordinate_1 = x, Coordinate_2 = y, Clusters = comb)
tiff("MDS_comb.tiff", units="in", width=7, height=5, res=300)
mds1 <- ggplot(mplot, aes(Coordinate_1, Coordinate_2, color=Clusters, label=row.names(mplot))) +
  geom_point(size=3) +
  scale_color_manual(values=c("Orange", "Purple")) +
  labs(title="MDS - Normalized counts", x="Coordinate 1", y="Coordinate 2")
mds1
dev.off()

##########################################################################
# 5. LINEAR MODEL AND DIFFERENTIAL EXPRESSION
##########################################################################

design <- model.matrix(~0 + comb)
fit <- lmFit(v, design)
contr <- makeContrasts(AD_vs_Control = combAD - combControl, levels = design)
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

results <- matrix()
list_sig_genes <- vector()
for(i in 1:ncol(contr)){
  res_cont <- topTable(tmp, coef=i, number=nrow(norm_counts), sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val < 0.05 & res_cont$logFC > 1)
  sig_genes_down <- which(res_cont$adj.P.Val < 0.05 & res_cont$logFC < -1)
  
  res_cont$sig_genes <- 0
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  
  colnames(res_cont) <- paste(colnames(contr)[i], colnames(res_cont), sep=":")
  results <- cbind(results, res_cont)
  list_sig_genes <- c(list_sig_genes, rownames(res_cont)[c(sig_genes_up, sig_genes_down)])
}

results_MM <- results[,-1]
list_sig_genes_MM <- unique(list_sig_genes)

##########################################################################
# 6. HEATMAP OF SIGNIFICANT GENES
##########################################################################

hm_sig <- norm_counts[list_sig_genes_MM, ]
hm_sig_scale <- as.matrix(t(scale(t(hm_sig))))

color_HM <- jdb_palette("cyan_violet")
color_corona <- jdb_palette("corona")

tiff("HeatMap.tiff", units="in", width=8, height=6, res=300)
ha1 <- HeatmapAnnotation(
  groups = design_matrix$combine,
  gender = design_matrix$gender,
  age = age,
  col = list(
    groups = c("AD"=color_corona[1], "Control"=color_corona[2]),
    age = setNames(color_HM, levels(age)),
    gender = c("Female"=color_corona[3], "Male"=color_corona[4])
  ),
  show_annotation_name = TRUE
)
Heatmap(hm_sig_scale, name = "value", top_annotation = ha1, show_column_dend = TRUE, show_row_names = FALSE)
dev.off()

##########################################################################
# 7. GENE SET ENRICHMENT ANALYSIS (GSEA)
##########################################################################

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
anno <- getBM(
  values = rownames(results_MM),
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  mart = human
)

top.table <- merge(cbind(Id=rownames(results_MM), results_MM), anno,
                   by.x="Id", by.y="ensembl_gene_id")

geneList <- top.table %>%
  dplyr::select(entrezgene_id, `AD_vs_Control:logFC`, `AD_vs_Control:adj.P.Val`) %>%
  filter(!is.na(entrezgene_id) & `AD_vs_Control:adj.P.Val` < 0.05) %>%
  arrange(desc(`AD_vs_Control:logFC`))

vector <- setNames(geneList$`AD_vs_Control:logFC`, geneList$entrezgene_id)

kk <- gseKEGG(
  geneList = vector,
  organism = "hsa",
  minGSSize = 30,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

ridgeplot(kk)
browseKEGG(kk, 'hsa04080')
pathview(gene.data = vector, pathway.id = "hsa05022", species = "hsa")

##########################################################################
# 8. CORRELATION WITH EXTERNAL MARKERS
##########################################################################

df_markers <- read.table("markers_for_external_d20.txt", header = TRUE, sep = ",", dec = ".")
rownames(df_markers) <- df_markers[,1]
df_markers <- df_markers[,-1]

top.table <- top.table %>% distinct(hgnc_symbol, .keep_all = TRUE)

# Merge with external marker data
i <- top.table %>% filter(`AD_vs_Control:adj.P.Val` < 0.05) %>% rename(genes = hgnc_symbol)
j <- df_markers %>% filter(p_val_adj < 0.05)
j$genes <- rownames(j)
ij <- merge(i, j, by = "genes")
ij <- ij %>% dplyr::select(`AD_vs_Control:logFC`, p_val_adj)
