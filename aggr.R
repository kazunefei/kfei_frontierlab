library(Seurat)
library(hdf5r)
library(dplyr)
library(tidyr) 
library(stringr)
library(patchwork)
setwd("C:/Users/okadalab/Documents/Mouse_scRNA-seq_data/aggr_0_50")
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
counts$`Gene Expression`
counts$Peaks

pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc@assays$RNA@counts
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 2500 & nFeature_RNA < 10500 & percent.mt < 45)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap",label=T)
saveRDS(pbmc, file = "neuron.rds")
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,  fc.threshold = 0.25)
clustermarker=pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

#adding a column for the sample number
sample_number <- c(1)
pbmc.test <- cbind(pbmc@meta.data, sample_number)
for (i in 1:nrow(pbmc.test)) {
  pbmc.test$sample_number <- replace(pbmc.test$sample_number, i, as.numeric((str_sub(row.names(pbmc.test[i,]), -1))))
}
pbmc@meta.data$sample_number=pbmc.test$sample_number

#adding a column for the sample number + cluster number
sample_cluster <- c("cluster")
pbmc.test <- cbind(pbmc@meta.data, sample_cluster)
for (i in 1:nrow(pbmc.test)) {
  pbmc.test$sample_cluster <- replace(pbmc.test$sample_cluster, i, paste(toString(pbmc$sample_number[i]), toString(pbmc.test$seurat_clusters[i]),  sep="_"))
}
pbmc@meta.data$sample_cluster=pbmc.test$sample_cluster

#improved version
gene_reg <- c()
for (i in 0:13) {
  marker <- (FindMarkers(pbmc, group.by = "sample_cluster", ident.1 = paste("1", toString(i), sep="_"), ident.2 = paste("2", toString(i), sep="_")))
  marker$cluster<-i
  gene_reg<-rbind(gene_reg, marker)
}
gene_reg <- dplyr::filter(gene_reg, p_val_adj<0.05)
for (i in 1:nrow(gene_reg)) {
  if(gene_reg$avg_log2FC[i]>0){
    gene_reg$regulation <- replace(gene_reg$regulation, i,"downregulation")
  }else{
    gene_reg$regulation <- replace(gene_reg$regulation,i,"upregulation")
  }
}
gene_regulation <- gene_reg
#end


nearest_cell<-FindNeighbors(pbmc,return.neighbor=TRUE)
near_cells<-TopNeighbors(nearest_cell[["RNA.nn"]],cell=colnames(pbmc)[1],n=10)

cell = c()
neighbors = c()

for (i in 1:nrow(pbmc@meta.data)) {
  near_cells<-TopNeighbors(nearest_cell[["RNA.nn"]],cell=colnames(pbmc)[i],n=20)
  if (sum(str_detect(neighbors, rownames(pbmc@meta.data)[i]))==0 | length(neighbors)==0){
    cell <- append(cell, rownames(pbmc@meta.data)[i])
    neighbors <- append(neighbors, toString(near_cells))
  }
}


cell_neighbors <- as.data.frame(cbind(cell, neighbors))

nearest5_1 <- c()
nearest5_2 <- c()

for (i in 1:nrow(cell_neighbors)){
  temp5_1 <- list()
  temp5_2 <- list()
  temp_list <- as.list(strsplit(cell_neighbors$neighbors, ", ")[[i]])
  for (j in 1:length(temp_list)){
    if (length(temp5_1)<5){
      if (grepl("-1", temp_list[j])){
        temp5_1 <- append(temp5_1, temp_list[j])
      }
    }
    if (length(temp5_2)<5){
      if (grepl("-2", temp_list[j])){
        temp5_2 <- append(temp5_2, temp_list[j])
      }
    }
  }
  nearest5_1 <- append(nearest5_1, toString(temp5_1))
  nearest5_2 <- append(nearest5_2, toString(temp5_2))
}

cell_neighbors <- cbind(cell_neighbors, nearest5_1, nearest5_2)
assayData <- as.data.frame(GetAssayData(object = pbmc, slot = "counts"))

list_of_df <- list()
for (i in 1:nrow(cell_neighbors)) {
  temp_list <- as.list(strsplit(cell_neighbors$nearest5_1, ", ")[[i]])
  
  temp_df <- assayData[,colnames(assayData) %in% temp_list]
  df <- as.data.frame(apply(as.matrix(temp_df),1,mean))
  colnames(df)<-cell_neighbors$cell[i]
  list_of_df <- append(list_of_df, df)
}
averages_1 <- bind_rows(list_of_df, .id = "column_label")
rownames(averages_1) <- rownames(assayData)

list_of_df <- list()
for (i in 1:nrow(cell_neighbors)) {
  temp_list <- as.list(strsplit(cell_neighbors$nearest5_2, ", ")[[i]])
  
  temp_df <- assayData[,colnames(assayData) %in% temp_list]
  df <- as.data.frame(apply(as.matrix(temp_df),1,mean))
  colnames(df)<-cell_neighbors$cell[i]
  list_of_df <- append(list_of_df, df)
}
averages_2 <- bind_rows(list_of_df, .id = "column_label")
rownames(averages_2) <- rownames(assayData)

#-----------------------------------

# using sc-type to try to annotate cell types
# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Brain") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# assign cell types
scRNAseqData = readRDS(gzcon(url('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/exampleData.RDS'))); #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)

#AUG
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

dev.off()

#gene ids are "symbol"

#"RNAseq analysis | Gene ontology (GO)in R - Sanbomics (youtube.com)"
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

#upregulated genes
for (i in 0:13) {
  cluster_GO <- dplyr::filter(gene_regulation, cluster == i)
  genes_to_test <- rownames(cluster_GO[cluster_GO$avg_log2FC<0,])
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1)
  if(nrow(as.data.frame(GO_results))!=0){
    png(paste(toString(i), "upregulated_results_plot.png", sep = "_"), res = 250, width = 2000, height = 2750)
    print(plot(barplot(GO_results, showCategory = 20)))
    dev.off()
    png(paste(toString(i), "upregulated_results_plot_GO.png", sep = "_"), res = 250, width = 5000, height = 5000)
    print((goplot(GO_results, showCategory = 20, )))
    dev.off()
  }else{
  }
}
#downregulated
for (i in 0:13) {
  cluster_GO <- dplyr::filter(gene_regulation, cluster == i)
  genes_to_test <- rownames(cluster_GO[cluster_GO$avg_log2FC>0,])
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1)
  if(nrow(as.data.frame(GO_results))!=0){
    png(paste(toString(i), "downregulated_results_plot.png", sep = "_"), res = 250, width = 2000, height = 2750)
    print(plot(barplot(GO_results, showCategory = 20)))
    dev.off()
    png(paste(toString(i), "downregulated_results_plot_GO.png", sep = "_"), res = 250, width = 5000, height = 5000)
    print((goplot(GO_results, showCategory = 20, )))
    dev.off()
  }else{
  }
}

cluster_GO <- dplyr::filter(gene_regulation, cluster == 1)
genes_to_test <- rownames(cluster_GO[cluster_GO$avg_log2FC<0,])
GOresults <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1)
plot(barplot(GOresults, showcategory=20,))

#gsea analysis
library(clusterProfiler)
library(org.Mm.eg.db)


#deseq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

#getting assay data (counts) and meta_data in a dataframe
assayData <- as.data.frame(GetAssayData(object = pbmc))
meta_data <- as.data.frame(pbmc@meta.data)
meta_data <- subset(meta_data, select = -c(id, name, chapters))

#-------------------------------------------------------------------------------------------------------------------------

#filtering pbmc to make smaller pbmc object
example <- GetAssayData(pbmc, slot="counts", assay="RNA")   
genes.percent.expressed <- rowMeans(example>0 )*100  

genes.filter <- names(genes.percent.expressed[genes.percent.expressed>5])  #select genes expressed in at least 15% of cells
counts.sub <- example[genes.filter,]
pbmc_filtered <- CreateSeuratObject(counts=counts.sub)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#clustering on filtered pbmc
pbmc_filtered[["percent.mt"]] <- PercentageFeatureSet(pbmc_filtered, pattern = "^mt-")
VlnPlot(pbmc_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc_filtered <- subset(pbmc_filtered, subset = nFeature_RNA > 2500 & nFeature_RNA < 10500 & percent.mt < 45)
pbmc_filtered <- NormalizeData(pbmc_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_filtered <- FindVariableFeatures(pbmc_filtered, selection.method = "vst", nfeatures = 1000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc_filtered), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc_filtered)
pbmc_filtered <- ScaleData(pbmc_filtered, features = all.genes)
pbmc_filtered <- RunPCA(pbmc_filtered, features = VariableFeatures(object = pbmc_filtered))
VizDimLoadings(pbmc_filtered, dims = 1:2, reduction = "pca")
DimPlot(pbmc_filtered, reduction = "pca")
pbmc_filtered <- JackStraw(pbmc_filtered, num.replicate = 100)
pbmc_filtered <- ScoreJackStraw(pbmc_filtered, dims = 1:20)
pbmc_filtered <- FindNeighbors(pbmc_filtered, dims = 1:10)
pbmc_filtered <- FindClusters(pbmc_filtered, resolution = 0.5)
head(Idents(pbmc_filtered), 5)
pbmc_filtered <- RunUMAP(pbmc_filtered, dims = 1:10)
DimPlot(pbmc_filtered, reduction = "umap",label=T)
saveRDS(pbmc_filtered, file = "neuron_filtered.rds")
cluster2.markers <- FindMarkers(pbmc_filtered, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
pbmc_filtered.markers <- FindAllMarkers(pbmc_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
clustermarker=pbmc_filtered.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

#adding a column for the sample number
sample_number <- c(1)
pbmc_filtered.test <- cbind(pbmc_filtered@meta.data, sample_number)
for (i in 1:nrow(pbmc_filtered.test)) {
  pbmc_filtered.test$sample_number <- replace(pbmc_filtered.test$sample_number, i, as.numeric((str_sub(row.names(pbmc_filtered.test[i,]), -1))))
}
pbmc_filtered@meta.data$sample_number=pbmc_filtered.test$sample_number

#adding a column for the sample number + cluster number
sample_cluster <- c("cluster")
pbmc_filtered.test <- cbind(pbmc_filtered@meta.data, sample_cluster)
for (i in 1:nrow(pbmc_filtered.test)) {
  pbmc_filtered.test$sample_cluster <- replace(pbmc_filtered.test$sample_cluster, i, paste(toString(pbmc_filtered$sample_number[i]), toString(pbmc_filtered.test$seurat_clusters[i]),  sep="_"))
}
pbmc_filtered@meta.data$sample_cluster=pbmc_filtered.test$sample_cluster

#improved version
gene_reg <- c()
for (i in 0:13) {
  marker <- (FindMarkers(pbmc_filtered, group.by = "sample_cluster", ident.1 = paste("1", toString(i), sep="_"), ident.2 = paste("2", toString(i), sep="_")))
  marker$cluster<-i
  gene_reg<-rbind(gene_reg, marker)
}
gene_reg <- dplyr::filter(gene_reg, p_val_adj<0.05)
for (i in 1:nrow(gene_reg)) {
  if(gene_reg$avg_log2FC[i]>0){
    gene_reg$regulation <- replace(gene_reg$regulation, i,"downregulation")
  }else{
    gene_reg$regulation <- replace(gene_reg$regulation,i,"upregulation")
  }
}
gene_regulation <- gene_reg
#end

meta_data_filtered <- pbmc_filtered@meta.data
meta_data_filtered_target <- subset(meta_data_filtered, select = -c(orig.ident, 
                                                                    nCount_RNA, 
                                                                    nFeature_RNA, 
                                                                    percent.mt, 
                                                                    RNA_snn_res.0.5, 
                                                                    sample_number, 
                                                                    sample_cluster))
assayData_filtered <- as.data.frame(GetAssayData(object = pbmc_filtered, slot = "counts"))

write.csv(meta_data_filtered_target, "C:/Users/okadalab/Documents/Mouse_scRNA-seq_data/aggr_0_50/target.counts/filtered_targets.csv")
write.csv(assayData_filtered, "C:/Users/okadalab/Documents/Mouse_scRNA-seq_data/aggr_0_50/target.counts/filtered_counts.csv")

#computing most similar cells
nearest_cell<-FindNeighbors(pbmc_filtered,return.neighbor=TRUE)
near_cells<-TopNeighbors(nearest_cell[["RNA.nn"]],cell=colnames(pbmc_filtered)[1],n=10)

cell = c()
neighbors = c()

for (i in 1:nrow(pbmc_filtered@meta.data)) {
  near_cells<-TopNeighbors(nearest_cell[["RNA.nn"]],cell=colnames(pbmc_filtered)[i],n=20)
  if (sum(str_detect(neighbors, rownames(pbmc_filtered@meta.data)[i]))==0 | length(neighbors)==0){
    cell <- append(cell, rownames(pbmc_filtered@meta.data)[i])
    neighbors <- append(neighbors, toString(near_cells))
  }
}


cell_neighbors <- as.data.frame(cbind(cell, neighbors))

nearest5_1 <- c()
nearest5_2 <- c()

for (i in 1:nrow(cell_neighbors)){
  temp5_1 <- list()
  temp5_2 <- list()
  temp_list <- as.list(strsplit(cell_neighbors$neighbors, ", ")[[i]])
  for (j in 1:length(temp_list)){
    if (length(temp5_1)<5){
      if (grepl("-1", temp_list[j])){
        temp5_1 <- append(temp5_1, temp_list[j])
      }
    }
    if (length(temp5_2)<5){
      if (grepl("-2", temp_list[j])){
        temp5_2 <- append(temp5_2, temp_list[j])
      }
    }
  }
  nearest5_1 <- append(nearest5_1, toString(temp5_1))
  nearest5_2 <- append(nearest5_2, toString(temp5_2))
}

cell_neighbors <- cbind(cell_neighbors, nearest5_1, nearest5_2)


list_of_df <- list()
for (i in 1:nrow(cell_neighbors)) {
  temp_list <- as.list(strsplit(cell_neighbors$nearest5_1, ", ")[[i]])
  
  temp_df <- assayData_filtered[,colnames(assayData_filtered) %in% temp_list]
  df <- as.data.frame(apply(as.matrix(temp_df),1,mean))
  colnames(df)<-cell_neighbors$cell[i]
  list_of_df <- append(list_of_df, df)
}
averages_1 <- bind_rows(list_of_df, .id = "column_label")
rownames(averages_1) <- rownames(assayData_filtered)

list_of_df <- list()
for (i in 1:nrow(cell_neighbors)) {
  temp_list <- as.list(strsplit(cell_neighbors$nearest5_2, ", ")[[i]])
  
  temp_df <- assayData_filtered[,colnames(assayData_filtered) %in% temp_list]
  df <- as.data.frame(apply(as.matrix(temp_df),1,mean))
  colnames(df)<-cell_neighbors$cell[i]
  list_of_df <- append(list_of_df, df)
}
averages_2 <- bind_rows(list_of_df, .id = "column_label")
rownames(averages_2) <- rownames(assayData_filtered)

seedf<-as.data.frame(cbind(rownames(assayData_filtered),assayData_filtered))
assayData_filtered[]

#calculating log2fc
averages_1_t <- as.data.frame(t(averages_1))
averages_2_t <- as.data.frame(t(averages_2))

Jun0min <- as.data.frame(averages_1_t$Jun)
Jun0min <- Jun0min+0.00000001
Jun50min <- as.data.frame(averages_2_t$Jun)
Jun50min <- Jun50min+0.00000001
target <- list()
for (i in 1:nrow(Jun50min)) {
  target <- append(target, log2(Jun50min[i,])-log2(Jun0min[i,]))
}

pbmc_filtered <- FindVariableFeatures(pbmc_filtered, selection.method = "vst", nfeatures = 5000)
variable_list = as.list(pbmc_filtered@assays$RNA@var.features)

counts5000 = subset(averages_1_t, select = names(averages_1_t) %in% variable_list)

#----------------
targetdf <- as.data.frame(target)
target_exp <- as.data.frame(t(targetdf))
rownames(target_exp) <- NULL



write.csv(target_exp, "C:/Users/okadalab/Documents/Mouse_scRNA-seq_data/aggr_0_50/target_log2fc.csv")
write.csv(counts5000, "C:/Users/okadalab/Documents/Mouse_scRNA-seq_data/aggr_0_50/counts_log2fc.csv")
