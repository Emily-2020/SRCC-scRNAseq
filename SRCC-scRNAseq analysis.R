#GSE241934分析
setwd("D:/wangsen/HLL SRCC/GSE241934/")
GSE241934=Read10X(data.dir = ("D:/wangsen/HLL SRCC/GSE241934/IIT"))
GSE241934_2=Read10X(data.dir = ("D:/wangsen/HLL SRCC/GSE241934/RWC"))
GSE241934 <- cbind(GSE241934, GSE241934_2)
GSE241934=CreateSeuratObject(counts = GSE241934,project = "GSE241934", min.cells = 3, min.features = 200)
rm("GSE241934_2")
gc()
GSE241934 <- NormalizeData(object = GSE241934)
GSE241934 <- FindVariableFeatures(object = GSE241934)
GSE241934 <- ScaleData(object = GSE241934)
GSE241934 <- RunPCA(object = GSE241934,npcs = 100)
GSE241934 <- FindNeighbors(object = GSE241934, dims = 1:30)
GSE241934 <- FindClusters(object = GSE241934,resolution = 0.05)
GSE241934 <- RunUMAP(object = GSE241934, dims = 1:30)
GSE241934 <- RunTSNE(object = GSE241934, dims = 1:30)
DimPlot(object = GSE241934)
Celltype = c("CD3D", "CD3E","CD79A", "IGHM","DCN", "COL1A1","CD68", "MARCO","PECAM1", "CLDN5", "EPCAM", "CDH1","KIT", "MS4A2"
)
DotPlot(GSE241934,features =Celltype,cols = c("grey","tomato"))+
  scale_size (range=c (0.1, 12))
new.cluster.ids <- c("T cell", "B cell","Fibrolast","Myeloid cell", "Endothelial", "Epithelial",
                     "B cell", "Epithelial", "Mast cell","Epithelial","Epithelial","Myeloid cell")
names(new.cluster.ids) <- levels(GSE241934)
GSE241934 <- RenameIdents(GSE241934, new.cluster.ids)
DimPlot(GSE241934, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DotPlot(GSE241934,features =Celltype,cols = c("grey","tomato"))+
  scale_size (range=c (0.1, 12))
##进行cnv寻找恶性上皮
B=subset(GSE241934,idents = c(1,2,5,6,7,9,10))
B <- RunPCA(object = B,npcs = 50)
B <- FindNeighbors(object = B, dims = 1:20)
B <- FindClusters(object = B,resolution = 0.5)
B <- RunUMAP(object = B, dims = 1:20)
DimPlot(B,label = T)
###逐步开始infercnv操作
a=GetAssayData(B,layer = "counts")
a=as(a, "dgCMatrix")
anno=data.frame(Idents(B))
gene_order= read.table("D:/wangsen/infercnv练习文件夹/gene_order_file.txt", header = TRUE, sep = "\t",row.names = 1)
row.names(gene_order)
a=as.data.frame(a)
a=a[rownames(a)%in%rownames(gene_order),]
B_infercnvobj= CreateInfercnvObject(raw_counts_matrix = a,
                                    annotations_file = anno,
                                    delim = "\t",
                                    gene_order_file= gene_order,
                                    ref_group_names = NULL)
B_infercnvobj= infercnv::run(B_infercnvobj,
                             cutoff = 0.1,
                             out_dir = "D:/wangsen/HLL SRCC/GSE241934/infercnv",
                             cluster_by_groups = T,
                             denoise = T,
                             write_expr_matrix = T,
                             HMM=F)
cnvO=read.table("D:/wangsen/HLL SRCC/GSE241934/infercnv/infercnv.observations.txt",header = TRUE)
cnvO <- t(cnvO)
cnvO <- as.data.frame(cnvO)
cnvO$cell <- rownames(cnvO)
cnvO$CNV_Score <- rowMeans(cnvO[, -which(names(cnvO) == "cell")])
result_matrix <- cnvO[, c("cell", "CNV_Score")]#得到每个细胞的CNV分数
rownames(result_matrix) <- gsub("\\.1", "-1", rownames(result_matrix))
colnames(result_matrix) <- c("cell", "CNV_Score")
result_matrix <- result_matrix[, "CNV_Score", drop = FALSE]
anno$cell <- rownames(anno)
result_matrix$cell <- rownames(result_matrix)
merged_data <- merge(anno, result_matrix, by = "cell", all.x = TRUE)#合并anno文件和CNV结果
cnv_scores <- merged_data[, c("cell", "CNV_Score")]
cell_ids <- rownames(B@meta.data)
cnv_scores <- cnv_scores[cnv_scores$cell %in% cell_ids,]# 将CNV_Score添加到Seurat对象B的meta data中，确保细胞ID匹配
B@meta.data$CNV_Score <- NA  # 初始化CNV_Score列
B@meta.data[cnv_scores$cell, "CNV_Score"] <- cnv_scores$CNV_Score#为每个细胞加上CNV分数
VlnPlot(B, features = "CNV_Score", group.by = "seurat_clusters",raster = F) +ggtitle("CNV Scores by Cluster") +
  theme_minimal()# 绘制小提琴图
meta_data <- B@meta.data %>% 
  select(CNV_Score, seurat_clusters)#导出CNV
meta_data$CNV_Score <- as.numeric(meta_data$CNV_Score)# 确保CNV_Score是数值型
cluster_means <- meta_data %>%
  group_by(seurat_clusters) %>%
  summarize(mean_CNV_Score = mean(CNV_Score, na.rm = TRUE))# 计算每个Cluster的平均CNV分数
DotPlot(B,features =c("EPCAM","KRT18","KRT8","DCN","CD79A") )
tumor=subset(B,idents = c(6,23,20,15,17,12,8))
saveRDS(tumor,"GSE241934_tumor.Rds")

#对tumor进行去双胞
sweep.res.list= paramSweep(tumor,PCs = 1:13,sct = F)
sweep.stats<-summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn<- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%as.character()%>% as.numeric()
DoubletRate = ncol(tumor)*8*1e-7 
DoubletRate = 0.074
homotypic.prop<- modelHomotypic(tumor$seurat_clusters)
nExp_poi<- round(DoubletRate*ncol(tumor))
nExp_poi.adj<- round(nExp_poi*(1-homotypic.prop))
ncol(tumor)
tumor <- doubletFinder(tumor,PCs = 1:13,pN = 0.25,pK = pK_bcmvn,nExp = nExp_poi.adj )
DimPlot(tumor,group.by = "DF.classifications_0.25_0.24_883")
tumor_nonDoubletFind=tumor
tumor=subset(tumor,subset = DF.classifications_0.25_0.24_883 == "Singlet")
#对得到的单个细胞的tumor进行分析
tumor <- NormalizeData(object = tumor)
tumor <- FindVariableFeatures(object = tumor,nfeatures = 2000)
marker_genes <- c("SLC52A2","NOTCH1")
variable_genes <- unique(c(VariableFeatures(tumor), marker_genes))
VariableFeatures(tumor) <- variable_genes
tumor <- ScaleData(object = tumor)
expr_matrix <- GetAssayData(tumor, slot = "data")
expr_matrix_sparse <- as(expr_matrix, "dgCMatrix")
library(Matrix)
library(RSpectra)
library(Matrix)
library(irlba)
weighted_pca_sparse <- function(data, weights, variableFeature, n_components) {
  # 提取variableFeature中的基因数据
  data_filtered <- data[variableFeature, , drop = FALSE]
  # 确保权重向量的长度与过滤后数据矩阵的行数相同
  if (length(weights) != nrow(data_filtered)) {
    stop("Length of weights vector must match the number of rows in the filtered data matrix.")
  }
  # 将权重向量变为对角稀疏矩阵
  weight_matrix <- Diagonal(x = weights)
  # 执行矩阵乘法，以保持稀疏性
  weighted_data <- weight_matrix %*% data_filtered
  # 计算中心化矩阵
  centered_weighted_data <- weighted_data - Matrix::rowMeans(weighted_data)
  # 执行稀疏SVD，k 为主成分数量
  svd_result <- svds(centered_weighted_data, k = n_components)
  # 提取SVD结果
  pca_result <- list()
  pca_result$sdev <- svd_result$d / sqrt(nrow(data_filtered) - 1)
  pca_result$rotation <- svd_result$v
  pca_result$x <- svd_result$u %*% diag(svd_result$d)
  # 设置行名和列名
  rownames(pca_result$rotation) <- colnames(data_filtered)
  colnames(pca_result$rotation) <- paste0("PC_", 1:ncol(pca_result$rotation))
  rownames(pca_result$x) <- rownames(data_filtered)
  colnames(pca_result$x) <- paste0("PC_", 1:ncol(pca_result$x))
  class(pca_result) <- "prcomp"
  return(pca_result)
}
weights <- rep(1, length(variable_genes))#设置权重向量
weights[which(variable_genes %in% marker_genes)] <- 10# 设置特定基因的权重为10
pca_result <- weighted_pca_sparse(data = expr_matrix_sparse, weights = weights,
                                  variableFeature = variable_genes,
                                  n_components = 15 )
pca_dimred <- CreateDimReducObject(
  embeddings = pca_result$rotation,
  loadings = pca_result$x,
  stdev = pca_result$sdev,
  assay = "RNA",
  key = "PC_"
)
tumor[["pca"]] <- pca_dimred
tumor <- FindNeighbors(object = tumor, dims = 1:10)
tumor <- FindClusters(object = tumor,resolution = 0.2)
tumor <- RunUMAP(object = tumor, dims = 1:10,)
DimPlot(object = tumor,label = T)
DotPlot(object = tumor,features = c("SLC52A2","NOTCH1","SLC52A1","SLC52A3"))
FeaturePlot(object = tumor,features = c("SLC52A2","NOTCH1"))
plot=FeaturePlot(object = tumor,features = c("SLC52A2","NOTCH1"))
cells=CellSelector(plot = plot)##选择SLC52A2+NOTCH1+细胞
Idents(tumor,cells=cells) <- "SRCC"
table(Idents(tumor))
cell=WhichCells(tumor,idents = 11)
Idents(tumor,cells=cell) <- "NSRCC"
DimPlot(tumor,label = T)
###characteristic analysis
CellCycling_Feature=list(Cycling=c("CDKN1B","CCNE1","CCND1","CDK6","CDK4","CDK2"),
                         G0arrest=c("PSAP","SERINC1","YPEL5","YPEL3","EPS8L2","TXNIP") )
NOTCH <- c("NRARP","MYC","HEYL","HEY1","HES5","HES4","HES1","NOTCH1")
antioxidants = c("MDH1","IDH1","GPX4","GSTP1","PRDX6","PRDX1")
EMT= c("VIM","CDH2","ZEB2","ZEB1","TWIST1","SNAI2","SNAI1","CDH1","EPCAM")

DotPlot(tumor,features =CellCycling_Feature,idents = c("SRCC","NSRCC"),cols = c("grey","tomato"))+
  scale_size (range=c (2, 12))+
  scale_y_discrete(,expand = c(0.1,1))+
  theme(plot.margin = unit(c(0.1,10, 0.1, 0.1), "cm"))+
  coord_flip()
DotPlot(tumor,features =antioxidants,idents = c("SRCC","NSRCC"),cols = c("grey","tomato"))+
  scale_size (range=c (2, 12))+
  scale_y_discrete(,expand = c(0.1,1))+
  theme(plot.margin = unit(c(0.1,10, 0.1, 0.1), "cm"))+
  coord_flip()
DotPlot(tumor,features =NOTCH,idents = c("SRCC","NSRCC"),cols = c("grey","tomato"))+
  scale_size (range=c (2, 12))+
  scale_y_discrete(,expand = c(0.1,1))+
  theme(plot.margin = unit(c(0.1,10, 0.1, 0.1), "cm"))+
  coord_flip()
DotPlot(tumor,features =EMT,idents = c("SRCC","NSRCC"),cols = c("grey","tomato"))+
  scale_size (range=c (2, 12))+
  scale_y_discrete(,expand = c(0.1,1))+
  theme(plot.margin = unit(c(0.1,10, 0.1, 0.1), "cm"))+
  coord_flip() 
