## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

## Ò»¡¢Load datasets and create Seurat objects with the raw (non-normalized data).
a <-  CreateSeuratObject(counts = Read10X(data.dir = "/public/home/yangjw28/data/BBB/GDR20070021-3sample10X-RNA/GDR20070021_order_1/Young/Young_count/outs/filtered_feature_bc_matrix"), project="a_Young", assay = "RNA")
b <-  CreateSeuratObject(counts = Read10X(data.dir = "/public/home/yangjw28/data/BBB/GDR20070021-3sample10X-RNA/GDR20070021_order_1/Middle/Middle_count/outs/filtered_feature_bc_matrix"), project="b_Middle",assay = "RNA")
c <-  CreateSeuratObject(counts = Read10X(data.dir = "/public/home/yangjw28/data/BBB/GDR20070021-3sample10X-RNA/GDR20070021_order_1/Old/Old__count/outs/filtered_feature_bc_matrix"), project="c_Old",assay = "RNA") 
#¼Ó3¸öprojectÇø·ÖÈý¸öÑùÆ·ÄÚµÄÏ¸°ûident
objList <- list(a,b,c)

##¶þ¡¢QC
#QC and filter out low-quality cells--for a
a[["percent.mt"]] <- PercentageFeatureSet(object =a, pattern = "^mt-")
#×¢Òâ²Î¿¼»ùÒò×éÀïÃæÏßÁ£ÌåÏà¹Ø»ùÒòµÄÃû³Æ
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/Young/Young_qc_plot.pdf")
p1=VlnPlot(a,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2=ggplot(data=a[["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
p3=ggplot(data=a[["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
p4=ggplot(data=a[["percent.mt"]],aes(x=percent.mt))+geom_density()
p5=FeatureScatter(a, feature1 = "nCount_RNA", feature2 = "percent.mt")
p6=FeatureScatter(a, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
dev.off()
# QC and filter out low-quality cells--for b
b[["percent.mt"]] <- PercentageFeatureSet(object =b, pattern = "^mt-")
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/Middle/Middle_qc_plot.pdf")
p1=VlnPlot(b,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2=ggplot(data=b[["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
p3=ggplot(data=b[["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
p4=ggplot(data=b[["percent.mt"]],aes(x=percent.mt))+geom_density()
p5=FeatureScatter(b, feature1 = "nCount_RNA", feature2 = "percent.mt")
p6=FeatureScatter(b, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
dev.off()
# QC and filter out low-quality cells--for c
c[["percent.mt"]] <- PercentageFeatureSet(object =c, pattern = "^mt-")
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/Old/Old_qc_plot.pdf")
p1=VlnPlot(c,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2=ggplot(data=c[["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
p3=ggplot(data=c[["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
p4=ggplot(data=c[["percent.mt"]],aes(x=percent.mt))+geom_density()
p5=FeatureScatter(c, feature1 = "nCount_RNA", feature2 = "percent.mt")
p6=FeatureScatter(c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
dev.off()

# Retain cells that express not more than 4000 genes (remove potential homotypic doublets ?) and had mitochondrial content <10%
a = subset(a,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )
b = subset(b,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )
c = subset(c,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )

# Retain cells that are both considered as real cells by cellranger and DropletUtils 
# Filter potential heterotypic doublets detected in silico modeling with 'scrublet'
#²»»á
#cell_level_data <- read.table("/public/home/shipy3/chinmo_project/input/LL3_testis_scRNA-seq/cell_level_data.txt",header=TRUE,sep="\t")

#for ( i in seq(1,3)){
#  retained_cells <- filter(cell_level_data,rep==str_c("rep",i),is_cell_used_in_study==TRUE)
#  retained_cells_barcodes <- str_c(gsub(pattern="rep\\d\\_","",x=retained_cells$cell_id),"-1")
#  objList[[i]] <- subset(objList[[i]],cells=retained_cells_barcodes)
#}


##Èý¡¢ÕûºÏÊý¾Ý¼¯

# Simply merge Seurat objects
merged_obj <- merge(x=a,y=c(b,c),add.cell.ids = paste("replicate",1:3,sep=""),project = "BBB_EC")
#¿ÉÒÔÓÃxxx[["yyy"]] <- zzz¸øxxxµÄµ¥Ï¸°ûÊý¾Ý¸³³ÆÎªyyyÊôÐÔ£¬ÖµÎªzzz
#Idents(merged_obj) <- gsub("_.*", "", colnames(merged_obj))

#split the combined object into a list, with each dataset as an element
BBB_EC.list <- SplitObject(merged_obj,split.by = "ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/YMO_highly_plot.pdf")
BBB_EC.list[[1]] <- NormalizeData(BBB_EC.list[[1]],verbose = FALSE)
BBB_EC.list[[1]] <- FindVariableFeatures(BBB_EC.list[[1]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(BBB_EC.list[[1]]),10)
plot1 <- VariableFeaturePlot(BBB_EC.list[[1]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
BBB_EC.list[[2]] <- NormalizeData(BBB_EC.list[[2]],verbose = FALSE)
BBB_EC.list[[2]] <- FindVariableFeatures(BBB_EC.list[[2]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(BBB_EC.list[[2]]),10)
plot3 <- VariableFeaturePlot(BBB_EC.list[[2]])
plot4 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot3
plot4
BBB_EC.list[[3]] <- NormalizeData(BBB_EC.list[[3]],verbose = FALSE)
BBB_EC.list[[3]] <- FindVariableFeatures(BBB_EC.list[[3]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(BBB_EC.list[[3]]),10)
plot5 <- VariableFeaturePlot(BBB_EC.list[[3]])
plot6 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot5
plot6
dev.off()

# identify anchors using the FindIntegrationAnchors function# Ê¶±ðÃªµã,Ê¹ÓÃFindIntegrationAnchorsº¯ÊýÀ´Ê¶±ðÃªµã£¨anchors£©£¬¸Ãº¯ÊýµÄÊäÈëÊý¾ÝÊÇSeurat¶ÔÏóµÄÁÐ±í
BBB_EC.anchors <- FindIntegrationAnchors(object.list = BBB_EC.list,anchor.features = 2000,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ¡®batch-corrected¡¯) expression matrix for all cells, enabling them to be jointly analyzed.
#½øÐÐÊý¾Ý¼¯ÕûºÏ
# ÒÑ¾­ÕûºÏºóµÄ±í´ï¾ØÕó´æ´¢ÔÚAssayÖÐ£¬Î´´¦ÀíµÄ±í´ï¾ÙÖ¤ÔÚRNA¶ÔÏóÖÐ
#È»ºó½«ÕâÐ©Ãª£¨anchors£©´«µÝ¸øIntegrateDataº¯Êý£¬¸Ãº¯Êý·µ»ØSeurat¶ÔÏó¡£
BBB_EC.integrated <- IntegrateData(anchorset = BBB_EC.anchors, dims = 1:30,features.to.integrate = rownames(BBB_EC.list[[1]]))
#head(BBB_EC.integrated@meta.data)
##ËÄ¡¢¿ÉÊÓ»¯
##switch to integrated assay
DefaultAssay(BBB_EC.integrated) <- "integrated"

# scale and center features in the dataset
BBB_EC.integrated <- ScaleData(BBB_EC.integrated, features =rownames(BBB_EC.integrated),verbose = FALSE)

# Perform linear dimensional reduction
BBB_EC.integrated <- RunPCA(BBB_EC.integrated, npcs = 50, verbose = FALSE)
#
#BBB_EC.integrated <- RunUMAP(BBB_EC.integrated, reduction = "pca", dims = 1:30)
#pdf("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/cell_cluster.pdf")
#p1 <- DimPlot(BBB_EC.integrated, reduction = "umap", group.by = "ident")
#p2 <- DimPlot(BBB_EC.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
#p1 
#p2
#dev.off()

# Determine the ¡®dimensionality¡¯ of the dataset
# JackStraw 
BBB_EC.integrated <- JackStraw(BBB_EC.integrated, num.replicate = 100, dims =50)
BBB_EC.integrated  <- ScoreJackStraw(BBB_EC.integrated, dims = 1:50)
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/integrated_pc.pdf")
JackStrawPlot(BBB_EC.integrated, dims = 1:50)
# ¡®Elbow plot¡¯
ElbowPlot(BBB_EC.integrated,ndims=50)
dev.off()

pdf("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/cell_cluster300-4000_15_20_0.4.pdf")
#¡ùÏ¸°û¾ÛÀà
BBB_EC.integrated <- FindNeighbors(object = BBB_EC.integrated, dims = 1:25)
BBB_EC.integrated <- FindClusters(object = BBB_EC.integrated, resolution = 0.6)
head(Idents(BBB_EC.integrated), 5)
table(BBB_EC.integrated@active.ident) # ²é¿´Ã¿Ò»ÀàÓÐ¶àÉÙ¸öÏ¸°û
# ÌáÈ¡Ä³Ò»ÀàÏ¸°û¡£
head(subset(as.data.frame(BBB_EC.integrated@active.ident),BBB_EC.integrated@active.ident=="2"))
#ÌáÈ¡Ä³Ò»clusterÏ¸°û
BBB_EC.integrated
subBBB_EC.integrated<-subset(x = BBB_EC.integrated,idents="0")
subBBB_EC.integrated
head(WhichCells(BBB_EC.integrated,idents="0"))
head(Idents(BBB_EC.integrated), 5)
#ÏµÍ³·¢ÓýÊ÷
#pbmc<-BuildClusterTree(BBB_EC.integrated)
#Tool(object = BBB_EC.integrated, slot = 'BuildClusterTree')
#PlotClusterTree(BBB_EC.integrated)

#Ö´ÐÐ·ÇÏßÐÔ½µÎ¬
BBB_EC.integrated <- RunUMAP(object = BBB_EC.integrated, dims = 1:25)
DimPlot(object = BBB_EC.integrated, group.by = "orig.ident",label = T)
DimPlot(object = BBB_EC.integrated, split.by = "orig.ident",label = T)
table(BBB_EC.integrated@meta.data$orig.ident,BBB_EC.integrated$seurat_clusters)
#²»Í¬ÄêÁäÏ¸°ûÔÚ²»Í¬clusterµÄ·Ö²¼
BBB_EC.integrated <- RunUMAP(object = BBB_EC.integrated, dims = 1:25)
DimPlot(object = BBB_EC.integrated, reduction = "umap",label = T)
BBB_EC.integrated <- RunTSNE(object = BBB_EC.integrated, dims = 1:25)
DimPlot(object = BBB_EC.integrated, reduction = "tsne",pt.size=0.1,label = T)
#group.by="seurat_clusters" group.by="celltype"

dev.off()
saveRDS(BBB_EC.integrated, file = "/public/home/yangjw28/projects/BBB/BBB_sc/YMO/BBB_EC.integrated350-4000_15_25_0.6_abc.rds")

#BBB_EC.integrated <- readRDS("BBB_EC.integrated300-4000_15_20_0.4_abc.rds")

#cluster annotation
#»ùÒòÕ¼±ÈÍ¼
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/cluster_annotation.pdf")
features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3","Tgfb2", "Glul","Slc26a10","Lypd1","Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc","Tfrc", "Car4","Itm2a","Chn2","Lcn2","Slc38a5","Nr2f2", "Sox12","Tbx1", "Tmbs10","Icam1","Vcam1", "Vwf", "P2ry1")
DotPlot(BBB_EC.integrated, features = features,dot.scale = 3) + RotatedAxis()
dev.off()

pdf("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/cluster_annotation2.pdf")
features <- c( "Plvap", "Plpp3","Esm1","Pdgfrb", "Cspg4","Kcnj8","Isg15", "Ifit1","Ifit3","Ifit3b","Prox1", "Pdpn","Acta2", "Pdlim3","Myh11","Aif1","Csf1r", "Cd68","Dcn", "Pdgfra","Mbp","Opalin", "Mobp","Syt1","Syp", "Eno2")
DotPlot(BBB_EC.integrated, features = features,dot.scale = 3) + RotatedAxis()
dev.off()

#Finding differentially expressed features (cluster biomarkers)

