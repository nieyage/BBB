## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

##Load datasets and create Seurat objects with the raw (non-normalized data).
Old   <-  CreateSeuratObject(counts = Read10X(data.dir = "/md01/nieyg/ori/project/BBB/data/Old/Old__count/outs/filtered_feature_bc_matrix"), 
	project="Old", assay = "RNA")
Young <-  CreateSeuratObject(counts = Read10X(data.dir = "/md01/nieyg/ori/project/BBB/data/Young/Young_count/outs/filtered_feature_bc_matrix"), 
	project="Young", assay = "RNA")

/md01/yangjw28/projects/BBB/scRNA_Mid1/Mid1-count



#####添加sample信息####
objList <- list(Young,Old)

##############################
#######pre-processing#########
#######      QC      #########
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#注意参考基因组里面线粒体相关基因的名称
Old[["percent.mt"]] <- PercentageFeatureSet(Old, pattern = "^mt-")
Young[["percent.mt"]] <- PercentageFeatureSet(Young, pattern = "^mt-")
# Visualize QC metrics as a violin plot
pdf("/md01/nieyg/ori/project/BBB/plot/QC/Old_Young_qc_plot.pdf")
VlnPlot(Old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Young, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(Old, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Old, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot3 <- FeatureScatter(Young, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(Young, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4

# Retain cells that express not more than 4000 genes (remove potential homotypic doublets ?) and had mitochondrial content <15%
Young = subset(Young,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )
Old = subset(Old,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )

##########################################
###########Normalizing the data###########
##########################################


# Simply merge Seurat objects
merged_obj <- merge(x=Young,y=Old,add.cell.ids = c("Young","Old"),project = "BBB")
Idents(merged_obj) <- gsub("_.*", "", colnames(merged_obj))
#split the combined object into a list, with each dataset as an element
BBB_EC.list <- SplitObject(merged_obj,split.by = "ident")
pdf("/md01/nieyg/ori/project/BBB/plot/Young_normalization_plot.pdf")
BBB_EC.list[[1]] <- NormalizeData(BBB_EC.list[[1]],verbose = FALSE)
BBB_EC.list[[1]] <- FindVariableFeatures(BBB_EC.list[[1]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(BBB_EC.list[[1]]),10)
plot1 <- VariableFeaturePlot(BBB_EC.list[[1]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off();

pdf("/md01/nieyg/ori/project/BBB/plot/Old_normalization_plot.pdf")
BBB_EC.list[[2]] <- NormalizeData(BBB_EC.list[[2]],verbose = FALSE)
BBB_EC.list[[2]] <- FindVariableFeatures(BBB_EC.list[[2]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(BBB_EC.list[[2]]),10)
plot1 <- VariableFeaturePlot(BBB_EC.list[[2]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off();

###Scaling the data
###Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function:
###
###Shifts the expression of each gene, so that the mean expression across cells is 0
###Scales the expression of each gene, so that the variance across cells is 1
###This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
###The results of this are stored in pbmc[["RNA"]]@scale.data

# identify anchors using the FindIntegrationAnchors function# 识别锚点,使用FindIntegrationAnchors函数来识别锚点（anchors），该函数的输入数据是Seurat对象的列表
BBB_EC.anchors <- FindIntegrationAnchors(object.list = BBB_EC.list,anchor.features = 2000,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
#进行数据集整合
# 已经整合后的表达矩阵存储在Assay中，未处理的表达举证在RNA对象中
#然后将这些锚（anchors）传递给IntegrateData函数，该函数返回Seurat对象。
BBB_EC.integrated <- IntegrateData(anchorset = BBB_EC.anchors, dims = 1:30,features.to.integrate = rownames(BBB_EC.list[[1]]))
head(BBB_EC.integrated@meta.data)

##四、可视化
##switch to integrated assay
DefaultAssay(BBB_EC.integrated) <- "integrated"

# scale and center features in the dataset
BBB_EC.integrated <- ScaleData(BBB_EC.integrated, features =rownames(BBB_EC.integrated),verbose = FALSE)

# Perform linear dimensional reduction
BBB_EC.integrated <- RunPCA(BBB_EC.integrated, npcs = 50, verbose = FALSE)
#
BBB_EC.integrated <- RunUMAP(BBB_EC.integrated, reduction = "pca", dims = 1:30)
pdf("/md01/nieyg/ori/project/BBB/plot/Umap-sample.pdf.pdf")
p1 <- DimPlot(BBB_EC.integrated, reduction = "umap", group.by = "ident")
p2 <- DimPlot(BBB_EC.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
p1 
p2
dev.off()
# Determine the ‘dimensionality’ of the dataset
# JackStraw 
BBB_EC.integrated <- JackStraw(BBB_EC.integrated, num.replicate = 100, dims =50)
BBB_EC.integrated  <- ScoreJackStraw(BBB_EC.integrated, dims = 1:50)
pdf("/md01/nieyg/ori/project/BBB/plot/seclect_pca.pdf")
JackStrawPlot(BBB_EC.integrated, dims = 1:50)
# ‘Elbow plot’
ElbowPlot(BBB_EC.integrated,ndims=50)
dev.off()

pdf("/md01/nieyg/ori/project/BBB/plot/cell_cluster.pdf")
#※细胞聚类
BBB_EC.integrated <- FindNeighbors(object = BBB_EC.integrated, dims = 1:25)
BBB_EC.integrated <- FindClusters(object = BBB_EC.integrated, resolution = 0.6)
head(Idents(BBB_EC.integrated), 5)
table(BBB_EC.integrated@active.ident) # 查看每一类有多少个细胞
# 提取某一类细胞。
head(subset(as.data.frame(BBB_EC.integrated@active.ident),BBB_EC.integrated@active.ident=="2"))
#提取某一cluster细胞
BBB_EC.integrated
subBBB_EC.integrated<-subset(x = BBB_EC.integrated,idents="0")
subBBB_EC.integrated
head(WhichCells(BBB_EC.integrated,idents="0"))
head(Idents(BBB_EC.integrated), 5)
#系统发育树

pbmc<-BuildClusterTree(BBB_EC.integrated)
Tool(object = BBB_EC.integrated, slot = 'BuildClusterTree')
PlotClusterTree(BBB_EC.integrated)

#执行非线性降维
BBB_EC.integrated <- RunUMAP(object = BBB_EC.integrated, dims = 1:25)
DimPlot(object = BBB_EC.integrated, group.by = "orig.ident",label = T)
DimPlot(object = BBB_EC.integrated, split.by = "orig.ident",label = T)
table(BBB_EC.integrated@meta.data$orig.ident,BBB_EC.integrated$seurat_clusters)
#不同年龄细胞在不同cluster的分布
BBB_EC.integrated <- RunUMAP(object = BBB_EC.integrated, dims = 1:25)
DimPlot(object = BBB_EC.integrated, reduction = "umap",label = T)
#BBB_EC.integrated <- RunTSNE(object = BBB_EC.integrated, dims = 1:25)
#DimPlot(object = BBB_EC.integrated, reduction = "tsne",pt.size=0.1,label = T)
#group.by="seurat_clusters" group.by="celltype"

dev.off()
saveRDS(BBB_EC.integrated, file = "BBB_YO.rds")

#BBB_EC.integrated <- readRDS("BBB_EC.integrated300-4000_15_20_0.4_abc.rds")

#cluster annotation
#基因占比图
pdf("/md01/nieyg/project/BBB/plot/cluster_annotation-all.pdf",width=16,height=8)

features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
	"Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
	"Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",####1
	"Tfrc", "Car4",	"Itm2a","Chn2",#####C-V
	"Lcn2","Slc38a5","Nr2f2", "Sox12",##Venous
	"Tbx1", "Tmbs10","Icam1","Vcam1", "Vwf", "P2ry1", ####A/V
	"Plvap", "Plpp3","Esm1",####choroid plexus
	"Pdgfrb", "Cspg4","Kcnj8",####Pericyte
	"Isg15", "Ifit1","Ifit3","Ifit3b",####Interferon
	"Prox1", "Pdpn",#####Lymphatics
	"Acta2", "Pdlim3","Myh11",#####SMC
	"Aif1","Csf1r", "Cd68",#####Microglia
	"Dcn", "Pdgfra",###Fibroblast
	"Mbp",	"Opalin", "Mobp",#####OLi
	"Syt1","Syp", "Eno2")#######Neuron
DotPlot(BBB_EC.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()
