## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
setwd("/public/home/yangjw28/projects/BBB/BBB_sc")
BBB_EC.integrated <- readRDS("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/BBB_EC.integrated350-4000_15_25_0.5_abc.rds")
#subset 取BBB亚群0-12
#subBBB_EC.integrated <-BBB_EC.integrated[,which(BBB_EC.integrated@meta.data$seurat_clusters==c("12","11"))]

DefaultAssay(BBB_EC.integrated) <- "RNA" # Create dummy new assay to demo switching default assays
EC_cells <- subset(BBB_EC.integrated,idents=c(0:10))
EC_cells <- DietSeurat(
  object=EC_cells,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
)


## split the integrated object into a list, with each dataset as an element
EC_cells.ls <- SplitObject(EC_cells,split.by = "orig.ident")


## identify variable features data里面的是已经标准化的，只需要找变异特征
for (i in 1:length(EC_cells.ls)){
  EC_cells.ls[[i]] <- FindVariableFeatures(EC_cells.ls[[i]],selection.method = "vst", nfeatures = 3000)
}


## identify anchors using the FindIntegrationAnchors function
EC_cells.anchors <- FindIntegrationAnchors(object.list = EC_cells.ls,anchor.features = 3000,dims = 1:30)


## pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
EC_cells.integrated <- IntegrateData(anchorset = EC_cells.anchors, dims = 1:30,features.to.integrate = rownames(EC_cells.ls[[1]]))


## switch to integrated assay
DefaultAssay(EC_cells.integrated) <- "integrated"


## scale and center features in the dataset
EC_cells.integrated <- ScaleData(EC_cells.integrated, features =rownames(EC_cells.integrated))

## Perform linear dimensional reduction
EC_cells.integrated <- RunPCA(EC_cells.integrated, npcs = 50, verbose = FALSE)


## Determine the ‘dimensionality’ of the dataset
EC_cells.integrated <- JackStraw(EC_cells.integrated, num.replicate = 100, dims =50)
EC_cells.integrated <- ScoreJackStraw(EC_cells.integrated, dims = 1:50)
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/2integrated_EC_cells_pc.pdf")
JackStrawPlot(EC_cells.integrated, dims = 1:50)
ElbowPlot(EC_cells.integrated,ndims=50)
dev.off()

pdf("/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/EC_2cluster_350-4000_15_25_0.4.pdf")
#※细胞聚类
EC_cells.integrated <- FindNeighbors(object = EC_cells.integrated, dims = 1:25)
EC_cells.integrated <- FindClusters(object = EC_cells.integrated, resolution = 0.4)

#系统发育树
EC_cells.integrated<-BuildClusterTree(EC_cells.integrated)
Tool(object = EC_cells.integrated, slot = 'BuildClusterTree')
PlotClusterTree(EC_cells.integrated)

#执行非线性降维
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:25)
DimPlot(object = EC_cells.integrated, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated, split.by = "orig.ident",label = T)
table(EC_cells.integrated@meta.data$orig.ident,EC_cells.integrated$seurat_clusters)
#不同年龄细胞在不同cluster的分布
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:25)
DimPlot(object = EC_cells.integrated, reduction = "umap",label = T)
EC_cells.integrated <- RunTSNE(object = EC_cells.integrated, dims = 1:25)
DimPlot(object = EC_cells.integrated, reduction = "tsne",pt.size=0.1,label = T)
#group.by="seurat_clusters" group.by="celltype"

dev.off()
saveRDS(EC_cells.integrated, file = "/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/EC_cell_2.integrated350-4000_15_25_0.4_abc.rds")
#EC_cell.integrated <- readRDS("EC_cell.integrated350-4000_15_25_0.5_abc.rds")
#分clustrer正常？
table(EC_cells.integrated@active.ident) # 查看每一cluster有多少个细胞
table(EC_cells.integrated@meta.data$orig.ident,EC_cell.integrated$seurat_clusters) #查看3个样品不同cluster有多少个细胞


#cluster annotation
#基因占比图
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/cluster_2annotation1.pdf")
features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3","Tgfb2", "Glul","Slc26a10","Lypd1","Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc","Tfrc", "Car4","Itm2a","Chn2","Lcn2","Slc38a5","Nr2f2", "Sox12","Tbx1", "Tmbs10","Icam1","Vcam1", "Vwf", "P2ry1")
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) + RotatedAxis()+ theme(axis.text.x=element_text(size=8))
dev.off()

pdf("/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/cluster_2annotation2.pdf")
features <- c( "Plvap", "Plpp3","Esm1","Pdgfrb", "Cspg4","Kcnj8","Isg15", "Ifit1","Ifit3","Ifit3b","Prox1", "Pdpn","Acta2", "Pdlim3","Myh11","Aif1","Csf1r", "Cd68","Dcn", "Pdgfra","Mbp","Opalin", "Mobp","Syt1","Syp", "Eno2")
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) + RotatedAxis() + theme(axis.text.x=element_text(size=8))
dev.off()

#3个年龄样本细胞数目占比的barplot
library("ggsci")
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/barplot_ECcell_2cluster.pdf")
ggplot(EC_cells.integrated@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar(position="fill") + scale_fill_igv()
#+coord_flip()
#ggplot(BBB_EC.integrated@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar()
dev.off()


#23个基因在不同年龄的表达量热图
pdf("/public/home/yangjw28/projects/BBB/BBB_sc/EC_cell/23_heatmap_2cluster.pdf")
for (i in seq(0,17)){
  subBBB_EC.integrated<-subset(x = EC_cells.integrated,idents=i)
  features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2",'Nt5c2','Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2','Palmd')
  p <-DoHeatmap(subBBB_EC.integrated, features = features, group.by = "orig.ident") 
  print(p)
}
dev.off()















