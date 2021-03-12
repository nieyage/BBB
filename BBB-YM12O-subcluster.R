## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
#subset 取BBB亚群0-9,12
#subBBB_EC <-BBB_all[,which(BBB_all@meta.data$seurat_clusters==c(0:9,12))]
BBB_EC.integrated<-readRDS("BBB_YM12O.rds")
DefaultAssay(BBB_EC.integrated) <- "RNA" # Create dummy new assay to demo switching default assays
EC_cells <- subset(BBB_EC.integrated,idents=c("C_V_1","C_A_1","Capillary_1","C_V_2","C_A_2","Venous","Interferon","Arterial_1","Capillary_2",
  "Arterial_2","Capillary_3","Choroid_plexus"))

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
pdf("/md01/nieyg/project/BBB/YM12O_results/plot/OnlyEC_selectPC.pdf")
JackStrawPlot(EC_cells.integrated, dims = 1:50)
ElbowPlot(EC_cells.integrated,ndims=50)
dev.off()
saveRDS(EC_cells.integrated,"YM12O_onlyEC_integrated.rds")

EC_cells.integrated<-readRDS("YM")
pdf("/md01/nieyg/project/BBB/YM12O_results/plot/EC_cluster_tree.pdf")
#※细胞聚类
EC_cells.integrated <- FindNeighbors(object = EC_cells.integrated, dims = 1:20)
EC_cells.integrated <- FindClusters(object = EC_cells.integrated, resolution = 0.4)

#系统发育树
EC_cells.integrated<-BuildClusterTree(EC_cells.integrated)
Tool(object = EC_cells.integrated, slot = 'BuildClusterTree')
PlotClusterTree(EC_cells.integrated)

#执行非线性降维
pdf("/md01/nieyg/project/BBB/YM12O_results/onlyEC/plot/EC_cluster_Umap.pdf")
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:22)
DimPlot(object = EC_cells.integrated, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated, split.by = "orig.ident",label = T)
table(EC_cells.integrated@meta.data$orig.ident,EC_cells.integrated$seurat_clusters)

             0    1    2    3    4    5    6    7    8    9   10   11
  Middle1 1255 1350 1677 1017  849  604  456 1133  417  332  107    3
  Middle2 1515 1487  953  432  445  561  653  493  392  239  151   14
  Old     3526 1466 1640 1361 1206  968  834  447  845  362  214   47
  Young   1961 1407 1362  821  526  657  645  442  507  358  157    9


pdf("/md01/nieyg/project/BBB/YM12O_results/onlyEC/plot/EC_test_Umap.pdf")
DimPlot(object = EC_cells.integrated, reduction = "umap",label = T)
dev.off()

#####Only Young/Old/Middle12####
####free Umap##########

## 提取 Umap 降维信息到Umap

Umap = EC_cells.integrated@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx =EC_cells.integrated@meta.data$orig.ident)
pdf(file = "/md01/nieyg/project/BBB/plot/Umap_onlyYoung_OnlyOld.pdf", width = 8, height = 8)
p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c("Young" = "deepskyblue", "Old" = "grey"))   
p+labs(title="Young cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))

p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c( "Old"= "red", "Young" = "grey"))

p+labs(title="Old cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
dev.off()

table(EC_cells.integrated@active.ident,EC_cells.integrated@meta.data$orig.ident) # 查看每一cluster有多少个细胞
table(EC_cells.integrated@meta.data$orig.ident,EC_cells.integrated$seurat_clusters) #查看3个时间点不同cluster有多少个细胞

#####把celltype改为cluster信息##############
new.cluster.ids <- factor(EC_cells.integrated$seurat_clusters)
head(new.cluster.ids)
EC_cells.integrated@active.ident <- new.cluster.ids
head(EC_cells.integrated@active.ident)


#cluster annotation
#基因占比图
pdf("/md01/nieyg/project/BBB/plot/BBB_onlyEC_cluster_annotation.pdf",width=12,height=8)

features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
	"Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
	"Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",####1
	"Tfrc", "Car4",	"Itm2a","Chn2",#####C-V
	"Lcn2","Slc38a5","Nr2f2", "Sox12",##Venous
	"Tbx1", "Tmbs10","Icam1","Vcam1", "Vwf", "P2ry1", ####A/V
	"Plvap", "Plpp3","Esm1",####choroid plexus
	"Pdgfrb", "Cspg4","Kcnj8",####Pericyte
	"Isg15", "Ifit1","Ifit3","Ifit3b"####Interferon
	)
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))

dev.off()


pdf("/md01/nieyg/project/BBB/plot/BBB_onlyEC_anno_by_other_celltype.pdf")
features <- c( 	"Prox1", "Pdpn",#####Lymphatics
	"Acta2", "Pdlim3","Myh11",#####SMC
	"Aif1","Csf1r", "Cd68",#####Microglia
	"Dcn", "Pdgfra",###Fibroblast
	"Mbp",	"Opalin", "Mobp",#####OLi
	"Syt1","Syp", "Eno2")#######Neuron)
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()
#######重新用注释标记cluster####
new.cluster.ids <- c("Capillary","C-V","C-V","Capillary","Arterial","Venous",
	    "Capillary","C-A","choroid plexus","Interferon")
names(new.cluster.ids) <- levels(EC_cells.integrated)
EC_cells.integrated <- RenameIdents(EC_cells.integrated, new.cluster.ids)
pdf("/md01/nieyg/project/BBB/plot/BBB_EC_cells_annotation.pdf")
DimPlot(EC_cells.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "Brain ECs subtypes")
dev.off();
saveRDS(EC_cells.integrated, file = "BBB_onlyEC.rds")

