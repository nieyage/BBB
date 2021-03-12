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
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/OnlyEC_selectPC.pdf")
JackStrawPlot(EC_cells.integrated, dims = 1:50)
ElbowPlot(EC_cells.integrated,ndims=50)
dev.off()
saveRDS(EC_cells.integrated,"YMO_onlyEC_integrated.rds")

EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/EC_cluster_tree.pdf")
#※细胞聚类
EC_cells.integrated <- FindNeighbors(object = EC_cells.integrated, dims = 1:20)
EC_cells.integrated <- FindClusters(object = EC_cells.integrated, resolution = 0.4)

#系统发育树
EC_cells.integrated<-BuildClusterTree(EC_cells.integrated)
Tool(object = EC_cells.integrated, slot = 'BuildClusterTree')
PlotClusterTree(EC_cells.integrated)
dev.off();

#执行非线性降维
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/EC_cluster_Umap.pdf")
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:22)
DimPlot(object = EC_cells.integrated, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated, split.by = "orig.ident",label = T)
table(EC_cells.integrated@meta.data$orig.ident,EC_cells.integrated$seurat_clusters)
dev.off();
             0    1    2    3    4    5    6    7    8    9   10
  Middle1 1703 1926  879 1184  858    7  462  893  412  178  107
  Old     3919 1691 1713 1155  776 1545  683  117  357  597  217
  Young   2238 1547 1049 1294  969   16  380  510  373  273  159


pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/EC_test_Umap.pdf")
DimPlot(object = EC_cells.integrated, reduction = "umap",label = T)
dev.off()

#####Only Young/Old/Middle12####
####free Umap##########

## 提取 Umap 降维信息到Umap

Umap = EC_cells.integrated@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx =EC_cells.integrated@meta.data$orig.ident)
pdf(file = "/md01/nieyg/project/BBB/YMO_results/ECplot/Umap_onlyYoung_OnlyOld.pdf", width = 8, height = 8)
p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c("Young" = "deepskyblue", "Old" = "grey"))   
p+labs(title="Young cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))

p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c( "Old"= "red", "Young" = "grey"))

p+labs(title="Aged cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
dev.off()

table(EC_cells.integrated@active.ident,EC_cells.integrated@meta.data$orig.ident) # 查看每一cluster有多少个细胞
table(EC_cells.integrated@meta.data$orig.ident,EC_cells.integrated$seurat_clusters) #查看3个时间点不同cluster有多少个细胞



EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")

#####把celltype改为cluster信息##############
new.cluster.ids <- factor(EC_cells.integrated$seurat_clusters)
head(new.cluster.ids)
EC_cells.integrated@active.ident <- new.cluster.ids
head(EC_cells.integrated@active.ident)


#cluster annotation
#基因占比图
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_onlyEC_cluster_annotation.pdf",width=12,height=8)

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


pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_onlyEC_anno_by_other_celltype.pdf")
features <- c( 	"Prox1", "Pdpn",#####Lymphatics
	"Acta2", "Pdlim3","Myh11",#####SMC
	"Aif1","Csf1r", "Cd68",#####Microglia
	"Dcn", "Pdgfra",###Fibroblast
	"Mbp",	"Opalin", "Mobp",#####OLi
	"Syt1","Syp", "Eno2")#######Neuron)
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()
#######重新用注释标记cluster####
new.cluster.ids <- c("Capillary_1","C_V_1","Interferon","C_V_2","C_A","Capillary_2","Arterial","Capillary_3",
  "Capillary_4","Venous","Choroid_plexus")
names(new.cluster.ids) <- levels(EC_cells.integrated)
EC_cells.integrated <- RenameIdents(EC_cells.integrated, new.cluster.ids)
EC_cells.integrated$celltype<-EC_cells.integrated@active.ident
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_EC_cells_annotation.pdf")
DimPlot(EC_cells.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "Brain ECs subtypes")
dev.off();
#saveRDS(EC_cells.integrated, file = "BBB_onlyEC.rds")

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_onlyEC_subtype_annotation.pdf",width=12,height=8)

features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
  "Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
  "Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",####1
  "Tfrc", "Car4", "Itm2a","Chn2",#####C-V
  "Lcn2","Slc38a5","Nr2f2", "Sox12","Icam1","Vcam1", "Vwf",##Venous
  "Plvap", "Plpp3","Esm1",####choroid plexus
  "Isg15", "Ifit1","Ifit3","Ifit3b"####Interferon
  )
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))

dev.off()

  
#Heatmap for BBB  genes(23) in different clusters 
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/23genes_heatmap_forsubtype.pdf")
for (i in levels(EC_cells.integrated)){
  subEC_cells.integrated<-subset(x = EC_cells.integrated,idents=i)
  features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
    'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
  p <-DoHeatmap(EC_cells.integrated, features = features, group.by = "orig.ident") 
  print(p)
}
dev.off()

DoHeatmap(EC_cells.integrated, features = features,size=3,angle = -50, hjust=0.8) 

dev.off()

markers <- FindMarkers(EC_cells.integrated, ident.1 = "Choroid_plexus", min.pct = 0.25,logfc.threshold=0)
head(markers, n=20)


pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/DEG_volcano_subtype_YO.pdf",width=5,height=5)
  mydata<-as.data.frame(markers)
  mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Up",  ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Down","normal"))
  ##对满足不同条件的数据给不同的标记，放入Condition列
  p <-ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+xlim(c(-2, 2)) +
              ggtitle("Choroid_plexus vs other ECs")+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Up'='red','Down'='deepskyblue','normal'='gray'));

dev.off()

#########DEGs GO and KEGG#################


pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus-up-GO.pdf",width=13,height=8)

gene.df <- bitr(up, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"UP-choroid_plexus-GO-BP.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"UP-choroid_plexus-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"UP-choroid_plexus-GO-CC.csv")


pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus-up-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"UP-choroid_plexus-KEGG.csv")


#####Only Young/Old/Middle12####
####free Umap##########

## 提取 Umap 降维信息到Umap
EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
DefaultAssay(EC_cells.integrated) <- "RNA"
Umap = EC_cells.integrated@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx =EC_cells.integrated@meta.data$orig.ident)
pdf(file = "/md01/nieyg/project/BBB/YMO_results/ECplot/Umap_YMA.pdf", width = 8, height = 8)
p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c("Young" = "#00A087", "Middle1"="grey","Old" = "grey"))   
p+labs(title="Young cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c("Young" = "grey", "Middle1"="#FFA040","Old" = "grey"))   
p+labs(title="Middle cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c( "Old"= "#E64B35", "Middle1"="grey","Young" = "grey"))
p+labs(title="Aged cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
dev.off()


