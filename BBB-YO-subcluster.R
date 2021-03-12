## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
BBB_all<-readRDS("BBB_YO.rds")
#########寻找DEG########
######only up genes#####
cluster12.markers <- FindMarkers(BBB_all, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n=20)
cluster12.markers<-cluster12.markers[order(cluster12.markers$avg_logFC,decreasing=T),]

new.cluster.ids <- c("C_V",rep("Capillary",3),"Venous","C_A","Interferon","Arterial","Capillary",
	"Choroid_plexus","Microglia","Pericyte_SMC","Erythroid_cell","Oligo_Lymp")
names(new.cluster.ids) <- levels(BBB_all)
BBB_all <- RenameIdents(BBB_all, new.cluster.ids)
pdf("/md01/nieyg/project/BBB/plot/ALLcelltype_annotation.pdf")
DimPlot(BBB_all, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "Brain ECs with other celltypes")


#subset 取BBB亚群0-9,12
#subBBB_EC <-BBB_all[,which(BBB_all@meta.data$seurat_clusters==c(0:9,12))]
DefaultAssay(BBB_all) <- "RNA" # Create dummy new assay to demo switching default assays
EC_cells <- subset(BBB_all,idents=c("Venous","Capillary","Interferon","Arterial","Choroid plexus"))

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
pdf("/md01/nieyg/project/BBB/plot/OnlyEC_selectPC.pdf")
JackStrawPlot(EC_cells.integrated, dims = 1:50)
ElbowPlot(EC_cells.integrated,ndims=50)
dev.off()

pdf("/md01/nieyg/project/BBB/plot/EC_cluster_tree.pdf")
#※细胞聚类
EC_cells.integrated <- FindNeighbors(object = EC_cells.integrated, dims = 1:22)
EC_cells.integrated <- FindClusters(object = EC_cells.integrated, resolution = 0.4)

#系统发育树
EC_cells.integrated<-BuildClusterTree(EC_cells.integrated)
Tool(object = EC_cells.integrated, slot = 'BuildClusterTree')
PlotClusterTree(EC_cells.integrated)

#执行非线性降维
pdf("/md01/nieyg/project/BBB/plot/EC_cluster_Umap.pdf")
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:22)
DimPlot(object = EC_cells.integrated, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated, split.by = "orig.ident",label = T)
table(EC_cells.integrated@meta.data$orig.ident,EC_cells.integrated$seurat_clusters)
pdf("/md01/nieyg/project/BBB/plot/EC_test_Umap.pdf")
DimPlot(object = EC_cells.integrated, reduction = "umap",label = T) + labs(title = "Brain ECs")

dev.off()

#####Only Young/Old####
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
new.cluster.ids <- c("Capillary","C_V","C_A","Capillary","Arterial","Venous",
	    "Capillary","C_A","choroid_plexus","Interferon")
names(new.cluster.ids) <- levels(EC_cells.integrated)
EC_cells.integrated <- RenameIdents(EC_cells.integrated, new.cluster.ids)
pdf("/md01/nieyg/project/BBB/plot/BBB_EC_cells_annotation.pdf")
DimPlot(EC_cells.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "Brain ECs subtypes")
dev.off();
saveRDS(EC_cells.integrated, file = "BBB_onlyEC.rds")

EC_cells.integrated<-readRDS("BBB_onlyEC.rds")
#####寻找Cluster3、6的top marker ######
DefaultAssay(EC_cells.integrated)<-"RNA"
levels(EC_cells.integrated)

cluster3.markers <- FindMarkers(EC_cells.integrated, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n=20)
cluster3.markers<-cluster3.markers[order(cluster3.markers$avg_logFC,decreasing=T),]

cluster6.markers <- FindMarkers(EC_cells.integrated, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n=20)
cluster6.markers<-cluster6.markers[order(cluster6.markers$avg_logFC,decreasing=T),]

library(clusterProfiler)
library(org.Mm.eg.db)
pdf("cluster6-logFC05-pval005-GO.pdf")
gene.df <- bitr(gene, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Mm.eg.db)

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

#两个时间点样本不同cell type细胞数目占比的barplot
######金字塔barplot#######
library("ggsci")
pdf("/md01/nieyg/project/BBB/plot/ECcell_barplot_2ways.pdf")
EC_cells.integrated<-readRDS("BBB_onlyEC.rds")
m<-table(EC_cells.integrated@meta.data$seurat_clusters,EC_cells.integrated@meta.data$orig.ident)
     Old Young
  0 4699  3034
  1 1588  1867
  2 1520  1392
  3 2430    95
  4 1148   872
  5  766   306
  6   14   643
  7  302   311
  8  217   158
  9  128   142

library(DescTools)
library(ggsci)
mypal =pal_npg("nrc", alpha =0.7)(7)
mypal
meta<-matrix(data=NA,nrow=7,ncol=2);
colnames(meta)<-c("Old","Young");
rownames(meta)<-c("Capillary","C_V","Venous","C_A","Arterial","choroid_plexus","Interferon");
meta[,1]<-c(7134,1588,766,1822,1148,217,128)
meta[,2]<-c(3772,1867,306,1703,872,158,142)
PlotPyramid(meta[,1],meta[,2],
            ylab =rownames(meta),space=0,
            #col=c('cornflowerblue','indianred'),
            xlim=c(-8000,4000),
            main = 'EC cells distribution at Old/Young',
            lxlab='Old',rxlab='Young',
            col=rep(mypal,each=2),
            gapwidth = 0,ylab.x = -8000)


m<-EC_cells.integrated@meta.data;
table(EC_cells.integrated@active.ident) # 查看每一cluster有多少个细胞
table(EC_cells.integrated@meta.data$orig.ident,EC_cell.integrated$seurat_clusters) #查看3个时间点不同cluster有多少个细胞
dev.off()


#23个基因在不同年龄的表达量热图
pdf("/md01/nieyg/project/BBB/plot/23genes_heatmap_forYO.pdf")
for (i in levels(EC_cells.integrated)){
  EC_cells.integrated<-subset(x = EC_cells.integrated,idents=i)
  features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
  	'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
  p <-DoHeatmap(EC_cells.integrated, features = features, group.by = "orig.ident") 
  print(p)
}
dev.off()

#####改变小提琴横坐标的顺序#######

my_levels <- c("Arterial","C_A","Capillary","C_V","Venous","Interferon","choroid_plexus")
Idents(EC_cells.integrated) <- factor(Idents(EC_cells.integrated), levels= my_levels)
sample_levels<-c("Young","Old")
EC_cells.integrated@meta.data$orig.ident<-factor(EC_cells.integrated@meta.data$orig.ident, levels= sample_levels)

#VlnPlot(seurat.object_copy, features = c("some_gene"))

pdf("/md01/nieyg/project/BBB/plot/23genes_Violin_forYO.pdf")
 ##Violin plots can also be split on some variable. Simply add the splitting variable to object
# metadata and pass it to the split.by argument
features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
  	'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
#VlnPlot(EC_cells.integrated, features = features)
#VlnPlot(EC_cells.integrated, features = features, split.by = "orig.ident")
#dev.off()
#library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature,split.by = "orig.ident", pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )
       return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}
#配色方案
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
StackedVlnPlot(EC_cells.integrated, features[1:4], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[5:8], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[9:12], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[13:16], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[17:20], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[21:23], pt.size=0, cols=my36colors)

pdf("/md01/nieyg/project/BBB/plot/23genes_Violin_forsubtype.pdf")
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )
       return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}

StackedVlnPlot(EC_cells.integrated, features[1:4], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[5:8], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[9:12], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[13:16], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[17:20], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, features[21:23], pt.size=0, cols=my36colors)
dev.off()

#########寻找DEG########
########将Old Young cell number downsample 一致##
#################################################
old<-sample(old,7820)
down<-c(old,Young)
EC_cells.integrated@meta.data<-cellorder

EC_cells<- EC_cells.integrated[ , which(EC_cells.integrated@meta.data$cellorder %in% down)]

######各个subtype的DEG 以及统计结果######
######only up genes#####
EC.markers <- FindAllMarkers(EC_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- EC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("/md01/nieyg/project/BBB/plot/top10gene_heatmap.pdf")
#######size,angle,hjust 防止下标出界#######
p<-DoHeatmap(EC_cells.integrated, features = top10$gene,size = 2.5, angle = -50, hjust=1) + NoLegend()
p+theme(axis.text.y= element_text(size=5, color="black", face= "bold"))
dev.off()
#####UP and Down genes#####
EC.markers <- FindAllMarkers(EC_cells.integrated, min.pct = 0.25, logfc.threshold = 0.25)
#针对各个subtype的火山图：
#log2FoldChange和padj合并到一个数据框

####火山图##########
pdf("/md01/nieyg/project/BBB/plot/DEG_volcano.pdf")
for (i in levels(EC_cells.integrated)){
	mydata<-EC.markers[which(EC.markers$cluster==i),];
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"up",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"down","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	p <-ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+
              ggtitle(i)+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('up'='red','down'='blue','normal'='gray'));
    print(p)
}
dev.off()
######barplot统计结果####
|avg_logFC|>=0.5&&p_val_adj<=0.05

new.cluster.ids <- c("Capillary","CV","Arterial","Venous",
	    "CA","choroid_plexus","Interferon")
names(new.cluster.ids) <- levels(EC_cells.integrated)
EC_cells.integrated <- RenameIdents(EC_cells.integrated, new.cluster.ids)
for (i in levels(EC_cells.integrated)){
	#subtype<-levels(EC_cells.integrated)[i];
	mydata<-EC.markers[which(EC.markers$cluster==i),];
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"up",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"down","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	assign(paste("upgene_", i, sep = ""), rownames(mydata[which(mydata$Condition=="up"),]) )
	assign(paste("downgene_", i, sep = ""), rownames(mydata[which(mydata$Condition=="down"),]) )
}
#######百分比########
subtype =rep(levels(EC_cells.integrated),2)
change =c(rep("up",7),rep("down",7))
num =c(52,46,76,63,113,218,64,134,22,78,24,117,253,0)
subtype_num = data.frame(subtype,change,num)

library(plyr)
subtype_per = ddply(subtype_num, "subtype", transform, percent = num / sum(num) * 100)
pdf("/md01/nieyg/project/BBB/plot/DEG_bar_percentage.pdf")
ggplot(subtype_per, aes(x = subtype, y = percent, fill = change)) +
  # 条形图函数：未将position参数显示设置为dodge，则绘制出的条形图为堆积型
  geom_bar(stat = "identity", colour = "black") +
  theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
  # 调色标尺
  scale_fill_brewer(palette = "Pastel1")

#######
subtype_per = arrange(subtype_num, subtype, change)
# 对不同Date分组内的数据进行累加求和
subtype_per = ddply(subtype_per, "subtype", transform, label_y = cumsum(num) - 0.5*num)
 
# 基函数
ggplot(subtype_per, aes(x = subtype, y = num, fill =change)) +
  # 条形图函数
  geom_bar(stat = "identity", colour = "black") +
  # 标签函数：paste和format方法对标签进行格式化
  geom_text(aes(y=label_y, label = paste(format(num, nsmall=2), "genes")), size = 2.5) +
  # 图例函数：反转图例
  guides(fill = guide_legend(reverse = TRUE)) +
  # 调色标尺
  scale_fill_brewer(palette = "Pastel1") +theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())

######每个subtype上调的基因的upsetR图####
upgene_CV<-gsub("[.].$", "", upgene_CV)
upgene_Arterial<-gsub("[.].$", "", upgene_Arterial)
upgene_Venous<-gsub("[.].$", "", upgene_Venous)
upgene_CA<-gsub("[.].$", "", upgene_CA)  
upgene_choroid_plexus<-gsub("[.].$", "", upgene_choroid_plexus)
upgene_Interferon<-gsub("[.].$", "", upgene_Interferon)

alltermgene<-union(upgene_Capillary,upgene_CV)
alltermgene<-union(alltermgene,upgene_Arterial)
alltermgene<-union(alltermgene,upgene_Venous)
alltermgene<-union(alltermgene,upgene_CA)
alltermgene<-union(alltermgene,upgene_choroid_plexus)
alltermgene<-union(alltermgene,upgene_Interferon)
#####构建upsetR  matrix#########
meta<-matrix(data=NA,nrow=length(alltermgene),ncol=7);
colnames(meta)<-c("Capillary","CV","Arterial","Venous",
	    "CA","choroid_plexus","Interferon");
rownames(meta)<-alltermgene;

for ( i in 1:length(alltermgene))
{ 
  if (alltermgene[i] %in% upgene_Capillary)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (alltermgene[i] %in% upgene_CV)
    meta[i,2] = 1
    else 
    meta[i,2] = 0
  if (alltermgene[i] %in% upgene_Arterial)
    meta[i,3] = 1
    else 
    meta[i,3] = 0
      if (alltermgene[i] %in% upgene_Venous)
    meta[i,4] = 1
    else 
    meta[i,4] = 0
  if (alltermgene[i] %in% upgene_CA)
    meta[i,5] = 1
    else 
    meta[i,5] = 0
  if (alltermgene[i] %in% upgene_choroid_plexus)
    meta[i,6] = 1
    else 
    meta[i,6] = 0
   if (alltermgene[i] %in% upgene_Interferon)
    meta[i,7] = 1
    else 
    meta[i,7] = 0

}
######plot heatmap####
pdf("/md01/nieyg/project/BBB/plot/DEG_upsetR_Up-regulated.pdf")
#meta = data.frame(meta)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
library(UpSetR)
m<-as.data.frame(meta)
c<-cbind(alltermgene,m)
upset(c, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      nsets = 7, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Up-regulated genes",
      sets.x.label = "Total number of Up-regulated genes in each cell type", text.scale = c(2, 2, 0.8, 2, 2, 2))
######每个subtype下调的基因的upsetR图####
downgene_CV<-gsub("[.].$", "", downgene_CV)
downgene_Arterial<-gsub("[.].$", "", downgene_Arterial)
downgene_Venous<-gsub("[.].$", "", downgene_Venous)
downgene_CA<-gsub("[.].$", "", downgene_CA)
downgene_choroid_plexus<-gsub("[.].$", "", downgene_choroid_plexus)
downgene_Interferon<-gsub("[.].$", "", downgene_Interferon)

alltermgene<-union(downgene_Capillary,downgene_CV)
alltermgene<-union(alltermgene,downgene_Arterial)
alltermgene<-union(alltermgene,downgene_Venous)
alltermgene<-union(alltermgene,downgene_CA)
alltermgene<-union(alltermgene,downgene_choroid_plexus)
alltermgene<-union(alltermgene,downgene_Interferon)
#####构建upsetR  matrix#########
meta<-matrix(data=NA,nrow=length(alltermgene),ncol=7);
colnames(meta)<-c("Capillary","CV","Arterial","Venous",
	    "CA","choroid_plexus","Interferon");
rownames(meta)<-alltermgene;

for ( i in 1:length(alltermgene))
{ 
  if (alltermgene[i] %in% downgene_Capillary)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (alltermgene[i] %in% downgene_CV)
    meta[i,2] = 1
    else 
    meta[i,2] = 0
  if (alltermgene[i] %in% downgene_Arterial)
    meta[i,3] = 1
    else 
    meta[i,3] = 0
      if (alltermgene[i] %in% downgene_Venous)
    meta[i,4] = 1
    else 
    meta[i,4] = 0
  if (alltermgene[i] %in% downgene_CA)
    meta[i,5] = 1
    else 
    meta[i,5] = 0
  if (alltermgene[i] %in% downgene_choroid_plexus)
    meta[i,6] = 1
    else 
    meta[i,6] = 0
   if (alltermgene[i] %in% downgene_Interferon)
    meta[i,7] = 1
    else 
    meta[i,7] = 0

}
pdf("/md01/nieyg/project/BBB/plot/DEG_upsetR_Down-regulated.pdf")
#meta = data.frame(meta)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
library(UpSetR)
m<-as.data.frame(meta)
c<-cbind(alltermgene,m)
upset(c, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      nsets = 7, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Down-regulated genes",
      sets.x.label = "Total number of Down-regulated genes in each cell type", text.scale = c(2, 2, 0.8, 2, 2, 2))
######GO analysis for all DEGs and  common DEGs#####

for (i in levels(EC_cells.integrated)){
	print(paste("downgene_", i, sep = ""))
}
downgene<-c("downgene_Capillary","downgene_CV","downgene_Arterial","downgene_Venous","downgene_CA","downgene_choroid_plexus","downgene_Interferon")
upgene<-c("upgene_Capillary","upgene_CV","upgene_Arterial","upgene_Venous","upgene_CA","upgene_choroid_plexus","upgene_Interferon")

ego_geneID <- bitr(downgene_Capillary, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
pdf("downgene_Capillary-GOandKEGG.pdf")
ego <- enrichKEGG(
  gene = ego_geneID$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"downgene_Capillary-KEGG.csv")
ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"downgene_Capillary-BP.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"downgene_Capillary-MF.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"downgene_Capillary-CC.csv")

dev.off()

########UMAP of commonDEGS######


######Young Old之间的DEG########

############Pseudotime trajectory##################
#library(monocle3)
library(Seurat)
#library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(monocle)
EC_cells.integrated<-readRDS("BBB_onlyEC.rds")
EC_cells.integrated@meta.data$subtype<-Idents(EC_cells.integrated)
DefaultAssay(EC_cells.integrated)<-"RNA"
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(EC_cells.integrated@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EC_cells.integrated@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
 #monocle_cds<-monocle_cds
# 归一化 
 monocle_cds<- estimateSizeFactors(monocle_cds)
 monocle_cds<- estimateDispersions(monocle_cds)
###Filtering low-quality cells
 monocle_cds <- detectGenes(monocle_cds, min_expr = 3 )
 print(head(fData(monocle_cds)))

 expressed_genes <- row.names(subset(fData(monocle_cds),
                                     num_cells_expressed >= 10))
 print(head(pData(monocle_cds)))

#Clustering cells without marker genes 

disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(monocle_cds)

# Plots the percentage of variance explained by the each component based on PCA from the normalized expression data using the same procedure used in reduceDimension function.
# monocle_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(monocle_cds, return_all = F) # norm_method='log'

monocle_cds <- reduceDimension(monocle_cds, max_components = 2, num_dim = 10,
                        reduction_method = 'tSNE', verbose = T)

monocle_cds <- clusterCells(monocle_cds, num_clusters = 7)
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## 不要也写0.1 ，而是要写0.01。

monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

#Trajectory step 2: reduce data dimensionality  
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                            method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")


The23genes <- row.names(subset(fData(monocle_cds),
                               gene_short_name %in% c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
  	'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")))
plot_genes_branched_pseudotime(monocle_cds[The23genes,],
                               branch_point = 1,
                               color_by = "subtype")
#,
                               #ncol = 1)

cds_subset <- monocle_cds[The23genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "subtype")

dev.off()