######Capillary VS Choroid plexus########
## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(ArchR)
EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
#####Only Young Sample######
Idents(EC_cells.integrated)<-EC_cells.integrated$orig.ident
EC_cells.integrated_Young<-subset(x=EC_cells.integrated,idents="Young")
EC_cells.integrated_Aged<-subset(x=EC_cells.integrated,idents="Old")

Idents(EC_cells.integrated_Young)<-EC_cells.integrated_Young$celltype

###downsample#
ch<-WhichCells(EC_cells.integrated_Young,idents="Choroid_plexus")
C<-WhichCells(EC_cells.integrated_Young,idents=c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"))
#ch<-which(EC_cells.integrated@meta.data$celltype =="Choroid_plexus")
#C<-which(EC_cells.integrated@meta.data$celltype ==c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"))
C<-sample(C,159)
downsample<-c(ch,C)
#EC_cells.integrated@meta.data$cellorder<-1:30187
DefaultAssay(EC_cells.integrated)<-"RNA"
#ChpandC<- EC_cells.integrated[ , which(EC_cells.integrated@meta.data$cellorder %in% downsample)]
ChpandC<-subset(EC_cells.integrated_Young,cells=downsample)

ChpvsC_down.markers <- FindMarkers(ChpandC, ident.1 = "Choroid_plexus", ident.2 = c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"), min.pct = 0.25,logfc.threshold = 0)
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Young_ChpvsC_volcano_downsample.pdf",width=5,height=5)
  mydata2<-as.data.frame(ChpvsC_down.markers)
  mydata2$Condition=ifelse(mydata2$avg_logFC>=0.5 & mydata2$p_val_adj<=0.05,"Choroid_plexus",  
  	ifelse(mydata2$avg_logFC<=-0.5 & mydata2$p_val_adj<=0.05,"Capillary","normal"))
  ##对满足不同条件的数据给不同的标记，放入Condition列
  p <-ggplot(data=mydata2, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+xlim(c(-2, 2)) +
              ggtitle("Choroid_plexus vs Capillary(downsample)")+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Choroid_plexus'='red','Capillary'='deepskyblue','normal'='gray'));
p
dev.off()
Choroid_plexus_downsample<-rownames(mydata2[which(mydata2$Condition=="Choroid_plexus"),])
Capillary_downsample<-rownames(mydata2[which(mydata2$Condition=="Capillary"),])
#####use downsample results######
cols=c("Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29","Choroid_plexus"="#31326F")
level<-c("Capillary_1","Capillary_2","Capillary_3","Capillary_4","Choroid_plexus")
ChpandC$celltype<-factor(ChpandC$celltype,levels=level)
ChpandC@active.ident<-ChpandC$celltype
ChpandC <- ScaleData(ChpandC, features =rownames(ChpandC))

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Young_ChpvsC_heatmap_downsample.pdf")
DoHeatmap(ChpandC, features = Choroid_plexus_downsample,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1) +scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(ChpandC, features = Capillary_downsample,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1) +scale_fill_gradientn(colors = c("blue", "black", "red"))
dev.off()

#########DEGs GO and KEGG#################
library(clusterProfiler)
library(org.Mm.eg.db)

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-GO.pdf",width=13,height=8)

gene.df <- bitr(Choroid_plexus_downsample, fromType = "SYMBOL",
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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-GO-BP.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-GO-MF.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-GO-CC.csv")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-KEGG.csv")
dev.off()

#####Capillary
pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-GO.pdf",width=13,height=8)

gene.df <- bitr(Capillary_downsample, fromType = "SYMBOL",
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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-GO-BP.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-GO-MF.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-GO-CC.csv")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1)
ego
barplot(ego, showCategory=20)
write.table(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-KEGG.csv")
dev.off()

######Violinplot#########

level<-c("Capillary_1","Capillary_2","Capillary_3","Capillary_4","Choroid_plexus")
ChpandC$celltype<-factor(ChpandC$celltype,levels=level)
ChpandC@active.ident<-ChpandC$celltype
new.cluster.ids <- c("Capillary","Capillary","Capillary","Capillary","Choroid_plexus")
names(new.cluster.ids) <- levels(ChpandC)
ChpandC <- RenameIdents(ChpandC, new.cluster.ids)
ChpandC$celltype<-ChpandC@active.ident
library(patchwork)
library(magrittr)
cols=c("Capillary"="#FFCC29","Choroid_plexus"="#31326F")
pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/choroid_plexus-violin.pdf")

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VllnPlot(obj, features = feature,assay="RNA",col=cols,ncol=3,group.by="celltype",pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +geom_boxplot(width=0.2,outlier.size=0)+
               #+geom_boxplot(width=0.2)+
               theme(legend.position = "none",
               axis.text.x = element_text(face="bold", color="black", size=6,angle=0),
               #axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               #axis.ticks.y = element_line(),
               #axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )
       return(p)
}
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
       return(p)
}
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[1:12], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[13:24], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[25:36], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[37:48], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[49:60], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[61:72], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[73:84], pt.size=0)
dev.off()

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Young/Capillary-violin.pdf")
StackedVlnPlot(ChpandC, Capillary_downsample[1:12], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[13:24], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[25:36], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[37:48], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[49:60], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[61:72], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[73:84], pt.size=0)
dev.off()
#################################################
#######Aged######################################
#################################################

Idents(EC_cells.integrated_Aged)<-EC_cells.integrated_Aged$celltype

###downsample#
ch<-WhichCells(EC_cells.integrated_Aged,idents="Choroid_plexus")
C<-WhichCells(EC_cells.integrated_Aged,idents=c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"))
C<-sample(C,217)
downsample<-c(ch,C)
ChpandC<-subset(EC_cells.integrated_Aged,cells=downsample)

ChpvsC_down.markers <- FindMarkers(ChpandC, ident.1 = "Choroid_plexus", ident.2 = c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"), min.pct = 0.25,logfc.threshold = 0)
pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Aged_ChpvsC_volcano_downsample.pdf",width=5,height=5)
  mydata2<-as.data.frame(ChpvsC_down.markers)
  mydata2$Condition=ifelse(mydata2$avg_logFC>=0.5 & mydata2$p_val_adj<=0.05,"Choroid_plexus",  
  	ifelse(mydata2$avg_logFC<=-0.5 & mydata2$p_val_adj<=0.05,"Capillary","normal"))
  ##对满足不同条件的数据给不同的标记，放入Condition列
  p <-ggplot(data=mydata2, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+xlim(c(-2, 2)) +
              ggtitle("Choroid_plexus vs Capillary(downsample)")+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Choroid_plexus'='red','Capillary'='deepskyblue','normal'='gray'));
p
dev.off()
Choroid_plexus_downsample<-rownames(mydata2[which(mydata2$Condition=="Choroid_plexus"),])
Capillary_downsample<-rownames(mydata2[which(mydata2$Condition=="Capillary"),])
#####use downsample results######
cols=c("Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29","Choroid_plexus"="#31326F")
level<-c("Capillary_1","Capillary_2","Capillary_3","Capillary_4","Choroid_plexus")
ChpandC$celltype<-factor(ChpandC$celltype,levels=level)
ChpandC@active.ident<-ChpandC$celltype
ChpandC <- ScaleData(ChpandC, features =rownames(ChpandC))

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Aged_ChpvsC_heatmap_downsample.pdf")
DoHeatmap(ChpandC, features = Choroid_plexus_downsample,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1) +scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(ChpandC, features = Capillary_downsample,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1) +scale_fill_gradientn(colors = c("blue", "black", "red"))
dev.off()

#########DEGs GO and KEGG#################
pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-GO.pdf",width=13,height=8)
gene.df <- bitr(Choroid_plexus_downsample, fromType = "SYMBOL",
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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-GO-BP.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-GO-MF.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-GO-CC.csv")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-KEGG.csv")
dev.off()

#####Capillary
pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-GO.pdf",width=13,height=8)

gene.df <- bitr(Capillary_downsample, fromType = "SYMBOL",
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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-GO-BP.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-GO-MF.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-GO-CC.csv")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1)
ego
barplot(ego, showCategory=20)
write.table(ego,"/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-KEGG.csv")
dev.off()

######Violinplot#########
new.cluster.ids <- c("Capillary","Capillary","Capillary","Capillary","Choroid_plexus")
names(new.cluster.ids) <- levels(ChpandC)
ChpandC <- RenameIdents(ChpandC, new.cluster.ids)
ChpandC$celltype<-ChpandC@active.ident
cols=c("Capillary"="#FFCC29","Choroid_plexus"="#31326F")
pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/choroid_plexus-violin.pdf")
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[1:12], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[13:24], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[25:36], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[37:48], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[49:60], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[61:72], pt.size=0)
StackedVlnPlot(ChpandC, Choroid_plexus_downsample[73:84], pt.size=0)
dev.off()

pdf("/md01/nieyg/project/BBB/YMO_results/choroid_plexus/Aged/Capillary-violin.pdf")
StackedVlnPlot(ChpandC, Capillary_downsample[1:12], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[13:24], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[25:36], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[37:48], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[49:60], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[61:72], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[73:84], pt.size=0)
dev.off()

