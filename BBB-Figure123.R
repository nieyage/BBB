######BBB-Figure1#######
## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(ArchR)
EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
BBB_EC.integrated<-readRDS("../BBB_YMO.rds")
DefaultAssay(BBB_EC.integrated) <- "RNA" 
pdf("/md01/nieyg/project/BBB/YMO_results/EC_marker_UMAP.pdf")
m<-c("Pecam1", "Cdh5","Tek")
FeaturePlot(BBB_EC.integrated, features =m, pt.size = 0.2)
dev.off()
cols=c("Arterial"="#206A5D","C_A"="#81B264","Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29",
    "C_V_1"="#F58634","C_V_2"="#F58634","Venous"="#BE0000","Choroid_plexus"="#31326F","Interferon"="#93329E")
level<-c("Young","Middle1","Old")
EC_cells.integrated$orig.ident<-factor(EC_cells.integrated$orig.ident,levels=level)
level<-c("Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4",
    "C_V_1","C_V_2","Venous","Choroid_plexus","Interferon")
EC_cells.integrated$celltype<-factor(EC_cells.integrated$celltype,levels=level)
EC_cells.integrated@active.ident<-EC_cells.integrated$celltype
#####The 23 BBB genes######
features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
    'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
EC_cells.integrated <- ScaleData(EC_cells.integrated, features =rownames(EC_cells.integrated))
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/23genes_heatmap_forsubtype-V8.pdf")
p<-DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1.5,disp.max=0.5)+scale_fill_gradientn(colors = c("blue", "white", "red"))
p
i=15
for(i in 15:30){
col.pal <- colorRampPalette(ArchRPalettes[[26]])
fill.colors <- col.pal(64)
p<-DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1.5,disp.max=0.5)+scale_fill_gradientn(colors = fill.colors)
print(p)
}
dev.off()
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/23genes_heatmap_forsubtype-V9.pdf")
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = fill.colors)
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#FA00FA", "black", "#F0F000"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("navy", "black", "firebrick3"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#1873CC", "#00CC00", "#FFFF00"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#352A86", "#2DB7A3", "#F8FA0D"))
dev.off()

#####violin plot#######
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VllnPlot(obj, features = feature,assay="RNA",col=cols,group.by="celltype",pt.size = pt.size, ... ) +
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
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(monocle)
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/23genes_and_4TFs_heatmap_forsubtype.pdf")
StackedVlnPlot(EC_cells.integrated, features[1:4], pt.size=0)
StackedVlnPlot(EC_cells.integrated, features[5:8], pt.size=0)
StackedVlnPlot(EC_cells.integrated, features[9:12], pt.size=0)
StackedVlnPlot(EC_cells.integrated, features[13:16], pt.size=0)
StackedVlnPlot(EC_cells.integrated, features[17:20], pt.size=0)
StackedVlnPlot(EC_cells.integrated, features[21:23], pt.size=0)
TF<-c("Tcf7","Lef1","Erg","Sox17")
StackedVlnPlot(EC_cells.integrated, TF[1:4], pt.size=0)
dev.off()

#######Figure2 ##########
Young_DEG<-c("Macf1","Hmcn1"  ,"Ivns1abp","Lars2" ,     
"Gm42418"      ,"Ndnf"   ,"AY036118","Zfp655",     
"Sacs"         ,"Filip1" ,"Atrx"    ,"Sptbn1",     
"C130074G19Rik","Tpm1"   ,"Sparc"   ,"Apcdd1",     
"Col4a1"       ,"Ccn3"   ,"Ptn"     ,"Btg2"  ,     
"Utrn"         ,"Zfp36l1","Dleu2"   ,"Cxcl1" ,     
"Csf1"         ,"Txnip"  ,"Ahnak"   ,"Dab2"  ,     
"H1f0"         ,"Hbb-bs")
Aged_DEG<-c(
"Ly6a"   ,"Acer2"  ,"Ttr"    ,"Fmo2"    ,"Pglyrp1","Xdh"    ,
"Depp1"  ,"Selenop","Tmem252","Hbb-bs"  ,"Fxyd5"  ,"Kbtbd11",
"Plat"   ,"Vwf"    ,"H2-K1"  ,"Cp"      ,"Grrp1"  ,"Map3k6" ,
"Apold1" ,"Stra6"  ,"Bnip3"  ,"Ier3"    ,"Tgm2"   ,"Actg1"  ,
"Slc2a1" ,"Ly6c1"  ,"B2m"    ,"Actb"    ,"Ccn1"   ,"Timp3"  ,
"Pmp22"  ,"Klf2"   ,"Mat2a"  ,"Cdkn1a"  ,"Clic4"  ,"Nr4a1"  ,
"Kdm6b"  ,"Adamts1","Bhlhe40","Lrrc8a"  ,"Snrk"   ,"Crip1"  ,
"Tinagl1","Id1"    ,"Ndrg1"  ,"Sgk1"    ,"Rhob"   ,"H2-T23" ,
"Ifitm2" ,"Per1"   ,"Rgcc"   ,"Itih5"   ,"Rhoc"   ,"Jund"   ,
"H2-D1"  ,"Ctla2a" ,"Aldoa"  ,"Malat1"  ,"H2-Q7"  ,"Litaf"  ,
"Plaur"  ,"Gadd45g","Mt1"    ,"Mt2"     ,"Scgb3a1","Cebpd"  ,
"Lypd1"  ,"H2-Q6"  ,"S100a6" ,"Ifi27l2a","Socs3"  ,"Edn1"   ,
"Igfbp7" ,"Lcn2"   ,"Cd200"  ,"Srarp"   ,"Podxl"  ,"Fas"    ,
"Lrg1"   ,"Ucp2"   ,"Xbp1"   ,"Adamts9" ,"Osmr"   ,"Gapdh"  ,
"Plscr1" ,"Cxcl12" ,"Sbno2"  ,"Il1r1"   ,"Kansl1l","Gm47283",
"Hif1a"  ,"Cd14"   ,"Slc16a1","Ldha"    ,"Cd9"    ,"Stmn2"  ,
"Mgp"    ,"Junb"   ,"Fos"    ,"Cyp26b1" ,"Inhbb"  ,"Car4"   ,
"Bsg"    ,"Por"    ,"Id4"    ,"Hba-a1"  ,"Hba-a2")
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/YA_DEG_heatmap_forYMA.pdf")
cols<-c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35")
DoHeatmap(EC_cells.integrated,Young_DEG,group.by = "orig.ident",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(EC_cells.integrated,Aged_DEG,group.by = "orig.ident",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
dev.off();
Caplillary2.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Capillary_2", min.pct = 0.25,logfc.threshold = 0)
  mydata<-as.data.frame(Caplillary2.markers)
  mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Up",  ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Down","normal"))
Capillary_2_top<-rownames(mydata[which(mydata$Condition=="Up"),])
> Capillary_2_top
 [1] "mt-Co3"  "mt-Co2"  "mt-Atp6" "mt-Cytb" "mt-Nd4"  "Hbb-bs"  "mt-Nd2" 
 [8] "Igf1r"   "Ttr"     "mt-Nd1"  "Ccdc141" "Xdh"     "mt-Nd3"  "Mcf2l"  
[15] "Hba-a1"  "Hba-a2"  "Kank3"   "Fnbp1l"  "Slc38a3" "Gm47283" "Rasgrp2"

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/C2vsotherallEC_heatmap.pdf")
cols=c("Arterial"="#206A5D","C_A"="#81B264","Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29",
    "C_V_1"="#F58634","C_V_2"="#F58634","Venous"="#BE0000","Choroid_plexus"="#31326F","Interferon"="#93329E")
DoHeatmap(EC_cells.integrated,Capillary_2_top,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
cols<-c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35")
DoHeatmap(EC_cells.integrated,Capillary_2_top,group.by = "orig.ident",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))


Caplillary2vs3.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Capillary_2", ident.2 = "Capillary_3", min.pct = 0.25,logfc.threshold = 0)

library(ggrepel)
pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/C2vsC3_volcano-3.pdf",width=5,height=5)
  mydata<-as.data.frame(Caplillary2vs3.markers)
  mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Capillary_2",  ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Capillary_3","normal"))
  ##对满足不同条件的数据给不同的标记，放入Condition列
  mydata$gene<-rownames(mydata)
  ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+xlim(c(-2, 2)) +
              ggtitle("Capillary_2 vs Capillary_3")+
              geom_text_repel(data=mydata[mydata$Condition %in% c("Capillary_2","Capillary_3"),],segment.color = "black", aes(label=gene),size=3,show.legend=FALSE)+
              theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Capillary_2'='red','Capillary_3'='deepskyblue','normal'='gray'));

color = "black",





dev.off()
Capillary_2<-rownames(mydata[which(mydata$Condition=="Capillary_2"),])
Capillary_3<-rownames(mydata[which(mydata$Condition=="Capillary_3"),])
Capillary_2<-c("Ly6c1"  ,"Slc2a1","Ly6a"  ,"Ttr"   ,"Xdh"    ,"Hbb-bs" ,"Igfbp7",
"Selenop","Cdkn1a","Acer2" ,"Plat"  ,"Pglyrp1","Tmem252","Fnbp1l",
"Grrp1"  ,"Ndrg1" ,"Hba-a1","Ctla2a","Mcf2l"  ,"Fas"    ,"Mal"   ,
"Vwf"    ,"Osmr"  ,"Hba-a2","Fth1"  ,"Mt1")    
Capillary_3<-c(
"Gm42418","Utrn","Macf1","Lars2","Syne1","Atrx",    "Dusp1"  ,
"Btg2"   ,"Plk2","Ccn1" ,"Tfrc"   )

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/C2vsC3_heatmap.pdf",width=5,height=5)
cols=c("Arterial"="#206A5D","C_A"="#81B264","Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29",
    "C_V_1"="#F58634","C_V_2"="#F58634","Venous"="#BE0000","Choroid_plexus"="#31326F","Interferon"="#93329E")
DoHeatmap(Capillary_2vs3,Capillary_2,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(Capillary_2vs3,Capillary_3,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
cols<-c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35")
DoHeatmap(Capillary_2vs3,Capillary_2,group.by = "orig.ident",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(Capillary_2vs3,Capillary_3,group.by = "orig.ident",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
dev.off()
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/C2vsC3_violin_Capillary_2vs3.pdf",width=5,height=5)

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VllnPlot(obj, features = feature,assay="RNA",col=cols,group.by="celltype",split.by="orig.ident",pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +geom_boxplot(width=0.2,outlier.size=0,position=position_dodge(width=0.9))+
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
StackedVlnPlot(Capillary_2vs3, Capillary_2[1:4], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_2[5:8], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_2[9:12], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_2[13:16], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_2[17:20], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_2[21:24], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_2[25:26], pt.size=0)

StackedVlnPlot(Capillary_2vs3, Capillary_3[1:4], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_3[5:8], pt.size=0)
StackedVlnPlot(Capillary_2vs3, Capillary_3[9:11], pt.size=0)
dev.off()


> intersect(Capillary_2_top,Capillary_2)
[1] "Hbb-bs" "Ttr"    "Xdh"    "Mcf2l"  "Hba-a1" "Hba-a2" "Fnbp1l"


#######Figure3########
ChpvsC.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Choroid_plexus", ident.2 = c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"), min.pct = 0.25,logfc.threshold = 0)

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/ChpvsC_volcano.pdf",width=5,height=5)
  mydata<-as.data.frame(ChpvsC.markers)
  mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Choroid_plexus",  ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Capillary","normal"))
  ##对满足不同条件的数据给不同的标记，放入Condition列
  p <-ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+xlim(c(-2, 2)) +
              ggtitle("Choroid_plexus vs Capillary")+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Choroid_plexus'='red','Capillary'='deepskyblue','normal'='gray'));
p
dev.off()
Choroid_plexus<-rownames(mydata[which(mydata$Condition=="Choroid_plexus"),])
Capillary<-rownames(mydata[which(mydata$Condition=="Capillary"),])

###downsample#
ch<-which(EC_cells.integrated@meta.data$celltype =="Choroid_plexus")
C<-which(EC_cells.integrated@meta.data$celltype ==c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"))
C<-sample(C,483)
downsample<-c(ch,C)
EC_cells.integrated@meta.data$cellorder<-1:30187
DefaultAssay(EC_cells.integrated)<-"RNA"
ChpandC<- EC_cells.integrated[ , which(EC_cells.integrated@meta.data$cellorder %in% downsample)]
ChpvsC_down.markers <- FindMarkers(ChpandC, ident.1 = "Choroid_plexus", ident.2 = c("Capillary_1","Capillary_2","Capillary_3","Capillary_4"), min.pct = 0.25,logfc.threshold = 0)
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/ChpvsC_volcano_downsample.pdf",width=5,height=5)
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

Chp<-intersect(Choroid_plexus_downsample,Choroid_plexus)
C<-intersect(Capillary_downsample,Capillary)
#####use downsample results######
cols=c("Arterial"="#206A5D","C_A"="#81B264","Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29",
    "C_V_1"="#F58634","C_V_2"="#F58634","Venous"="#BE0000","Choroid_plexus"="#31326F","Interferon"="#93329E")
level<-c("Arterial","C_A","C_V_1","C_V_2","Venous","Interferon","Capillary_1","Capillary_2","Capillary_3","Capillary_4",
    "Choroid_plexus")
EC_cells.integrated$celltype<-factor(EC_cells.integrated$celltype,levels=level)
EC_cells.integrated@active.ident<-EC_cells.integrated$celltype

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/ChpvsC_heatmap_downsample.pdf")
DoHeatmap(ChpandC, features = Choroid_plexus_downsample,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1) +scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(ChpandC, features = Capillary_downsample,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1) +scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
dev.off()

#########DEGs GO and KEGG#################
library(clusterProfiler)
library(org.Mm.eg.db)

Choroid_plexus_downsample<-rownames(mydata2[which(mydata2$Condition=="Choroid_plexus"),])
Capillary_downsample<-rownames(mydata2[which(mydata2$Condition=="Capillary"),])

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-GO.pdf",width=13,height=8)

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-GO-BP.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-GO-MF.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-GO-CC.csv")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-KEGG.csv")
dev.off()

#####Capillary
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-GO.pdf",width=13,height=8)

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-GO-BP.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-GO-MF.csv")

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
write.csv(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-GO-CC.csv")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1)
ego
barplot(ego, showCategory=20)
write.table(ego,"/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-KEGG.csv")
dev.off()


######Violinplot#########

ChpandC
new.cluster.ids <- c("Capillary","Capillary","Capillary","Capillary","Choroid_plexus")
names(new.cluster.ids) <- levels(ChpandC)
ChpandC <- RenameIdents(ChpandC, new.cluster.ids)
ChpandC$celltype<-ChpandC@active.ident

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/choroid_plexus_GOKEGG/choroid_plexus-violin.pdf")
cols=c("Capillary"="#FFCC29","Choroid_plexus"="#31326F")

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



pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_GOKEGG/Capillary-violin.pdf")
cols=c("Capillary"="#FFCC29","Choroid_plexus"="#31326F")

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
StackedVlnPlot(ChpandC, Capillary_downsample[1:12], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[13:24], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[25:36], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[37:48], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[49:60], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[61:72], pt.size=0)
StackedVlnPlot(ChpandC, Capillary_downsample[73:84], pt.size=0)
dev.off()

C2top<-c("mt-Co3","mt-Co2","mt-Atp6","mt-Cytb","mt-Nd4" ,"Hbb-bs" ,"mt-Nd2" ,
"Igf1r" ,"Ttr"   ,"mt-Nd1" ,"Ccdc141","Xdh"    ,"mt-Nd3" ,"Mcf2l"  ,
"Hba-a1","Hba-a2","Kank3"  ,"Fnbp1l" ,"Slc38a3","Gm47283","Rasgrp2")
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Capillary_2_top-violin.pdf")
StackedVlnPlot(EC_cells.integrated, C2top[1:4], pt.size=0)
StackedVlnPlot(EC_cells.integrated, C2top[5:8], pt.size=0)
StackedVlnPlot(EC_cells.integrated, C2top[9:12], pt.size=0)
StackedVlnPlot(EC_cells.integrated, C2top[13:16], pt.size=0)
StackedVlnPlot(EC_cells.integrated, C2top[17:20], pt.size=0)
StackedVlnPlot(EC_cells.integrated, C2top[21], pt.size=0)
dev.off();


Aged_DEG<-c("Bsg","Slc16a1","Slc2a1")
Young_DEG<-c("Ptn","Utrn","Ndnf","Ahnak")

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Aged_DEG_BBB-violin.pdf")
StackedVlnPlot(EC_cells.integrated, Aged_DEG[1:4], pt.size=0)

dev.off();
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Young_DEG_BBB-violin.pdf")
StackedVlnPlot(EC_cells.integrated, Young_DEG[1:4], pt.size=0)
dev.off();