## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

########将Old Young cell number downsample 一致##
#################################################
EC_cells.integrated<-readRDS("BBB_onlyEC.rds")
old<-which(EC_cells.integrated@meta.data$orig.ident =="Old")
young<-which(EC_cells.integrated@meta.data$orig.ident =="Young")
old<-sample(old,8820)
down<-c(young,old)
EC_cells.integrated@meta.data$cellorder<-1:21632
DefaultAssay(EC_cells.integrated)<-"RNA"
DefaultAssay(EC_cells)<-"RNA"
EC_cells<- EC_cells.integrated[ , which(EC_cells.integrated@meta.data$cellorder %in% down)]

EC_cells<-readRDS("EC_cells_downsample.rds")


######各个subtype的DEG 以及统计结果######
######only up genes#####
EC.markers <- FindAllMarkers(EC_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
EC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- EC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("/md01/nieyg/project/BBB/plot/top10gene_heatmap_downsample.pdf")
#######size,angle,hjust 防止下标出界#######
p<-DoHeatmap(EC_cells, features = top10$gene,size = 2.5, angle = -50, hjust=1) + NoLegend()
p+theme(axis.text.y= element_text(size=5, color="black", face= "bold"))
dev.off()
#####UP and Down genes#####
EC.markers <- FindAllMarkers(EC_cells, min.pct = 0.25, logfc.threshold = 0)


#######提取metadata矩阵,更改active.ident######

new.cluster.ids <- c("Capillary","C_V","C_A","Capillary","Arterial","Venous",
	    "Capillary","C_A","choroid_plexus","Interferon")
names(new.cluster.ids) <- levels(EC_cells)
EC_cells<- RenameIdents(EC_cells, new.cluster.ids)
EC_cells@meta.data$subtype<-Idents(EC_cells)
metadata<-EC_cells@meta.data[,c(1,9,6)]
sample<-as.character(metadata$subtype)
time<-as.character(metadata$orig.ident)
sub<-as.character(metadata$seurat_clusters)
for (i in 1:17640){
if(sample[i]=="Capillary"&&time[i]=="Young"){sub[i]="Capillary_Young"}
if(sample[i]=="C_V"&&time[i]=="Young"){sub[i]="C_V_Young"}
if(sample[i]=="Arterial"&&time[i]=="Young"){sub[i]="Arterial_Young"}
if(sample[i]=="Venous"&&time[i]=="Young"){sub[i]="Venous_Young"}
if(sample[i]=="C_A"&&time[i]=="Young"){sub[i]="C_A_Young"}
if(sample[i]=="choroid_plexus"&&time[i]=="Young"){sub[i]="choroid_plexus_Young"}
if(sample[i]=="Interferon"&&time[i]=="Young"){sub[i]="Interferon_Young"}
if(sample[i]=="Capillary"&&time[i]=="Old"){sub[i]="Capillary_Old"}
if(sample[i]=="C_V"&&time[i]=="Old"){sub[i]="C_V_Old"}
if(sample[i]=="Arterial"&&time[i]=="Old"){sub[i]="Arterial_Old"}
if(sample[i]=="Venous"&&time[i]=="Old"){sub[i]="Venous_Old"}
if(sample[i]=="C_A"&&time[i]=="Old"){sub[i]="C_A_Old"}
if(sample[i]=="choroid_plexus"&&time[i]=="Old"){sub[i]="choroid_plexus_Old"}
if(sample[i]=="Interferon"&&time[i]=="Old"){sub[i]="Interferon_Old"}
}
sub<-as.factor(sub)
names(sub)<-rownames(metadata)
EC_cells@active.ident <- sub
head(EC_cells@active.ident)
levels(metadata$subtype)
#########寻找各个subtype内部，YO之间的DEG#################
"Capillary"      "C-V"            "Arterial"       "Venous"        
"C-A"            "choroid plexus" "Interferon"   
        Arterial_Old       Arterial_Young              C_A_Old 
                 818                  872                 1232 
           C_A_Young              C_V_Old            C_V_Young 
                1703                 1115                 1867 
       Capillary_Old      Capillary_Young   choroid_plexus_Old 
                4884                 3772                  152 
choroid_plexus_Young       Interferon_Old     Interferon_Young 
                 158                   83                  142 
          Venous_Old         Venous_Young 
                 536                  306 

Caplillary.markers <- FindMarkers(EC_cells, ident.1 = "Capillary_Young", ident.2 = "Capillary_Old", min.pct = 0.25,logfc.threshold = 0)
CV.markers <- FindMarkers(EC_cells, ident.1 = "C_V_Young", ident.2 = "C_V_Old", min.pct = 0.25,logfc.threshold = 0)
CA.markers <- FindMarkers(EC_cells, ident.1 = "C_A_Young", ident.2 = "C_A_Old", min.pct = 0.25,logfc.threshold = 0)
Arterial.markers <- FindMarkers(EC_cells, ident.1 = "Arterial_Young", ident.2 = "Arterial_Old", min.pct = 0.25,logfc.threshold = 0)
Venous.markers <- FindMarkers(EC_cells, ident.1 = "Venous_Young", ident.2 = "Venous_Old", min.pct = 0.25,logfc.threshold = 0)
choroid.markers <- FindMarkers(EC_cells, ident.1 = "choroid_plexus_Young", ident.2 = "choroid_plexus_Old", min.pct = 0.25,logfc.threshold = 0)
Interferon.markers <- FindMarkers(EC_cells, ident.1 = "Interferon_Young", ident.2 = "Interferon_Old", min.pct = 0.25,logfc.threshold = 0)

mylist<-list(Caplillary.markers,CV.markers,Arterial.markers,Venous.markers,CA.markers,choroid.markers,Interferon.markers)
#####画火山图#####
pdf("/md01/nieyg/project/BBB/plot/DEG_volcano_subtype_YO_downsample.pdf",width=5,height=5)
for (i in 1:7){
	mydata<-mylist[i]
	mydata<-as.data.frame(mydata)
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Young",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Old","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	p <-ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("ln fold change") + ylab("-log10 padj")+
              xlim(c(-2, 2)) +
              ggtitle(levels(metadata$subtype)[i])+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Old'='red','Young'='blue','normal'='gray'));
    print(p)
}
dev.off()

#######百分比########

m<-c("Caplillary","C_V","Arterial","Venous","C_A","choroid_plexus","Interferon")
for (i in 1:7){
	mydata<-mylist[i]
	mydata<-as.data.frame(mydata)
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Young",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Old","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	assign(paste("Young_",m[i], sep = ""), rownames(mydata[which(mydata$Condition=="Young"),]) )
	assign(paste("Old_",m[i], sep = ""), rownames(mydata[which(mydata$Condition=="Old"),]) )
}


Young_Capillary
Young_C_V
Young_Arterial
Young_Venous
Young_C_A
Young_choroid_plexus
Young_Interferon
Old_Capillary
Old_C_V
Old_Arterial
Old_Venous
Old_C_A
Old_choroid_plexus
Old_Interferon


subtype =rep(levels(metadata$subtype),2)
change =c(rep("Young",7),rep("Old",7))
num =c(9,3,5,7,3,3,6,36,21,21,15,19,13,12)
subtype_num = data.frame(subtype,change,num)

library(plyr)
subtype_per = ddply(subtype_num, "subtype", transform, percent = num / sum(num) * 100)
pdf("/md01/nieyg/project/BBB/plot/DEG_bar_percentage_YO_downsample.pdf")
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


alltermgene<-union(Old_Capillary,Old_C_V)
alltermgene<-union(alltermgene,Old_Arterial)
alltermgene<-union(alltermgene,Old_Venous)
alltermgene<-union(alltermgene,Old_C_A)
alltermgene<-union(alltermgene,Old_choroid_plexus)
alltermgene<-union(alltermgene,Old_Interferon)

all_Young_gene:
 [1] "Macf1"    "AY036118" "Gm42418"  "Ptn"      "Utrn"     "Btg2"    
 [7] "Tpm1"     "Irf1"     "Ndnf"     "Csf1"     "Dusp1"    "Cxcl1"   
[13] "Ivns1abp" "Hmcn1"    "Lars2"    "Sparc"    "Apcdd1"   "Col4a1"  
[19] "Txnip"    "Arhgap29" "Ahnak"    "Dab2"     "Hlx"     

all_Old_gene:

 [1] "H2-Q7"    "Ttr"      "H2-Q6"    "Depp1"    "H2-K1"    "Hbb-bs"  
 [7] "Lypd1"    "Tmem252"  "Mgp"      "Gadd45g"  "Adamts9"  "Clu"     
[13] "Litaf"    "S100a6"   "Socs3"    "Pglyrp1"  "Ifi27l2a" "Edn1"    
[19] "Hba-a1"   "Mt1"      "Hba-a2"   "Igfbp7"   "Ly6a"     "Ly6c1"   
[25] "Fxyd5"    "Acer2"    "Scgb3a1"  "Mt2"      "Ctla2a"   "Srarp"   
[31] "Grrp1"    "Ifitm2"   "Lcn2"     "Aldoa"    "Fas"      "Cd200"   
[37] "Xdh"      "Ucp2"     "Podxl"    "Xbp1"     "Vwf"      "Cp"      
[43] "Cxcl12"   "Lrg1"     "Kansl1l"  "Sbno2"    "Map3k6"   "Fmo2"    
[49] "Galnt15"  "Slc16a1"  "Plat"     "Kbtbd11"  "Selenop"  "Tgm2"    
[55] "Ier3"     "Stra6"    "Bnip3"    "Cd9"      "Stmn2"    "Plaur"   
[61] "Cystm1"   "Ccn1"     "H2-Q4"    "Bhlhe40"  "Slfn2"    "Cyp26b1" 
[67] "Inhbb"    "Car4"     "Cdkn1a"   "Ubc"      "Por"      "Bsg" 

#####构建upsetR  matrix#########
meta<-matrix(data=NA,nrow=length(alltermgene),ncol=7);
colnames(meta)<-c("Capillary","C_V","Arterial","Venous",
	    "C_A","choroid_plexus","Interferon");
rownames(meta)<-alltermgene;

for ( i in 1:length(alltermgene))
{ 
  if (alltermgene[i] %in% Old_Capillary)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (alltermgene[i] %in% Old_C_V)
    meta[i,2] = 1
    else 
    meta[i,2] = 0
  if (alltermgene[i] %in% Old_Arterial)
    meta[i,3] = 1
    else 
    meta[i,3] = 0
      if (alltermgene[i] %in% Old_Venous)
    meta[i,4] = 1
    else 
    meta[i,4] = 0
  if (alltermgene[i] %in% Old_C_A)
    meta[i,5] = 1
    else 
    meta[i,5] = 0
  if (alltermgene[i] %in% Old_choroid_plexus)
    meta[i,6] = 1
    else 
    meta[i,6] = 0
   if (alltermgene[i] %in% Old_Interferon)
    meta[i,7] = 1
    else 
    meta[i,7] = 0
}
pdf("/md01/nieyg/project/BBB/plot/DEG_upsetR_Old.pdf")
#meta = data.frame(meta)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
library(UpSetR)
m<-as.data.frame(meta)
c<-cbind(alltermgene,m)
upset(c, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      nsets = 7, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Old genes",
      sets.x.label = "Total number of Old genes in each cell type", text.scale = c(2, 2, 0.8, 2, 2, 2))


mart<-read.csv("gene-description-mart_export.csv")

YC<-mart[which(mart$Gene.name%in%Old_choroid_plexus),]
write.csv(YC,"Old_choroid_plexus_annotation.csv")






pdf("Old_choroid_plexus-GO-BP.pdf",width=13,height=8)

gene.df <- bitr(Old_choroid_plexus, fromType = "SYMBOL",
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
write.csv(ego,"Old_C_A-GO-BP.csv")

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
write.csv(ego,"Old_C_A-GO-MF.csv")

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
write.csv(ego,"Old_C_A-GO-CC.csv")








> Young_Capillary
[1] "Ivns1abp" "Hmcn1"    "Lars2"    "Gm42418"  "AY036118"
> Young_C_V
[1] "Macf1"    "Col4a1"   "Hmcn1"    "Apcdd1"   "Ivns1abp" "Ndnf"     "Gm42418" 
> Young_Arterial
[1] "Ptn"   "Utrn"  "Btg2"  "Tpm1"  "Irf1"  "Ndnf"  "Csf1"  "Dusp1" "Cxcl1"
> Young_Venous
[1] "Macf1"  "Sparc"  "Apcdd1"
> Young_C_A
[1] "Macf1"    "AY036118" "Gm42418" 
> Young_choroid_plexus
[1] "Col4a1"   "Txnip"    "Arhgap29" "Ahnak"    "Dab2"     "Hlx"     
> Young_Interferon
[1] "Macf1"  "Sparc"  "Col4a1"
> Old_Capillary
"Ttr"     "Ly6a"    "Acer2"   "Fmo2"    "Depp1"   "Hbb-bs"  "Xdh" "Pglyrp1" "Tmem252" "Plat"    "Kbtbd11" "Selenop" "Vwf"     "H2-K1"  "Fxyd5"   "Grrp1"   "Map3k6"  "Cp"      "Tgm2"    "Ier3"    "Stra6"  
> Old_C_V
"Ttr"     "Ly6a"    "Acer2"   "Depp1"   "Fmo2"    "Tmem252" "Aldoa"  "Plat"    "Hbb-bs"  "H2-K1"   "Xdh"     "Cystm1"  "Vwf"     "Stra6"  "Ccn1"   
> Old_Arterial
"Igfbp7"  "Ttr"     "Ly6a"    "Ly6c1"   "Fxyd5"   "Acer2"   "Scgb3a1" "Mt2"     "Ctla2a"  "Srarp"   "Grrp1"   "H2-Q7"   "Litaf"   "Mt1" "Ifitm2"  "Lcn2"    "Aldoa"   "Fas"     "Cd200"   "Tmem252" "Xdh"  "Ucp2"    "Podxl"   "Xbp1"    "Vwf"     "Adamts9" "H2-Q6"   "Cp"     "Cxcl12"  "Lrg1"  "Kansl1l" "Sbno2"   "Map3k6"  "Fmo2"    "Galnt15" "Slc16a1"
> Old_Venous
"Ttr"     "Acer2"   "Ly6a"    "Tmem252" "Fmo2"    "Aldoa"   "Hbb-bs" "Litaf"   "Depp1"   "Plat"    "H2-K1"   "Bnip3"   "Cd9"     "Stmn2"  "Vwf"     "Plaur"   "Gadd45g" "Mt1"     "Mt2"    
> Old_C_A
"H2-Q7"    "Ttr"      "H2-Q6"    "Depp1"    "H2-K1"    "Hbb-bs"  "Lypd1"    "Tmem252"  "Mgp"      "Gadd45g"  "Adamts9"  "Clu"  "Litaf"    "S100a6"   "Socs3"    "Pglyrp1"  "Ifi27l2a" "Edn1"    "Hba-a1"   "Mt1"      "Hba-a2"  
_choroid_plexus
"Cyp26b1" "H2-Q7"   "Ttr"     "H2-Q6"   "Acer2"   "Grrp1"   "Inhbb"  "Car4"    "Cdkn1a"  "Ubc"     "Por"     "Bsg"    
> Old_Interferon
"H2-Q7"   "H2-Q6"   "Ttr"     "Depp1"   "Ly6a"    "H2-Q4"   "Hbb-bs" "Acer2"   "Fmo2"    "Bhlhe40" "Slfn2"   "H2-K1"   "Xdh"    
> 
