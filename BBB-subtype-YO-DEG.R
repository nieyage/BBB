#######提取metadata矩阵,更改active.ident######
metadata<-EC_cells.integrated@meta.data[,c(1,8,6)]
sample<-as.character(metadata$subtype)
time<-as.character(metadata$orig.ident)
sub<-as.character(metadata$seurat_clusters)
for (i in 1:21632){
if(sample[i]=="Capillary"&&time[i]=="Young"){sub[i]="Capillary-Young"}
if(sample[i]=="C-V"&&time[i]=="Young"){sub[i]="C-V-Young"}
if(sample[i]=="Arterial"&&time[i]=="Young"){sub[i]="Arterial-Young"}
if(sample[i]=="Venous"&&time[i]=="Young"){sub[i]="Venous-Young"}
if(sample[i]=="C-A"&&time[i]=="Young"){sub[i]="C-A-Young"}
if(sample[i]=="choroid plexus"&&time[i]=="Young"){sub[i]="choroid-Young"}
if(sample[i]=="Interferon"&&time[i]=="Young"){sub[i]="Interferon-Young"}
if(sample[i]=="Capillary"&&time[i]=="Old"){sub[i]="Capillary-Old"}
if(sample[i]=="C-V"&&time[i]=="Old"){sub[i]="C-V-Old"}
if(sample[i]=="Arterial"&&time[i]=="Old"){sub[i]="Arterial-Old"}
if(sample[i]=="Venous"&&time[i]=="Old"){sub[i]="Venous-Old"}
if(sample[i]=="C-A"&&time[i]=="Old"){sub[i]="C-A-Old"}
if(sample[i]=="choroid plexus"&&time[i]=="Old"){sub[i]="choroid-Old"}
if(sample[i]=="Interferon"&&time[i]=="Old"){sub[i]="Interferon-Old"}
}
sub<-as.factor(sub)
names(sub)<-rownames(metadata)
EC_cells.integrated@active.ident <- sub
head(EC_cells.integrated@active.ident)
levels(metadata$subtype)
#########寻找各个subtype内部，YO之间的DEG#################
"Capillary"      "C-V"            "Arterial"       "Venous"        
"C-A"            "choroid plexus" "Interferon"   

Caplillary.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Capillary-Young", ident.2 = "Capillary-Old", min.pct = 0.25)
CV.markers <- FindMarkers(EC_cells.integrated, ident.1 = "C-V-Young", ident.2 = "C-V-Old", min.pct = 0.25)
CA.markers <- FindMarkers(EC_cells.integrated, ident.1 = "C-A-Young", ident.2 = "C-A-Old", min.pct = 0.25)
Arterial.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Arterial-Young", ident.2 = "Arterial-Old", min.pct = 0.25)
Venous.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Venous-Young", ident.2 = "Venous-Old", min.pct = 0.25)
choroid.markers <- FindMarkers(EC_cells.integrated, ident.1 = "choroid-Young", ident.2 = "choroid-Old", min.pct = 0.25)
Interferon.markers <- FindMarkers(EC_cells.integrated, ident.1 = "Interferon-Young", ident.2 = "Interferon-Old", min.pct = 0.25)

mylist<-list(Caplillary.markers,CV.markers,Arterial.markers,Venous.markers,CA.markers,choroid.markers,Interferon.markers)
#####画火山图#####
pdf("/md01/nieyg/project/BBB/plot/DEG_volcano_subtype_YO.pdf",width=5,height=5)
for (i in 1:7){
	mydata<-mylist[i]
	mydata<-as.data.frame(mydata)
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Young",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Old","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	p <-ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+
              ggtitle(levels(metadata$subtype)[i])+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Old'='red','Young'='deepskyblue','normal'='gray'));
    print(p)
}
dev.off()

m<-levels(metadata$subtype)
m[2]<-"CV";m[5]<-"CA";m[6]<-"choroid"
for (i in 1:7){
	mydata<-mylist[i]
	mydata<-as.data.frame(mydata)
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Young",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Old","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	assign(paste("Young_",m[i], sep = ""), rownames(mydata[which(mydata$Condition=="Young"),]) )
	assign(paste("Old_",m[i], sep = ""), rownames(mydata[which(mydata$Condition=="Old"),]) )
}
#######百分比########
Young_Capillary
Young_CV
Young_Arterial
Young_Venous
Young_CA
Young_choroid
Young_Interferon
Old_Capillary
Old_CV
Old_Arterial
Old_Venous
Old_CA
Old_choroid
Old_Interferon


subtype =rep(levels(metadata$subtype),2)
change =c(rep("Young",7),rep("Old",7))
num =c(5,8,3,9,2,5,2,26,19,19,35,10,15,16)
subtype_num = data.frame(subtype,change,num)

library(plyr)
subtype_per = ddply(subtype_num, "subtype", transform, percent = num / sum(num) * 100)
pdf("/md01/nieyg/project/BBB/plot/DEG_bar_percentage_YO.pdf")
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


alltermgene<-union(Old_Capillary,Old_CV)
alltermgene<-union(alltermgene,Old_Arterial)
alltermgene<-union(alltermgene,Old_Venous)
alltermgene<-union(alltermgene,Old_CA)
alltermgene<-union(alltermgene,Old_choroid)
alltermgene<-union(alltermgene,Old_Interferon)

all_Young_gene:
 [1] "Ivns1abp" "Hmcn1"    "Lars2"    "Gm42418"  "AY036118" "Macf1"   
 [7] "Sparc"    "Apcdd1"   "Ndnf"     "Ptn"      "Btg2"     "Utrn"    
[13] "Tpm1"     "Irf1"     "Dusp1"    "Csf1"     "Cxcl1"    "Sptbn1"  
[19] "Col4a1"   "Txnip"    "Ahnak"    "H1f0"     "Hbb-bs"  

all_Old_gene:
 [1] "Ttr"      "Ly6a"     "Acer2"    "Fmo2"     "Hbb-bs"   "Depp1"   
 [7] "Xdh"      "Pglyrp1"  "Tmem252"  "Plat"     "Kbtbd11"  "Selenop" 
[13] "Vwf"      "Fxyd5"    "H2-K1"    "Grrp1"    "Map3k6"   "Cp"      
[19] "Ier3"     "Stra6"    "Tgm2"     "Aldoa"    "Bnip3"    "Cd9"     
[25] "Plaur"    "Gadd45g"  "Mt1"      "Scgb3a1"  "Mt2"      "H2-Q7"   
[31] "H2-Q6"    "Lypd1"    "Mgp"      "Adamts9"  "S100a6"   "Clu"     
[37] "Socs3"    "Ifi27l2a" "Edn1"     "Hba-a1"   "Hba-a2"   "Igfbp7"  
[43] "Ly6c1"    "Ctla2a"   "Litaf"    "Srarp"    "Lcn2"     "Cd200"   
[49] "Podxl"    "Fas"      "Ifitm2"   "Ucp2"     "Xbp1"     "Gapdh"   
[55] "Lrg1"     "Cxcl12"   "Kansl1l"  "Plscr1"   "Galnt15"  "Slc16a1" 
[61] "H2-D1"    "Mat2a"    "Cyp26b1"  "Inhbb"    "Car4"     "Cdkn1a"  
[67] "Bsg"      "Por"      "Id4"      "H2-Q4"    "Slfn2"    "Bhlhe40" 

#####构建upsetR  matrix#########
meta<-matrix(data=NA,nrow=length(alltermgene),ncol=7);
colnames(meta)<-c("Capillary","CV","Arterial","Venous",
	    "CA","choroid_plexus","Interferon");
rownames(meta)<-alltermgene;

for ( i in 1:length(alltermgene))
{ 
  if (alltermgene[i] %in% Old_Capillary)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (alltermgene[i] %in% Old_CV)
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
  if (alltermgene[i] %in% Old_CA)
    meta[i,5] = 1
    else 
    meta[i,5] = 0
  if (alltermgene[i] %in% Old_choroid)
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

######在Young中的common gene在umap上的分布#######
o<-intersect(Young_Capillary,Young_CV)
o<-intersect(o,Young_Arterial)
pdf("/md01/nieyg/project/BBB/plot/DEG_overlap_Young_umap.pdf")
FeaturePlot(EC_cells.integrated, features = o[1])
FeaturePlot(EC_cells.integrated, features = o[2])
dev.off()
o<-intersect(Old_Capillary,Old_CV)
o<-intersect(o,Old_Arterial)
o<-intersect(o,Old_Venous)
o<-intersect(o,Old_CA)
o<-intersect(o,Old_choroid)
o<-intersect(o,Old_Interferon)
pdf("/md01/nieyg/project/BBB/plot/DEG_overlap_Old_umap.pdf")
FeaturePlot(EC_cells.integrated, features = o[1])
FeaturePlot(EC_cells.integrated, features = o[2])













###################
> Young_Capillary
[1] "Ivns1abp" "Hmcn1"    "Lars2"    "Gm42418"  "AY036118"
> Young_CV
[1] "Macf1"    "Sparc"    "Hmcn1"    "Apcdd1"   "Ivns1abp" "Ndnf"     "Gm42418" 
[8] "AY036118"
> Young_Arterial
[1] "Macf1"    "AY036118" "Gm42418" 
> Young_Venous
[1] "Ptn"   "Btg2"  "Utrn"  "Tpm1"  "Irf1"  "Ndnf"  "Dusp1" "Csf1"  "Cxcl1"
> Young_CA
[1] "Sptbn1" "Tpm1"  
> Young_choroid
[1] "Col4a1" "Txnip"  "Ahnak"  "H1f0"   "Hbb-bs"
> Young_Interferon
[1] "Sparc" "Btg2" 
> Old_Capillary
"Ttr"     "Ly6a"    "Acer2"   "Fmo2"    "Hbb-bs"  "Depp1"   "Xdh"    
"Pglyrp1" "Tmem252" "Plat"    "Kbtbd11" "Selenop" "H2-Q7"   "Vwf"    
"Fxyd5"   "H2-K1"   "Grrp1"   "Map3k6"  "Cp"      "Cd9"     "Ier3"   
"Stra6"   "Adamts9" "Tgm2"    "Mgp"     "Scgb3a1"
> Old_CV
  "Ttr"     "Acer2"   "Ly6a"    "Tmem252" "Aldoa"   "Depp1"   "Hbb-bs" 
  "Fmo2"    "Plat"    "H2-K1"   "Xdh"     "Bnip3"   "Cd9"     "Plaur"  
 "Vwf"     "Gadd45g" "Mt1"     "Scgb3a1" "Mt2"    
> Old_Arterial
"H2-Q7"    "Ttr"      "H2-Q6"    "Depp1"    "Hbb-bs"   "Lypd1"   
"Tmem252"  "Mgp"      "Gadd45g"  "Adamts9"  "S100a6"   "Clu"     
"Socs3"    "Pglyrp1"  "Ifi27l2a" "Edn1"     "Hba-a1"   "Mt1"     
"Hba-a2"  
> Old_Venous
"Igfbp7"  "Ly6a"    "Ttr"     "Ly6c1"   "Fxyd5"   "Scgb3a1" "Acer2"  
"Ctla2a"  "Mt2"     "H2-Q7"   "Grrp1"   "Litaf"   "Mt1"     "Srarp"  
"Aldoa"   "Lcn2"    "Xdh"     "Cd200"   "Podxl"   "Fas"     "Ifitm2" 
"Ucp2"    "Vwf"     "Xbp1"    "Adamts9" "H2-Q6"   "Cp"      "Gapdh"  
"Lrg1"    "Cxcl12"  "Kansl1l" "Fmo2"    "Plscr1"  "Galnt15" "Slc16a1"
> Old_CA
"Ttr"     "Acer2"   "H2-D1"   "Ly6a"    "Plat"    "H2-K1"   "Igfbp7" 
"Tmem252" "Mat2a"   "Pglyrp1"
> Old_choroid
"Cyp26b1" "H2-Q7"   "Ttr"     "H2-Q6"   "Acer2"   "Inhbb"   "Grrp1"  
"Car4"    "Cdkn1a"  "Gadd45g" "Depp1"   "Bsg"     "Por"     "Id4"    
> Old_Interferon
"H2-Q7"   "H2-Q6"   "Ly6a"    "Ttr"     "Depp1"   "H2-Q4"   "Hbb-bs" 
"Fmo2"    "Xdh"     "H2-K1"   "Acer2"   "Slfn2"   "Pglyrp1" "Bhlhe40"
"Tmem252" "Fxyd5"  

