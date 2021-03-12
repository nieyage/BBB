## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
#subset 取BBB亚群0-9,12
#subBBB_EC <-BBB_all[,which(BBB_all@meta.data$seurat_clusters==c(0:9,12))]
EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
DefaultAssay(EC_cells.integrated) <- "RNA" # Create dummy new assay to demo switching default assays
Idents(EC_cells.integrated)<-"orig.ident"
EC_cells <- subset(EC_cells.integrated,idents=c("Young","Old"))

DefaultAssay(EC_cells)<-"RNA"

######各个subtype的DEG 以及统计结果######
Idents(EC_cells)<-"celltype"
EC.markers <- FindAllMarkers(EC_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
EC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- EC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
1、ls()来查看当前所有对象名，对于每一个对象，可以通过object.size(x)来查看其占用内存的大小。

　　如果是因为当前对象占用内存过多，那么可以通过处理对象来获取更大的可用内存。一个很有用的方法是改变对象的存储模式，通过storage.mode（x）可以看到某个对象的存储模式，比如某个矩阵默认就是“double”的，如果这个矩阵的数值都是整数甚至0-1，完全没必要使用double来占用空间，可以使用storage.mode(x） <- "integer"将其改为整数型，可以看到该对象的大小会变为原来的一半。


2、object.size()看每个变量占多大内存。
3、memory.size()查看现在的work space的内存使用
4、memory.limit()查看系统规定的内存使用上限。如果现在的内存上限不够用，可以通过memory.limit(newLimit)更改到一个新的上限。注意，在32位的R中，封顶上限为4G，无法在一个程序上使用超过4G （数位上限）。这种时候，可以考虑使用64位的版本。

对于一些很大的但无用的中间变量，养成清理的习惯：

可以使用rm(object)删除变量，但是记住，rm后记得使用gc()做Garbage collection，否则内存是不会自动释放的，相当于你没做rm.

#####UP and Down genes#####
EC.markers <- FindAllMarkers(EC_cells, min.pct = 0.25, logfc.threshold = 0)
metadata<-EC_cells@meta.data
sample<-as.character(metadata$celltype)
time<-as.character(metadata$orig.ident)
sub<-as.character(metadata$seurat_clusters)
for (i in 1:21578){
if(sample[i]=="Capillary_1"&&time[i]=="Young"){sub[i]="Capillary_1_Young"}
if(sample[i]=="Capillary_2"&&time[i]=="Young"){sub[i]="Capillary_2_Young"}
if(sample[i]=="Capillary_3"&&time[i]=="Young"){sub[i]="Capillary_3_Young"}
if(sample[i]=="Capillary_4"&&time[i]=="Young"){sub[i]="Capillary_4_Young"}
if(sample[i]=="C_V_1"&&time[i]=="Young"){sub[i]="C_V_1_Young"}
if(sample[i]=="C_V_2"&&time[i]=="Young"){sub[i]="C_V_2_Young"}
if(sample[i]=="Arterial"&&time[i]=="Young"){sub[i]="Arterial_Young"}
if(sample[i]=="Venous"&&time[i]=="Young"){sub[i]="Venous_Young"}
if(sample[i]=="C_A"&&time[i]=="Young"){sub[i]="C_A_Young"}
if(sample[i]=="Choroid_plexus"&&time[i]=="Young"){sub[i]="Choroid_plexus_Young"}
if(sample[i]=="Interferon"&&time[i]=="Young"){sub[i]="Interferon_Young"}

if(sample[i]=="Capillary_1"&&time[i]=="Old"){sub[i]="Capillary_1_Aged"}
if(sample[i]=="Capillary_2"&&time[i]=="Old"){sub[i]="Capillary_2_Aged"}
if(sample[i]=="Capillary_3"&&time[i]=="Old"){sub[i]="Capillary_3_Aged"}
if(sample[i]=="Capillary_4"&&time[i]=="Old"){sub[i]="Capillary_4_Aged"}
if(sample[i]=="C_V_1"&&time[i]=="Old"){sub[i]="C_V_1_Aged"}
if(sample[i]=="C_V_2"&&time[i]=="Old"){sub[i]="C_V_2_Aged"}
if(sample[i]=="Arterial"&&time[i]=="Old"){sub[i]="Arterial_Aged"}
if(sample[i]=="Venous"&&time[i]=="Old"){sub[i]="Venous_Aged"}
if(sample[i]=="C_A"&&time[i]=="Old"){sub[i]="C_A_Aged"}
if(sample[i]=="Choroid_plexus"&&time[i]=="Old"){sub[i]="Choroid_plexus_Aged"}
if(sample[i]=="Interferon"&&time[i]=="Old"){sub[i]="Interferon_Aged"}
}
       Arterial_Aged       Arterial_Young             C_A_Aged 
                 683                  380                  776 
           C_A_Young           C_V_1_Aged          C_V_1_Young 
                 969                 1691                 1547 
          C_V_2_Aged          C_V_2_Young     Capillary_1_Aged 
                1155                 1294                 3919 
   Capillary_1_Young     Capillary_2_Aged    Capillary_2_Young 
                2238                 1545                   16 
    Capillary_3_Aged    Capillary_3_Young     Capillary_4_Aged 
                 117                  510                  357 
   Capillary_4_Young  Choroid_plexus_Aged Choroid_plexus_Young 
                 373                  217                  159 
     Interferon_Aged     Interferon_Young          Venous_Aged 
                1713                 1049                  597 
        Venous_Young 
                 273 


sub<-as.factor(sub)
names(sub)<-rownames(metadata)
EC_cells@active.ident <- sub
head(EC_cells@active.ident)
levels(metadata$celltype)
#########寻找各个subtype内部，YO之间的DEG#################
"Capillary"      "C-V"            "Arterial"       "Venous"        
"C-A"            "choroid plexus" "Interferon"   

Caplillary_1.markers <- FindMarkers(EC_cells, ident.1 = "Capillary_1_Young", ident.2 = "Capillary_1_Aged", min.pct = 0.25,logfc.threshold = 0)
Caplillary_2.markers <- FindMarkers(EC_cells, ident.1 = "Capillary_2_Young", ident.2 = "Capillary_2_Aged", min.pct = 0.25,logfc.threshold = 0)
Caplillary_3.markers <- FindMarkers(EC_cells, ident.1 = "Capillary_3_Young", ident.2 = "Capillary_3_Aged", min.pct = 0.25,logfc.threshold = 0)
Caplillary_4.markers <- FindMarkers(EC_cells, ident.1 = "Capillary_4_Young", ident.2 = "Capillary_4_Aged", min.pct = 0.25,logfc.threshold = 0)

CV_1.markers <- FindMarkers(EC_cells, ident.1 = "C_V_1_Young", ident.2 = "C_V_1_Aged", min.pct = 0.25,logfc.threshold = 0)
CV_2.markers <- FindMarkers(EC_cells, ident.1 = "C_V_2_Young", ident.2 = "C_V_2_Aged", min.pct = 0.25,logfc.threshold = 0)

CA.markers <- FindMarkers(EC_cells, ident.1 = "C_A_Young", ident.2 = "C_A_Aged", min.pct = 0.25,logfc.threshold = 0)
Arterial.markers <- FindMarkers(EC_cells, ident.1 = "Arterial_Young", ident.2 = "Arterial_Aged", min.pct = 0.25,logfc.threshold = 0)
Venous.markers <- FindMarkers(EC_cells, ident.1 = "Venous_Young", ident.2 = "Venous_Aged", min.pct = 0.25,logfc.threshold = 0)
Choroid.markers <- FindMarkers(EC_cells, ident.1 = "Choroid_plexus_Young", ident.2 = "Choroid_plexus_Aged", min.pct = 0.25,logfc.threshold = 0)
Interferon.markers <- FindMarkers(EC_cells, ident.1 = "Interferon_Young", ident.2 = "Interferon_Aged", min.pct = 0.25,logfc.threshold = 0)

mylist<-list(Caplillary_1.markers,Caplillary_2.markers,Caplillary_3.markers,Caplillary_4.markers,CV_1.markers,CV_2.markers,
	Arterial.markers,Venous.markers,CA.markers,Choroid.markers,Interferon.markers)
#####画火山图#####
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/DEG_volcano_subtype.pdf",width=5,height=5)
for (i in 1:11){
	mydata<-mylist[i]
	mydata<-as.data.frame(mydata)
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Young",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Aged","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	p <-ggplot(data=mydata, aes(x=avg_logFC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("ln fold change") + ylab("-log10 padj")+
              xlim(c(-2, 2)) +
              ggtitle(levels(metadata$subtype)[i])+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('Aged'='#E64B35','Young'='#00A087','normal'='gray'));
    print(p)
}
dev.off()

#######百分比########

m<-c("Caplillary_1","Caplillary_2","Caplillary_3","Caplillary_4","C_V_1","C_V_2","Arterial","Venous","C_A","Choroid_plexus","Interferon")
for (i in 1:11){
	mydata<-mylist[i]
	mydata<-as.data.frame(mydata)
	mydata$Condition=ifelse(mydata$avg_logFC>=0.5 & mydata$p_val_adj<=0.05,"Young",	
	ifelse(mydata$avg_logFC<=-0.5 & mydata$p_val_adj<=0.05,"Aged","normal"))##对满足不同条件的数据给不同的标记，放入Condition列
	assign(paste("Young_",m[i], sep = ""), rownames(mydata[which(mydata$Condition=="Young"),]) )
	assign(paste("Aged_",m[i], sep = ""), rownames(mydata[which(mydata$Condition=="Aged"),]) )
}
Young_Arterial
Young_C_A
Young_Caplillary_1
Young_Caplillary_2
Young_Caplillary_3
Young_Caplillary_4
Young_C_V_1
Young_C_V_2
Young_Venous
Young_Choroid_plexus
Young_Interferon

Aged_Arterial
Aged_C_A
Aged_Caplillary_1
Aged_Caplillary_2
Aged_Caplillary_3
Aged_Caplillary_4
Aged_C_V_1
Aged_C_V_2
Aged_Venous
Aged_Choroid_plexus
Aged_Interferon

 Young_Arterial
"Macf1"  ,  "AY036118" "Gm42418"  "Ccn3"    
 Young_C_A
] "Macf1"  ,  "Ivns1abp"
Young_Caplillary_1
] "Macf1"  ,  "Hmcn1"    "Ivns1abp" "Lars2"    "Gm42418"  "Ndnf"     "AY036118"
Young_Caplillary_2
] "Zfp655" ,"Sacs"   "Filip1"
Young_Caplillary_3
] "Gm42418",  "Lars2"    "Macf1"    "Atrx"     "Ivns1abp"
Young_Caplillary_4
] "Sptbn1" ,       "C130074G19Rik" "Tpm1"         
Young_C_V_1
] "Macf1"  ,"Sparc"  "Apcdd1"
Young_C_V_2
] "Macf1"  ,  "Hmcn1"    "Col4a1"   "Ivns1abp" "Apcdd1"   "Ndnf"     "Gm42418" 
Young_Venous
] "Ptn"    , "Btg2"    "Tpm1"    "Utrn"    "Zfp36l1" "Dleu2"   "Ndnf"   
] "Cxcl1"  , "Csf1"   
Young_Choroid_plexus
] "Txnip"  ,"Ahnak"  "Dab2"   "H1f0"   "Hbb-bs"
Young_Interferon
] "Macf1"  ,  "Btg2"     "Ivns1abp" "Gm42418"  "AY036118"
> 
> Aged_Arterial
 [1] "H2-Q7"    "Ttr"      "Lypd1"    "H2-Q6"    "Depp1"    "Gadd45g" 
 [7] "Hbb-bs"   "Pglyrp1"  "Mt1"      "S100a6"   "Tmem252"  "Ifi27l2a"
[13] "Socs3"    "Edn1"    
> Aged_C_A
 [1] "Ttr"     "H2-Q7"   "Depp1"   "Tmem252" "Acer2"   "Fxyd5"   "H2-K1"  
 [8] "Hbb-bs"  "Aldoa"   "Plat"    "Ldha"    "Cd9"     "Stmn2"   "Fmo2"   
[15] "Kbtbd11" "Pmp22"   "Litaf"   "Stra6"   "Xdh"     "Vwf"     "Lypd1"  
[22] "Mgp"     "Gadd45g" "Junb"    "Fos"    
> Aged_Caplillary_1
 [1] "Ly6a"    "Acer2"   "Ttr"     "Fmo2"    "Pglyrp1" "Xdh"     "Depp1"  
 [8] "Selenop" "Tmem252" "Hbb-bs"  "Fxyd5"   "Kbtbd11" "Plat"    "Vwf"    
[15] "H2-K1"   "Cp"      "Grrp1"   "Map3k6"  "Apold1"  "Stra6"   "Bnip3"  
[22] "Ier3"    "Tgm2"   
> Aged_Caplillary_2
character(0)
> Aged_Caplillary_3
 [1] "Ttr"     "Fmo2"    "Actg1"   "Apold1"  "Slc2a1"  "Ly6a"    "Ly6c1"  
 [8] "B2m"     "Actb"    "Ccn1"    "Timp3"   "Pmp22"   "Depp1"   "Klf2"   
[15] "Acer2"   "Mat2a"   "Plat"    "Cdkn1a"  "Clic4"   "Nr4a1"   "Fxyd5"  
[22] "Kdm6b"   "Pglyrp1" "Adamts1" "Bhlhe40" "Lrrc8a"  "Snrk"    "Crip1"  
[29] "Tinagl1" "Id1"     "Cp"      "Ndrg1"   "Sgk1"    "Rhob"    "H2-K1"  
[36] "Selenop" "Vwf"     "H2-T23"  "Ifitm2"  "Per1"    "Tmem252" "Rgcc"   
[43] "Itih5"   "Rhoc"    "Jund"   
> Aged_Caplillary_4
 [1] "Ttr"     "Acer2"   "Ly6a"    "H2-D1"   "Plat"    "H2-K1"   "Tmem252"
 [8] "Ctla2a"  "Mat2a"   "Aldoa"   "Malat1"  "Pglyrp1"
> Aged_C_V_1
 [1] "Ttr"     "Acer2"   "H2-Q7"   "Tmem252" "Aldoa"   "Fmo2"    "Plat"   
 [8] "Pglyrp1" "Xdh"     "Litaf"   "Bnip3"   "Depp1"   "Vwf"     "Plaur"  
[15] "Gadd45g" "Mt1"     "Mt2"     "Scgb3a1"
> Aged_C_V_2
 [1] "Ttr"     "Ly6a"    "Acer2"   "Fmo2"    "Depp1"   "Tmem252" "Aldoa"  
 [8] "H2-K1"   "Xdh"     "Plat"    "Vwf"     "Bnip3"   "Cebpd"   "Mt1"    
[15] "Scgb3a1" "Mt2"    
> Aged_Venous
 [1] "Igfbp7"  "Ttr"     "Ly6a"    "H2-Q7"   "Mt2"     "Ly6c1"   "Lcn2"   
 [8] "Fxyd5"   "Mt1"     "Acer2"   "Scgb3a1" "Cd200"   "Srarp"   "Grrp1"  
[15] "H2-Q6"   "Ifitm2"  "Podxl"   "Litaf"   "Ctla2a"  "Fas"     "Xdh"    
[22] "Lrg1"    "Vwf"     "Aldoa"   "Ucp2"    "Cp"      "Xbp1"    "Adamts9"
[29] "Osmr"    "Gapdh"   "Plscr1"  "Cxcl12"  "Sbno2"   "Il1r1"   "Kansl1l"
[36] "Gm47283" "Hif1a"   "Fmo2"    "Cd14"    "Slc16a1"
> Aged_Choroid_plexus
 [1] "Cyp26b1" "H2-Q7"   "Ttr"     "H2-Q6"   "Acer2"   "Inhbb"   "Grrp1"  
 [8] "Car4"    "Cdkn1a"  "Gadd45g" "Bsg"     "Depp1"   "Por"     "Id4"    
[15] "Apold1" 
> Aged_Interferon
 [1] "Ttr"     "H2-Q7"   "Depp1"   "Acer2"   "Hbb-bs"  "Tmem252" "H2-K1"  
 [8] "Pglyrp1" "Selenop" "Xdh"     "Fmo2"    "Kbtbd11" "Adamts9" "Bhlhe40"
[15] "Grrp1"   "Vwf"     "Cd9"     "Ier3"    "Stmn2"   "Stra6"   "Mgp"    
[22] "Hba-a1"  "Hba-a2" 


Young_Arterial<-c("Macf1"  ,  "AY036118", "Gm42418" , "Ccn3"    )
Young_C_A<-c( "Macf1"  ,  "Ivns1abp")
Young_Caplillary_1<-c("Macf1"  ,  "Hmcn1" ,   "Ivns1abp" ,"Lars2"   , "Gm42418"  ,"Ndnf"  ,   "AY036118")
Young_Caplillary_2<-c( "Zfp655" ,"Sacs"  , "Filip1")
Young_Caplillary_3<-c( "Gm42418",  "Lars2",    "Macf1"  ,  "Atrx"   ,  "Ivns1abp")
Young_Caplillary_4<-c("Sptbn1" ,       "C130074G19Rik" ,"Tpm1"   )
Young_C_V_1<-c( "Macf1"  ,"Sparc" , "Apcdd1")
Young_C_V_2<-c("Macf1"  ,  "Hmcn1"  ,  "Col4a1" ,  "Ivns1abp" ,"Apcdd1" ,  "Ndnf"  ,   "Gm42418")
Young_Venous<-c( "Ptn"    , "Btg2"  ,  "Tpm1",    "Utrn"   , "Zfp36l1", "Dleu2" ,  "Ndnf"  , "Cxcl1"  , "Csf1" )
Young_Choroid_plexus<-c("Txnip"  ,"Ahnak", "Dab2"  , "H1f0"  , "Hbb-bs")
Young_Interferon<-c("Macf1"  ,  "Btg2"   , "Ivns1abp" ,"Gm42418" , "AY036118")
 
Aged_Arterial<-c(
 "H2-Q7"  ,  "Ttr"    ,  "Lypd1",    "H2-Q6"  ,  "Depp1" ,   "Gadd45g" ,
 "Hbb-bs" ,  "Pglyrp1",  "Mt1"  ,    "S100a6"  , "Tmem252" , "Ifi27l2a",
 "Socs3"  ,  "Edn1"   ) 
Aged_C_A<-c(
 "Ttr"    , "H2-Q7"   ,"Depp1"  , "Tmem252", "Acer2"  , "Fxyd5" , "H2-K1"  ,
 "Hbb-bs" , "Aldoa"   ,"Plat"   , "Ldha"   , "Cd9"    , "Stmn2" , "Fmo2"   ,
 "Kbtbd11", "Pmp22"   ,"Litaf"  , "Stra6"  , "Xdh"    , "Vwf"   , "Lypd1"  ,
 "Mgp"    , "Gadd45g" ,"Junb"   , "Fos"  )
Aged_Caplillary_1<-c(
 "Ly6a"   , "Acer2"   ,"Ttr"    , "Fmo2"   , "Pglyrp1", "Xdh"   , "Depp1"  ,
 "Selenop", "Tmem252" ,"Hbb-bs" , "Fxyd5"  , "Kbtbd11", "Plat"  , "Vwf"    ,
 "H2-K1"  , "Cp"      ,"Grrp1"  , "Map3k6" , "Apold1" , "Stra6" , "Bnip3"  ,
 "Ier3"   , "Tgm2"   )
Aged_Caplillary_2<-c()

Aged_Caplillary_3<-c(
 "Ttr"     ,"Fmo2"   , "Actg1"  , "Apold1" , "Slc2a1" , "Ly6a"    ,"Ly6c1"  ,
 "B2m"     ,"Actb"   , "Ccn1"   , "Timp3"  , "Pmp22"  , "Depp1"   ,"Klf2"   ,
 "Acer2"   ,"Mat2a"  , "Plat"   , "Cdkn1a" , "Clic4"  , "Nr4a1"   ,"Fxyd5"  ,
 "Kdm6b"   ,"Pglyrp1", "Adamts1", "Bhlhe40", "Lrrc8a" , "Snrk"    ,"Crip1"  ,
 "Tinagl1" ,"Id1"    , "Cp"     , "Ndrg1"  , "Sgk1"   , "Rhob"    ,"H2-K1"  ,
 "Selenop" ,"Vwf"    , "H2-T23" , "Ifitm2" , "Per1"   , "Tmem252" ,"Rgcc"   ,
 "Itih5"   ,"Rhoc"   , "Jund")
Aged_Caplillary_4<-c(
 "Ttr"     ,"Acer2"  , "Ly6a"   , "H2-D1"  , "Plat"   , "H2-K1"   ,"Tmem252",
 "Ctla2a"  ,"Mat2a"  , "Aldoa"  , "Malat1" , "Pglyrp1")
Aged_C_V_1<-c(
 "Ttr"     ,"Acer2"  , "H2-Q7"  , "Tmem252", "Aldoa"  , "Fmo2"    ,"Plat"   ,
 "Pglyrp1" ,"Xdh"    , "Litaf"  , "Bnip3"  , "Depp1"  , "Vwf"     ,"Plaur"  ,
 "Gadd45g" ,"Mt1"    , "Mt2"    , "Scgb3a1")
Aged_C_V_2<-c(
 "Ttr"     ,"Ly6a"   , "Acer2"  , "Fmo2"   , "Depp1"  , "Tmem252" ,"Aldoa"  ,
 "H2-K1"   ,"Xdh"    , "Plat"   , "Vwf"    , "Bnip3"  , "Cebpd"   ,"Mt1"    ,
 "Scgb3a1" ,"Mt2" )
Aged_Venous<-c(
 "Igfbp7"  ,"Ttr"    , "Ly6a"   , "H2-Q7"  , "Mt2"    , "Ly6c1"   ,"Lcn2"   ,
 "Fxyd5"   ,"Mt1"    , "Acer2"  , "Scgb3a1", "Cd200"  , "Srarp"   ,"Grrp1"  ,
 "H2-Q6"   ,"Ifitm2" , "Podxl"  , "Litaf"  , "Ctla2a" , "Fas"     ,"Xdh"    ,
 "Lrg1"    ,"Vwf"    , "Aldoa"  , "Ucp2"   , "Cp"     , "Xbp1"    ,"Adamts9",
 "Osmr"    ,"Gapdh"  , "Plscr1" , "Cxcl12" , "Sbno2"  , "Il1r1"   ,"Kansl1l",
 "Gm47283" ,"Hif1a"  , "Fmo2"   , "Cd14"   , "Slc16a1")
Aged_Choroid_plexus<-c(
 "Cyp26b1" ,"H2-Q7"  , "Ttr"    , "H2-Q6"  , "Acer2"  , "Inhbb"   ,"Grrp1"  ,
 "Car4"    ,"Cdkn1a" , "Gadd45g", "Bsg"    , "Depp1"  , "Por"     ,"Id4"    ,
 "Apold1")
Aged_Interferon<-c(
 "Ttr"     ,"H2-Q7"  , "Depp1"  , "Acer2"  , "Hbb-bs" , "Tmem252" ,"H2-K1"  ,
 "Pglyrp1" ,"Selenop", "Xdh"    , "Fmo2"   , "Kbtbd11", "Adamts9" ,"Bhlhe40",
 "Grrp1"   ,"Vwf"    , "Cd9"    , "Ier3"   , "Stmn2"  , "Stra6"   ,"Mgp"    ,
 "Hba-a1"  ,"Hba-a2" )


m<-c("Arterial","C_A","Caplillary_1","Caplillary_2","Caplillary_3","Caplillary_4","C_V_1","C_V_2","Venous","Choroid_plexus","Interferon")

metadata$celltype <- factor(metadata$celltype, levels= m)

celltype =rep(levels(metadata$celltype),2)

change =c(rep("Young",11),rep("Aged",11))

num =c(length(Young_Arterial),length(Young_C_A),length(Young_Caplillary_1),length(Young_Caplillary_2),length(Young_Caplillary_3),length(Young_Caplillary_4),
	length(Young_C_V_1),length(Young_C_V_2),length(Young_Venous),length(Young_Choroid_plexus),length(Young_Interferon),length(Aged_Arterial),length(Aged_C_A),
	length(Aged_Caplillary_1),length(Aged_Caplillary_2),length(Aged_Caplillary_3),length(Aged_Caplillary_4),length(Aged_C_V_1),length(Aged_C_V_2),
	length(Aged_Venous),length(Aged_Choroid_plexus),length(Aged_Interferon))
celltype_num = data.frame(celltype,change,num)

         celltype change num
1        Arterial  Young   4
2             C_A  Young   2
3    Caplillary_1  Young   7
4    Caplillary_2  Young   3
5    Caplillary_3  Young   5
6    Caplillary_4  Young   3
7           C_V_1  Young   3
8           C_V_2  Young   7
9          Venous  Young   9
10 Choroid_plexus  Young   5
11     Interferon  Young   5
12       Arterial   Aged  14
13            C_A   Aged  25
14   Caplillary_1   Aged  23
15   Caplillary_2   Aged   0
16   Caplillary_3   Aged  45
17   Caplillary_4   Aged  12
18          C_V_1   Aged  18
19          C_V_2   Aged  16
20         Venous   Aged  40
21 Choroid_plexus   Aged  15
22     Interferon   Aged  23


library(plyr)
celltype_per = ddply(celltype_num, "celltype", transform, percent = num / sum(num) * 100)
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/DEG_bar_percentage_YA.pdf")
ggplot(celltype_per, aes(x = celltype, y = percent, fill = change)) +
  # 条形图函数：未将position参数显示设置为dodge，则绘制出的条形图为堆积型
  geom_bar(stat = "identity", colour = "black") +
  theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle = 90),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
  # 调色标尺
  scale_fill_brewer(palette = "Pastel1")

#######
celltype_per = arrange(celltype_num, celltype, change)
# 对不同Date分组内的数据进行累加求和
celltype_per = ddply(celltype_per, "celltype", transform, label_y = cumsum(num) - 0.5*num)
 
# 基函数
ggplot(celltype_per, aes(x = celltype, y = num, fill =change)) +
  # 条形图函数
  geom_bar(stat = "identity", colour = "black") +
  # 标签函数：paste和format方法对标签进行格式化
  geom_text(aes(y=label_y, label = paste(format(num, nsmall=0), "genes")), size = 3) +
  # 调色标尺
  scale_fill_brewer(palette = "Pastel1") +theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle = 90),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())
dev.off();

#####构建upsetR  matrix#########
Young_DEG<-c(Young_Caplillary_1,Young_Caplillary_2,Young_Caplillary_3,Young_Caplillary_4,Young_C_V_1,Young_C_V_2,
	Young_Arterial,Young_Venous,Young_C_A,Young_Choroid_plexus,Young_Interferon)
Young_DEG<-Young_DEG[!duplicated(Young_DEG)]
Aged_DEG<-c(	Aged_Caplillary_1,Aged_Caplillary_2,Aged_Caplillary_3,Aged_Caplillary_4,Aged_C_V_1,Aged_C_V_2,
	Aged_Arterial,Aged_Venous,Aged_C_A,Aged_Choroid_plexus,Aged_Interferon)
Aged_DEG<-Aged_DEG[!duplicated(Aged_DEG)]

> Young_DEG
 [1] "Macf1"         "Hmcn1"         "Ivns1abp"      "Lars2"        
 [5] "Gm42418"       "Ndnf"          "AY036118"      "Zfp655"       
 [9] "Sacs"          "Filip1"        "Atrx"          "Sptbn1"       
[13] "C130074G19Rik" "Tpm1"          "Sparc"         "Apcdd1"       
[17] "Col4a1"        "Ccn3"          "Ptn"           "Btg2"         
[21] "Utrn"          "Zfp36l1"       "Dleu2"         "Cxcl1"        
[25] "Csf1"          "Txnip"         "Ahnak"         "Dab2"         
[29] "H1f0"          "Hbb-bs"       
> Aged_DEG
  [1] "Ly6a"     "Acer2"    "Ttr"      "Fmo2"     "Pglyrp1"  "Xdh"     
  [7] "Depp1"    "Selenop"  "Tmem252"  "Hbb-bs"   "Fxyd5"    "Kbtbd11" 
 [13] "Plat"     "Vwf"      "H2-K1"    "Cp"       "Grrp1"    "Map3k6"  
 [19] "Apold1"   "Stra6"    "Bnip3"    "Ier3"     "Tgm2"     "Actg1"   
 [25] "Slc2a1"   "Ly6c1"    "B2m"      "Actb"     "Ccn1"     "Timp3"   
 [31] "Pmp22"    "Klf2"     "Mat2a"    "Cdkn1a"   "Clic4"    "Nr4a1"   
 [37] "Kdm6b"    "Adamts1"  "Bhlhe40"  "Lrrc8a"   "Snrk"     "Crip1"   
 [43] "Tinagl1"  "Id1"      "Ndrg1"    "Sgk1"     "Rhob"     "H2-T23"  
 [49] "Ifitm2"   "Per1"     "Rgcc"     "Itih5"    "Rhoc"     "Jund"    
 [55] "H2-D1"    "Ctla2a"   "Aldoa"    "Malat1"   "H2-Q7"    "Litaf"   
 [61] "Plaur"    "Gadd45g"  "Mt1"      "Mt2"      "Scgb3a1"  "Cebpd"   
 [67] "Lypd1"    "H2-Q6"    "S100a6"   "Ifi27l2a" "Socs3"    "Edn1"    
 [73] "Igfbp7"   "Lcn2"     "Cd200"    "Srarp"    "Podxl"    "Fas"     
 [79] "Lrg1"     "Ucp2"     "Xbp1"     "Adamts9"  "Osmr"     "Gapdh"   
 [85] "Plscr1"   "Cxcl12"   "Sbno2"    "Il1r1"    "Kansl1l"  "Gm47283" 
 [91] "Hif1a"    "Cd14"     "Slc16a1"  "Ldha"     "Cd9"      "Stmn2"   
 [97] "Mgp"      "Junb"     "Fos"      "Cyp26b1"  "Inhbb"    "Car4"    
[103] "Bsg"      "Por"      "Id4"      "Hba-a1"   "Hba-a2"  

meta<-matrix(data=NA,nrow=length(Young_DEG),ncol=11);
colnames(meta)<-c("Arterial","C_A","Caplillary_1","Caplillary_2","Caplillary_3","Caplillary_4","C_V_1","C_V_2","Venous","Choroid_plexus","Interferon");
rownames(meta)<-Young_DEG;

for ( i in 1:length(Young_DEG))
{ 
  if (Young_DEG[i] %in% Young_Arterial)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (Young_DEG[i] %in% Young_C_A)
    meta[i,2] = 1
    else 
    meta[i,2] = 0  
  if (Young_DEG[i] %in% Young_Caplillary_1)
    meta[i,3] = 1
    else 
    meta[i,3] = 0
  if (Young_DEG[i] %in% Young_Caplillary_2)
    meta[i,4] = 1
    else 
    meta[i,4] = 0
  if (Young_DEG[i] %in% Young_Caplillary_3)
    meta[i,5] = 1
    else 
    meta[i,5] = 0
  if (Young_DEG[i] %in% Young_Caplillary_4)
    meta[i,6] = 1
    else 
    meta[i,6] = 0

  if (Young_DEG[i] %in% Young_C_V_1)
    meta[i,7] = 1
    else 
    meta[i,7] = 0
  if (Young_DEG[i] %in% Young_C_V_2)
    meta[i,8] = 1
    else 
    meta[i,8] = 0
  if (Young_DEG[i] %in% Young_Venous)
    meta[i,9] = 1
    else 
    meta[i,9] = 0
  if (Young_DEG[i] %in% Young_Choroid_plexus)
    meta[i,10] = 1
    else 
    meta[i,10] = 0
   if (Young_DEG[i] %in% Young_Interferon)
    meta[i,11] = 1
    else 
    meta[i,11] = 0
}

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/DEG_upsetR_Ageing_up.pdf")
#meta = data.frame(meta)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
library(UpSetR)
m<-as.data.frame(meta)
c<-cbind(Young_DEG,m)
upset(c, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      nsets = 7, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Aging UP gene genes",
      sets.x.label = "Aging UP gene", text.scale = c(2, 2, 0.8, 2, 2, 2))

mylist<-list(data.frame(Young_Arterial),data.frame(Young_C_A),data.frame(Young_Caplillary_1),data.frame(Young_Caplillary_2),
	data.frame(Young_Caplillary_3),	data.frame(Young_Caplillary_4),data.frame(Young_C_V_1),data.frame(Young_C_V_2),
	data.frame(Young_Venous),data.frame(Young_Choroid_plexus),data.frame(Young_Interferon))
write.csv(do.call(rbind.fill,mylist),"Young_DEG.csv")
mylevels<-c("Arterial","C_A","Caplillary_1","Caplillary_2","Caplillary_3","Caplillary_4","C_V_1","C_V_2","Venous","Choroid_plexus","Interferon");
Idents(EC_cells.integrated) <- factor(Idents(EC_cells.integrated), levels= my_levels)
EC_cells.integrated@meta.data$orig.ident<-factor(EC_cells.integrated@meta.data$orig.ident,levels=c("Young","Middle1","Aged"))



StackedVlnPlot(EC_cells.integrated, Aged_DEG[10:12], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[13:15], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[16:18], pt.size=0, cols=my36colors)

for (j in 1:10){
	i=i+3;
	StackedVlnPlot(EC_cells.integrated, Aged_DEG[i:i+2], pt.size=0, cols=my36colors)
}

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Young_DEG_Violin.pdf")
StackedVlnPlot(EC_cells.integrated, Young_DEG[1:3], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[4:6], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[7:9], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[10:12], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[13:15], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[16:18], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[19:21], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[22:24], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[25:27], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[28:30], pt.size=0, cols=my36colors)

dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Young_DEG_Violin.pdf")

StackedVlnPlot(EC_cells.integrated, Young_DEG[1:4], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[5:8], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[9:12], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[13:16], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[17:20], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[21:24], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[25:28], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Young_DEG[29:30], pt.size=0, cols=my36colors)

dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/Aged_DEG_Violin.pdf")

StackedVlnPlot(EC_cells.integrated, Aged_DEG[1:4], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[5:8], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[9:12], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[13:16], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[17:20], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[21:24], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[25:28], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[29:32], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[33:36], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[37:40], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[41:44], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[45:48], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[49:52], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[53:56], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[57:60], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[61:64], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[65:68], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[69:72], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[73:76], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[77:80], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[81:84], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[85:88], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[89:92], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[93:96], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[97:100], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[101:104], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Aged_DEG[105:107], pt.size=0, cols=my36colors)
dev.off();


