##############BBB pseudotime######
#####monocle downsample###########
## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(ArchR)
library(monocle)
EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
EC_five<-subset(EC_cells.integrated, idents = c( "Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4","C_V_1","C_V_2","Venous"))
###########pesudotime rm                                                                EC_five@meta.data$subtype<-Idents(EC_five)
DefaultAssay(EC_five)<-"RNA"
set.seed(2)
pdf("monocle_markergene-23.pdf")

all<-1:26063
down<-sample(all,3000)
EC_five@meta.data$cellorder<-1:26063
down_monocle<- EC_five[ , which(EC_five@meta.data$cellorder %in% down)]
EC_five<-down_monocle
EC_five<-readRDS("EC_five.rds")
new.cluster.ids <- c( "Arterial","C_A","Capillary","Capillary","Capillary","Capillary","C_V","C_V","Venous")
names(new.cluster.ids) <- levels(EC_five)
EC_five <- RenameIdents(EC_five, new.cluster.ids)
EC_five$celltype<-EC_five@active.ident

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(EC_five@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EC_five@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,phenoData = pd,	featureData = fd,lowerDetectionLimit = 0.5,
	expressionFamily = negbinomial.size())

	 ####使用RNA count
#当数据的表达量类型是FPKM值，分布类型设为tobit()使用expressionFamily=VGAM:::tobit(Lower=0.1)。
#当数据为UMIs, Transcript counts时，数据分布需要设为负二项分布，即negbinomial.size()。
#由于做差异分析通常用到mRNA couts，常常按照官方的方法转成分子数，使用rpc_matrix <- relative2abs(HSMM, method = "num_genes")。

# 归一化 
 monocle_cds<- estimateSizeFactors(monocle_cds)
 monocle_cds<- estimateDispersions(monocle_cds)
###Filtering low-quality cells
 monocle_cds <- monocle::detectGenes(monocle_cds, min_expr = 3 )
head(featureData(monocle_cds)@data)
expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
                                     num_cells_expressed >= 10))
features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
  "Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
  "Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",####C
  "Tfrc", "Car4", "Itm2a","Chn2",#####C-V
  "Lcn2","Slc38a5","Nr2f2", "Sox12",##Venous
  "Vcam1", "Vwf" ####A/V
  )

disp_table <- monocle::dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

f<-intersect(intersect(unsup_clustering_genes$gene_id,features),intersect(expressed_genes,features))
monocle_cds <- monocle::setOrderingFilter(monocle_cds, f)
#pdf("monocle_test.pdf")
#plot_ordering_genes(monocle_cds)

#Trajectory step 2: reduce data dimensionality  
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                            method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)

pdf("BBB-trajectory.pdf")
plot_cell_trajectory(monocle_cds,  color_by = "celltype",cell_size=1)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

plot_cell_trajectory(monocle_cds,  color_by = "celltype",cell_size=0.5) +
    facet_wrap(~celltype, nrow =4)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

plot_cell_trajectory(monocle_cds,  color_by = "State",cell_size=1)
plot_cell_trajectory(monocle_cds,  color_by = "Pseudotime",cell_size=1)


plot_cell_trajectory(monocle_cds,  color_by = "seurat_clusters",cell_size=1)

The23genes <- row.names(subset(fData(monocle_cds),
                               gene_short_name %in% c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
                                       'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd",
                                       "Tcf7","Lef1","Sox17","Erg")))
plot_genes_branched_pseudotime(monocle_cds[The23genes,],
                               branch_point = 1,ncol=2,
                               color_by = "celltype")

The23genes<-c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
                                       'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd",
                                       "Tcf7","Lef1","Sox17","Erg")
pdf("BBB-genes-TFs-zonation.pdf")
cds_subset <- monocle_cds[The23genes,]
plot_genes_in_pseudotime(cds_subset[1:9], nrow=3,ncol=3,color_by = "celltype",cell_size=0.50)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

plot_genes_in_pseudotime(cds_subset[10:18], nrow=3,ncol=3,color_by = "celltype",cell_size=0.50)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

plot_genes_in_pseudotime(cds_subset[19:27], nrow=3,ncol=3,color_by = "celltype",cell_size=0.50)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")
dev.off();

all<-rownames(monocle_cds@phenoData@data)
b<-all[which(monocle_cds@phenoData@data$State %in% c(2,4,7))]
a<-setdiff(all,b)
m<-t(monocle_cds@reducedDimS)
m<-m[a,]
m<-as.data.frame(m)

md<-EC_five@meta.data
md<-md[a,]
m$celltype<-md$celltype
colnames(m)<-c("Component1","Component2","celltype")
data<-monocle_cds@phenoData@data
data<-data[a,]
mdata<-cbind(m,data)
mdata<-mdata[,c(1:3,15:16)]

pdf("linear-BBB-traj-span0.15.pdf")
ggplot(data=mdata, aes(x=Component1, y=Component2,)) +geom_point(alpha=1, size=1,aes(color=celltype))  +  theme_bw()+
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
      scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264","Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000"))+
       stat_smooth(color="black",se=FALSE,span=0.15)

ggplot(data=mdata, aes(x=Component1, y=Component2,)) +geom_point(alpha=1, size=1,aes(color=Pseudotime))  +  theme_bw()+
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
scale_fill_gradientn(colors = c("#54AEF2", "#326997", "#152F48"))+
       stat_smooth(color="black",se=FALSE,span=0.15)

ggplot(data=mdata, aes(x=Component1, y=Component2,)) +geom_point(alpha=1, size=1,aes(color=State))  +  theme_bw()+
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
       stat_smooth(color="black",se=FALSE,span=0.15)
dev.off();
######The BBBgenes zonation########

count<-t(EC_five@assays$RNA[features,a])
count<-as.data.frame(count)

mcount<-cbind(mdata,count)
pdf("Themarkergenes-zonation.pdf")
i=1;
for(i in 1:31){
p<-ggplot(data=mcount, aes(x=Pseudotime, y=mcount[,i+5],)) +geom_point(alpha=1, size=1,aes(color=celltype))  +  
theme_bw()+ylab("expression") +ggtitle(features[i]) +
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
      scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264","Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000"))+
       stat_smooth(color="black",se=FALSE);
       print(p)
     }
dev.off();

features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
    'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
count<-t(EC_five@assays$RNA[features,a])
count<-as.data.frame(count)

mcount<-cbind(mdata,count)
pdf("TheBBBgenes-zonation.pdf")
i=1;
for(i in 1:23){
p<-ggplot(data=mcount, aes(x=Pseudotime, y=mcount[,i+5],)) +geom_point(alpha=1, size=1,aes(color=celltype))  +  
theme_bw()+ylab("expression") +ggtitle(features[i]) +
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
      scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264","Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000"))+
       stat_smooth(color="black",se=FALSE);
       print(p)
     }
dev.off();





The23genes<-c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
                                       'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd",
                                       "Tcf7","Lef1","Sox17","Erg")
pdf("BBB-genes-TFs-zonation.pdf")
cds_subset <- monocle_cds[The23genes,]
plot_genes_in_pseudotime(cds_subset[1:9], relative_expr = FALSE,nrow=3,ncol=3,color_by = "celltype",cell_size=0.50)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

plot_genes_in_pseudotime(cds_subset[10:18], relative_expr = FALSE,nrow=3,ncol=3,color_by = "celltype",cell_size=0.50)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

plot_genes_in_pseudotime(cds_subset[19:27], relative_expr = FALSE,nrow=3,ncol=3,color_by = "celltype",cell_size=0.50)+ 
 scale_color_manual(breaks = c("Arterial","C_A","Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")
dev.off();


pdf("BBB-genes-TFs-groupbycelltype-colorbyage.pdf")
for(i in 1:9){
  j<-3*i-2
  k<-3*i
Thegene<-The23genes[j:k]
print(Thegene)
cds_subset <- monocle_cds[Thegene,]
levels(cds_subset$orig.ident)[3]<-"Aged"
levels(cds_subset$orig.ident)[2]<-"Middle"
cds_subset$orig.ident<-factor(cds_subset$orig.ident,levels=c("Young","Middle","Aged"))
p<-plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow= 3,ncol=1,plot_trend = TRUE,cell_size=0.50)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
print(p)
}
dev.off();

pdf("BBB-genes-TFs-patternheatmap.pdf")
cds_subset <- monocle_cds[The23genes,]
plot_pseudotime_heatmap(cds_subset,
                #num_clusters = 3,
                cores = 1,
                show_rownames = T)
dev.off();

pdf("Aging-up-gene-groupbycelltype-colorbyage.pdf")
for(i in 1:36){
  j<-3*i-2
  k<-3*i
Thegene<-Aged_DEG[j:k]
print(Thegene)
cds_subset <- monocle_cds[Thegene,]
levels(cds_subset$orig.ident)[3]<-"Aged"
levels(cds_subset$orig.ident)[2]<-"Middle"
cds_subset$orig.ident<-factor(cds_subset$orig.ident,levels=c("Young","Middle","Aged"))
p<-plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow= 3,ncol=1,plot_trend = TRUE,cell_size=0.50)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
print(p)
}
dev.off();


pdf("Aging-down-gene-groupbycelltype-colorbyage.pdf")
for(i in 1:10){
  j<-3*i-2
  k<-3*i
Thegene<-Young_DEG[j:k]
print(Thegene)
cds_subset <- monocle_cds[Thegene,]
levels(cds_subset$orig.ident)[3]<-"Aged"
levels(cds_subset$orig.ident)[2]<-"Middle"
cds_subset$orig.ident<-factor(cds_subset$orig.ident,levels=c("Young","Middle","Aged"))
p<-plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow= 3,ncol=1,plot_trend = TRUE,cell_size=0.50)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
print(p)
}
dev.off();


all<-rownames(monocle_cds@phenoData@data)
Young<-all[which(monocle_cds@phenoData@data$orig.ident=="Young")]
The23genes
count<-t(EC_five@assays$RNA[The23genes,Young])
count<-as.data.frame(count)
data<-monocle_cds@phenoData@data
data<-data[Young,]
Youngcount<-cbind(data,count)

Middle<-all[which(monocle_cds@phenoData@data$orig.ident=="Middle1")]
count<-t(EC_five@assays$RNA[The23genes,Middle])
count<-as.data.frame(count)
data<-monocle_cds@phenoData@data
data<-data[Middle,]
Middlecount<-cbind(data,count)

Aged<-all[which(monocle_cds@phenoData@data$orig.ident=="Old")]
count<-t(EC_five@assays$RNA[The23genes,Aged])
count<-as.data.frame(count)
data<-monocle_cds@phenoData@data
data<-data[Aged,]
Agedcount<-cbind(data,count)


require(grid)
pdf("TheBBBgenes-zonation.pdf")
for(i in 1:27){
p<-ggplot(data=Youngcount, aes(x=Pseudotime, y=Youngcount[,i+13],)) +
geom_point(alpha=1, size=1,aes(color=celltype))  +  
theme_bw()+ylab("expression") +ggtitle(features[i]) +
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
      scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264",
        "Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000"))+
       stat_smooth(color="#009F86",se=FALSE);

q<-ggplot(data=Middlecount, aes(x=Pseudotime, y=Middlecount[,i+13],)) +
geom_point(alpha=1, size=1,aes(color=celltype))  +  
theme_bw()+ylab("expression") +ggtitle(features[i]) +
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
      scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264",
        "Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000"))+
       stat_smooth(color="#FF9F40",se=FALSE);
k<-ggplot(data=Agedcount, aes(x=Pseudotime, y=Agedcount[,i+13],)) +
geom_point(alpha=1, size=1,aes(color=celltype))  +  
theme_bw()+ylab("expression") +ggtitle(features[i]) +
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
      scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264",
        "Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000"))+
       stat_smooth(color="#E64A35",se=FALSE);
####排版图片#####
grid.newpage()  ###新建图表版面
pushViewport(viewport(layout = grid.layout(2,2))) ####将版面分成2*2矩阵
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
print(p, vp = vplayout(1,1))
print(q, vp = vplayout(1,2))
print(k, vp = vplayout(2,1))
     }
dev.off();

par(new=TRUE) 


all<-rownames(monocle_cds@phenoData@data)
count<-t(EC_five@assays$RNA[The23genes,])
count<-as.data.frame(count)
data<-monocle_cds@phenoData@data
data<-data[all,]
count<-cbind(data,count)
levels(count$orig.ident)[3]<-"Aged"
levels(count$orig.ident)[2]<-"Middle"
aes(color=celltype)
pdf("TheBBBgenes-zonation.pdf")
for(i in 14:40){
p<-ggplot(data=count, aes(x=Pseudotime, y=count[,i],group=orig.ident,col=orig.ident)) +
geom_point(alpha=1, size=1,aes(color=celltype))  +  
     scale_color_manual(values=c("Arterial"="#206A5D","C_A"="#81B264",
       "Capillary"="#FFCC29","C_V"="#F58634","Venous"="#BE0000","Young"="#009F86","Middle"="#FF9F40","Aged"="#E64A35"))+
theme_bw()+ylab("expression") +ggtitle(The23genes[i]) +
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+ 
geom_smooth(se=FALSE);
      print(p)
     }
dev.off();

