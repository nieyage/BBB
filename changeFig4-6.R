###change Fig4-6#####

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(ArchR)
library(patchwork)
library(magrittr)

EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
DefaultAssay(EC_cells.integrated) <- "RNA" 

BBBgene<-c("Tfrc","Bsg","Apoe")
library(ggplot2)
cols=c("Arterial"="#206A5D","C_A"="#81B264","Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29",
    "C_V_1"="#F58634","C_V_2"="#F58634","Venous"="#BE0000","Choroid_plexus"="#31326F","Interferon"="#93329E")
level<-c("Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4",
    "C_V_1","C_V_2","Venous","Interferon","Choroid_plexus")
EC_cells.integrated$celltype<-factor(EC_cells.integrated$celltype,levels=level)
EC_cells.integrated@active.ident<-EC_cells.integrated$celltype

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VllnPlot(obj, features = feature,assay="RNA",group.by="celltype", col=cols,pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +geom_boxplot(width=0.2,outlier.size=0,position=position_dodge(width=0.9))+
               #+geom_boxplot(width=0.2)+
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               #axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               #axis.ticks.y = element_line(),
               #axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
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

pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/BBB_gene_Violin.pdf")
StackedVlnPlot(EC_cells.integrated, BBBgene, pt.size=0)
dev.off();

#####downsample celltype cell number for 23 genes heatmap ############
#####The 23 BBB genes######
features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
    'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
EC_cells.integrated <- ScaleData(EC_cells.integrated, features =rownames(EC_cells.integrated))
pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/23genes_heatmap_forsubtype-nodownsample.pdf")

DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = fill.colors)
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#FA00FA", "black", "#F0F000"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("navy", "black", "firebrick3"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#1873CC", "#00CC00", "#FFFF00"))
DoHeatmap(EC_cells.integrated,features,group.by = "celltype",group.colors=cols,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#352A86", "#2DB7A3", "#F8FA0D"))

dev.off();




cells =
table(EC_cells.integrated$celltype)

   Capillary_1          C_V_1     Interferon          C_V_2            C_A 
          7860           5164           3641           3633           2603 
   Capillary_2       Arterial    Capillary_3    Capillary_4         Venous 
          1568           1525           1520           1142           1048 
Choroid_plexus 
           483 
Capillary_1<-which(EC_cells.integrated$celltype =="Capillary_1")
Capillary_2<-which(EC_cells.integrated$celltype =="Capillary_2")
Capillary_3<-which(EC_cells.integrated$celltype =="Capillary_3")
Capillary_4<-which(EC_cells.integrated$celltype =="Capillary_4")
Arterial<-which(EC_cells.integrated$celltype =="Arterial")
C_A     <-which(EC_cells.integrated$celltype =="C_A")
Venous  <-which(EC_cells.integrated$celltype =="Venous")
C_V_1   <-which(EC_cells.integrated$celltype =="C_V_1")
C_V_2   <-which(EC_cells.integrated$celltype =="C_V_2")
Interferon    <-which(EC_cells.integrated$celltype =="Interferon")
Choroid_plexus<-which(EC_cells.integrated$celltype =="Choroid_plexus")

Capillary_1<-sample(Capillary_1,438)
Capillary_2<-sample(Capillary_2,438)
Capillary_3<-sample(Capillary_3,438)
Capillary_4<-sample(Capillary_4,438)
Arterial<-sample(Arterial,438)
C_A     <-sample(C_A     ,438)
Venous  <-sample(Venous  ,438)
C_V_1   <-sample(C_V_1   ,438)
C_V_2   <-sample(C_V_2   ,438)
Interferon<-sample(Interferon,438)
down<-c(Capillary_1,Capillary_2,Capillary_3,Capillary_4,Arterial,C_A     ,
Venous  ,C_V_1   ,C_V_2   ,Interferon,Choroid_plexus)

EC_cells.integrated@meta.data$cellorder<-1:30187
DefaultAssay(EC_cells.integrated)<-"RNA"
EC_cells<- EC_cells.integrated[ , which(EC_cells.integrated@meta.data$cellorder %in% down)]
DefaultAssay(EC_cells)<-"RNA"

EC_cells <- ScaleData(EC_cells features =rownames(EC_cells))
pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/23genes_heatmap_forsubtype-downsample.pdf")

DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = fill.colors)
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#FA00FA", "black", "#F0F000"))
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "black", "red"))
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("navy", "black", "firebrick3"))
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#1873CC", "#00CC00", "#FFFF00"))
DoHeatmap(EC_cells,features,group.by = "celltype",group.colors=cols,size=4,disp.min=-1,disp.max=1)+scale_fill_gradientn(colors = c("#352A86", "#2DB7A3", "#F8FA0D"))

dev.off();





##
level<-c("Young","Middle1","Old")
EC_cells.integrated$orig.ident<-factor(EC_cells.integrated$orig.ident,levels=level)
cols<-c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35")

pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/Umap-sample.pdf.pdf")
DimPlot(EC_cells.integrated, reduction = "umap", group.by = "orig.ident",cols=cols)



EC_five<-subset(EC_cells.integrated, idents = c( "Capillary_2","Capillary_3"))
data <- as(as.matrix(EC_five@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EC_five@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,phenoData = pd,	featureData = fd,lowerDetectionLimit = 0.5,
	expressionFamily = negbinomial.size())


 monocle_cds<- estimateSizeFactors(monocle_cds)
 monocle_cds<- estimateDispersions(monocle_cds)
###Filtering low-quality cells
 monocle_cds <- monocle::detectGenes(monocle_cds, min_expr = 3 )
head(featureData(monocle_cds)@data)
expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
                                     num_cells_expressed >= 10))

disp_table <- monocle::dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

f<-intersect(unsup_clustering_genes$gene_id,expressed_genes)
monocle_cds <- monocle::setOrderingFilter(monocle_cds, f)
#pdf("monocle_test.pdf")
#plot_ordering_genes(monocle_cds)

#Trajectory step 2: reduce data dimensionality  
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                            method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)



library(Hmisc)
plot_genes_jitter<-function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0, 
    nrow = NULL, ncol = 2, panel_order = NULL, color_by = NULL, 
    plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE) 
{
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
        "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- exprs(cds_subset)
        cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- cds_exprs$gene_short_name
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
            levels = panel_order)
    }
    q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q 
    }
    else {
        q <- q 
    }
    if (plot_trend == TRUE) {
        q <- q + stat_summary(aes_string(color = color_by), fun.data = "mean_cl_boot", 
            size = 0.35)
        q <- q + stat_summary(aes_string(x = grouping, y = "expression", 
            color = color_by, group = color_by), fun.data = "mean_cl_boot", 
            size = 0.35, geom = "line")
    }
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
        ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    q <- q + ylab("Expression") + xlab(grouping)
    q <- q+theme_bw()+ theme(panel.grid=element_blank())
    q
}


data=mydata[mydata$Condition %in% c("Capillary_2","Capillary_3"),]
gene<-data$gene


levels(monocle_cds$orig.ident)[3]<-"Aged"
levels(monocle_cds$orig.ident)[2]<-"Middle"
monocle_cds$orig.ident<-factor(monocle_cds$orig.ident,levels=c("Young","Middle","Aged"))
pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/C3vsC3-gene-groupbycelltype-colorbyage.pdf")
cds_subset <- monocle_cds[gene[1:8],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",nrow= 4,ncol=2,cell_size =0, plot_trend = TRUE)+
scale_color_manual(breaks = c("Capillary_2","Capillary_3"), values=c("#fed049","#007580")) + theme(legend.position = "right")
cds_subset <- monocle_cds[gene[9:16],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",nrow= 4,ncol=2,cell_size =0, plot_trend = TRUE)+
scale_color_manual(breaks = c("Capillary_2","Capillary_3"), values=c("#fed049","#007580")) + theme(legend.position = "right")
cds_subset <- monocle_cds[gene[17:24],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",nrow= 4,ncol=2,cell_size =0, plot_trend = TRUE)+
scale_color_manual(breaks = c("Capillary_2","Capillary_3"), values=c("#fed049","#007580")) + theme(legend.position = "right")
cds_subset <- monocle_cds[gene[25:32],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",nrow= 4,ncol=2,cell_size =0, plot_trend = TRUE)+
scale_color_manual(breaks = c("Capillary_2","Capillary_3"), values=c("#fed049","#007580")) + theme(legend.position = "right")
cds_subset <- monocle_cds[gene[30:37],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",nrow= 4,ncol=2,cell_size =0, plot_trend = TRUE)+
scale_color_manual(breaks = c("Capillary_2","Capillary_3"), values=c("#fed049","#007580")) + theme(legend.position = "right")
dev.off()




EC_five<-subset(EC_cells.integrated, idents = c( "Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4","C_V_1","C_V_2","Venous"))
new.cluster.ids <- c( "Arterial","C_A","Capillary","Capillary","Capillary","Capillary","C_V","C_V","Venous")
names(new.cluster.ids) <- levels(EC_five)
EC_five <- RenameIdents(EC_five, new.cluster.ids)
EC_five$celltype<-EC_five@active.ident
data <- as(as.matrix(EC_five@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EC_five@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,phenoData = pd,	featureData = fd,lowerDetectionLimit = 0.5,
	expressionFamily = negbinomial.size())

 monocle_cds<- estimateSizeFactors(monocle_cds)
 monocle_cds<- estimateDispersions(monocle_cds)
###Filtering low-quality cells
# monocle_cds <- monocle::detectGenes(monocle_cds, min_expr = 3 )
#head(featureData(monocle_cds)@data)
#expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
#                                     num_cells_expressed >= 10))
#features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
#  "Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
#  "Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",####C
#  "Tfrc", "Car4", "Itm2a","Chn2",#####C-V
#  "Lcn2","Slc38a5","Nr2f2", "Sox12",##Venous
#  "Vcam1", "Vwf" ####A/V
#  )
#
#disp_table <- monocle::dispersionTable(monocle_cds)
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
#
#f<-intersect(union(unsup_clustering_genes$gene_id,features),union(expressed_genes,features))
#monocle_cds <- monocle::setOrderingFilter(monocle_cds, f)
#monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
#                            method = 'DDRTree')
#monocle_cds <- orderCells(monocle_cds)
monocle_cds$orig.ident<-factor(monocle_cds$orig.ident)
levels(monocle_cds$orig.ident)[2]<-"Aged"
levels(monocle_cds$orig.ident)[1]<-"Middle"
monocle_cds$orig.ident<-factor(monocle_cds$orig.ident,levels=c("Young","Middle","Aged"))
pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/Aging-up-gene-groupbycelltype-colorbyage.pdf")
cds_subset <- monocle_cds[gene[1:8],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",nrow= 4,ncol=2,cell_size =0, plot_trend = TRUE)+
scale_color_manual(breaks = c("Capillary_2","Capillary_3"), values=c("#fed049","#007580")) + theme(legend.position = "right")


pdf("Aging-up-gene-groupbycelltype-colorbyage-col2row4-changecolandrow.pdf")
for(i in 1:13){
  j<-8*i-7
  k<-8*i
Thegene<-Aged_DEG[j:k]
print(Thegene)
cds_subset <- monocle_cds[Thegene,]
p<-plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Arterial",  "C_A" , "Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")
print(p)
}
cds_subset <- monocle_cds[Aged_DEG[100:107],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Arterial",  "C_A" , "Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

dev.off();

pdf("Aging-down-gene-groupbycelltype-colorbyage-col2row4-changecolandrow.pdf")
for(i in 1:3){
  j<-8*i-7
  k<-8*i
Thegene<-Young_DEG[j:k]
print(Thegene)
cds_subset <- monocle_cds[Thegene,]
p<-plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Arterial",  "C_A" , "Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")
print(p)
}
cds_subset <- monocle_cds[Young_DEG[23:30],]
plot_genes_jitter(cds_subset,grouping = "orig.ident",color_by = "celltype",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Arterial",  "C_A" , "Capillary","C_V","Venous"), 
  values=c("#206A5D","#81B264","#FFCC29","#F58634","#BE0000")) + theme(legend.position = "right")

dev.off();

pdf("Aging-down-gene-groupbycelltype-colorbyage-col2row4.pdf")
for(i in 1:3){
  j<-8*i-7
  k<-8*i
Thegene<-Young_DEG[j:k]
print(Thegene)
cds_subset <- monocle_cds[Thegene,]
p<-plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
print(p)
}
cds_subset <- monocle_cds[Young_DEG[23:30],]
plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")

dev.off();


#####Paper gene
Young<-c("Tfrc","Lepr","Insr","Ager","Ldlr","Scarb1","Lrp1a","Lrp8","Lrp10","Slc7a5b","Slc3a2b","Slc2a1","Igf1r","Bsg","Clta","Cltb","Cltc",
	"Ap2a1","Ap2a2","Ap2b1","Ap2m1","Ap2s1","Picalm","Dab2","Hip1","Hip1r","Epn1","Numb","Ldlrap1","Eps15","Rin3","Mfsd2a")
Aged<-c("Cav1","Cav2","Cavin1","Cavin2","Cavin3")
all<-rownames(monocle_cds)
Young<-Young[which(Young%in% all)]

pdf("paper-Aging-gene.pdf")

cds_subset <- monocle_cds[Young[1:8],]
plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
cds_subset <- monocle_cds[Young[9:16],]
plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
cds_subset <- monocle_cds[Young[17:24],]
plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")
cds_subset <- monocle_cds[Young[25:29],]
plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")

cds_subset <- monocle_cds[Aged[1:5],]
plot_genes_jitter(cds_subset,grouping = "celltype",color_by = "orig.ident",
 nrow=4,ncol=2,plot_trend = TRUE,cell_size=0)+
scale_color_manual(breaks = c("Young","Middle","Aged"), 
  values=c("#009F86","#FF9F40","#E64A35")) + theme(legend.position = "right")

dev.off();



level<-c("Young","Middle1","Old")
EC_cells.integrated$orig.ident<-factor(EC_cells.integrated$orig.ident,levels=level)
cols<-c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35")

pdf("/md01/nieyg/project/BBB/YMO_results/fig4_6/Agegenes_heatmap_forsample.pdf")
DoHeatmap(EC_cells.integrated,Aged_DEG,group.by = "orig.ident",group.colors=cols,disp.min=-0.5,disp.max=0.5)+
scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(EC_cells.integrated,Young_DEG,group.by = "orig.ident",group.colors=cols,disp.min=-0.5,disp.max=0.5)+
scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))


dev.off()

c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29")