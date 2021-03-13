## Load required packages
library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)
library(magrittr)
########将Old  Young cell number downsample 一致##
#################################################
EC_cells.integrated<-readRDS("YMO_onlyEC_integrated.rds")
DefaultAssay(EC_cells.integrated)<-"RNA"
EC_five<-subset(EC_cells.integrated, idents = c( "Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4","C_V_1","C_V_2","Venous"))
###########pesudotime rm                                                                EC_five@meta.data$subtype<-Idents(EC_five)
DefaultAssay(EC_five)<-"RNA"
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(EC_five@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EC_five@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- monocle::newCellDataSet(data,phenoData = pd,
	featureData = fd,lowerDetectionLimit = 0.5，
	expressionFamily = negbinomial.size()) ####使用RNA count
#当数据的表达量类型是FPKM值，分布类型设为tobit()使用expressionFamily=VGAM:::tobit(Lower=0.1)。
#当数据为UMIs, Transcript counts时，数据分布需要设为负二项分布，即negbinomial.size()。
#由于做差异分析通常用到mRNA couts，常常按照官方的方法转成分子数，使用rpc_matrix <- relative2abs(HSMM, method = "num_genes")。

######MONOCLE3
#monocle_cds = preprocess_cds(monocle_cds, num_dim = 100)
#cds <- as.cell_data_set(EC_five)
#cds <- cluster_cells(cds)
#p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
#p2 <- plot_cells(cds, color_cells_by = "partition", 
#	show_trajectory_graph = FALSE)
#wrap_plots(p1, p2)

# 归一化 
 monocle_cds<- estimateSizeFactors(monocle_cds)
 monocle_cds<- estimateDispersions(monocle_cds)
###Filtering low-quality cells
 monocle_cds <- monocle::detectGenes(monocle_cds, min_expr = 3 )
head(featureData(monocle_cds)@data)
expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
                                     num_cells_expressed >= 10))
pdf("monocle-BB.pdf")
print(head(phenoData(monocle_cds)@data))

#Clustering cells without marker genes 

disp_table <- monocle::dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_cds <- monocle::setOrderingFilter(monocle_cds, 
	unsup_clustering_genes$gene_id)
plot_ordering_genes(monocle_cds)

# Plots the percentage of variance explained by the each component based on PCA from the normalized expression data using the same procedure used in reduceDimension function.
# monocle_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
monocle::plot_pc_variance_explained(monocle_cds, return_all = F) # norm_method='log'
monocle_cds <- reduceDimension(monocle_cds)
monocle_cds <- clusterCells(monocle_cds, num_clusters = 6)#6/10
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

#Trajectory step 2: reduce data dimensionality  
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                            method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")


pdf("monocle-fivesubs-V1.pdf")
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size=0.1)
plot_cell_trajectory(monocle_cds, color_by = "subtype",cell_size=0.1)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size=0.1)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size=0.1)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size=0.1) +
    facet_wrap(~seurat_clusters, nrow =4)
dev.off();
pdf("monocle-fivesubs-V2.pdf")
The23genes <- row.names(subset(fData(monocle_cds),
                               gene_short_name %in% c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1')))
##,'Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
  ###   'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")))
plot_genes_branched_pseudotime(monocle_cds[The23genes,],
                               branch_point = 1,ncol=2,
                               color_by = "subtype")

cds_subset <- monocle_cds[The23genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "subtype",cell_size=0.20)

dev.off();





