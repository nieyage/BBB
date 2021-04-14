library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
sdata <- readRDS(file = "C:/Users/YAGE NIE/Desktop/YMO_onlyEC_integrated.rds")
#添加raw counts和normalised expression
#seurat的矩阵需要进行行列转换，以使行为细胞，列为基因
dataset <- wrap_expression(
  counts = t(sdata@assays$RNA@counts),
  expression = t(sdata@assays$RNA@data)
)
dataset <- add_prior_information(
  dataset,
  start_id = "Young_AAACCCACACAAATGA-1"
)

#添加数据的cluster信息，这里我们直接用“seurat_clusters”即可
dataset <- add_grouping(
   dataset,
   sdata$seurat_clusters
)
model_paga_tree <- infer_trajectory(dataset, "paga_tree")
model_paga <- infer_trajectory(dataset, "paga")
model_slingshot <- infer_trajectory(dataset, "slingshot")
model_tscan <- infer_trajectory(dataset, "tscan")
#####slingshot####
model <- model_slingshot
pdf("test-four-methods.pdf")
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
model <- model_paga_tree
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
model <- model_paga
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
dev.off();