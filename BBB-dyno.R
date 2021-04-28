docker run -it --memory=20G --memory-swap=20G --oom-kill-disable -v /home/zanyuan:/zanyuan -v /tmp:/tmp -v /usr/bin/:/dockerbin -v /var/run/docker.sock:/var/run/docker.sock seqyuan/seqyuan-r:v0.0.1 /bin/bash


docker run -it --memory=20G --memory-swap=20G --oom-kill-disable -v /home/zanyuan:/zanyuan -v /tmp:/tmp -v /usr/bin/:/dockerbin -v /var/run/docker.sock:/var/run/docker.sock seqyuan/seqyuan-r:v0.0.1 


ERROR: error during connect: This error may indicate that the docker daemon is not running.: Get http://%2F%2F.%2Fpipe%2Fdocker_engine/v1.24/info: open //./pipe/docker_engine: The system cannot find the file specified.
errors pretty printing info

docker run -it seqyuan/seqyuan-r:v0.0.1 /bin/bash

docker cp YMO_onlyEC_integrated.rds  cdd7edee58be:/ 

library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
sdata <- readRDS(file = "EC_five.rds")
#添加raw counts和normalised expression
#seurat的矩阵需要进行行列转换，以使行为细胞，列为基因
dataset <- wrap_expression(
  counts = t(sdata@assays$RNA@counts),
  expression = t(sdata@assays$RNA@data)
)
dataset <- add_prior_information(
  dataset,
  start_id = "Young_AAACGCTCACAGTCCG-1"
)

#添加数据的cluster信息，这里我们直接用“seurat_clusters”即可
dataset <- add_grouping(
   dataset,
   sdata$celltype
)
model_paga_tree <- infer_trajectory(dataset, "paga_tree")
model_wishbone <- infer_trajectory(dataset, "wishbone")
model_scorpius <- infer_trajectory(dataset, "scorpius")
model_wanderlust <- infer_trajectory(dataset, "wanderlust")
#####slingshot####
model <- model_scorpius
pdf("test_model.pdf")
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
model<-model_wishbone
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
model<-model_wanderlust
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)

dev.off();
