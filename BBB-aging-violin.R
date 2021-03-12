library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(monocle)
EC_cells.integrated<-readRDS("BBB_onlyEC.rds")
########aging相关的基因在young和old的单细胞中是否差异#######

Decrease<-c("Tfrc","Mfsd2a","Tjp1","Slco1c1","Slco1a4","Slc29a1","Slco2b1")
Increase<-c("B2m","Tspo","Vwf","Ckb","Il18","Alpl")
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VllnPlot(obj, features = feature,split.by = "orig.ident", pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               #+geom_boxplot(width=0.2)+
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
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
#配色方案
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
pdf("/md01/nieyg/project/BBB/plot/Aging_Decrease_genes_Violin.pdf")
StackedVlnPlot(EC_cells.integrated, Decrease[1:3], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Decrease[4:6], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Decrease[7], pt.size=0, cols=my36colors)
dev.off()
pdf("/md01/nieyg/project/BBB/plot/Aging_Increase_genes_Violin.pdf")
StackedVlnPlot(EC_cells.integrated, Increase[1:3], pt.size=0, cols=my36colors)
StackedVlnPlot(EC_cells.integrated, Increase[4:6], pt.size=0, cols=my36colors)
dev.off();





