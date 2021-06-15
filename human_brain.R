###########human brain#########
#####GSE97930#########
## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

##Load datasets and create Seurat objects with the raw (non-normalized data).
CerebellarHem<-read.table("GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt",header=TRUE)
FrontalCortex<-read.table("GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt",header=TRUE)
VisualCortex<-read.table("GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt",header=TRUE)
fetalbrain<-read.table("GSE120046_brain_normlized.txt",header=TRUE,row.names=1)
CerebellarHem <-  CreateSeuratObject(counts = CerebellarHem, min.cells=3,
	project="CerebellarHem", assay = "RNA")
FrontalCortex <-  CreateSeuratObject(counts = FrontalCortex, min.cells=3,
	project="FrontalCortex", assay = "RNA")
VisualCortex  <-  CreateSeuratObject(counts = VisualCortex, min.cells=3,
	project="VisualCortex", assay = "RNA")
fetalbrain  <-  CreateSeuratObject(counts = fetalbrain, min.cells=3,
	project="fetalbrain", assay = "RNA")
m<-read.table("GSE120046_metadata.txt",header=T,row.names=1)
saveRDS(fetalbrain,"fetalbrain.rds")
#####添加sample信息####
objList <- list(CerebellarHem,FrontalCortex,VisualCortex)
#######pre-processing#########
#######      QC      #########
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#注意参考基因组里面线粒体相关基因的名称
CerebellarHem[["percent.mt"]] <- PercentageFeatureSet(CerebellarHem, pattern = "^MT")
VisualCortex[["percent.mt"]] <- PercentageFeatureSet(VisualCortex, pattern = "^MT")
FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT")

# Visualize QC metrics as a violin plot
pdf("/md01/nieyg/project/BBB/human_brain/plot/QC/human_adult_brain_qc_plot_nopoint.pdf")
VlnPlot(VisualCortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,)
VlnPlot(FrontalCortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,)
VlnPlot(CerebellarHem, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,)

# Retain cells that express not more than 4000 genes (remove potential homotypic doublets ?) and had mitochondrial content <15%
VisualCortex = subset(VisualCortex,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )
FrontalCortex = subset(FrontalCortex,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )
CerebellarHem = subset(CerebellarHem,nFeature_RNA >350 & nFeature_RNA < 4000 & percent.mt < 15 )
CerebellarHem@meta.data$pos<-"CerebellarHem"
VisualCortex@meta.data$pos<-"VisualCortex"
FrontalCortex@meta.data$pos<-"FrontalCortex"

# Simply merge Seurat objects
merged_obj <- merge(x=VisualCortex,y=c(FrontalCortex,CerebellarHem),add.cell.ids = c("VisualCortex","FrontalCortex","CerebellarHem"),project = "human_adult_brain")
#Idents(merged_obj) <- gsub("_.*", "", colnames(merged_obj))
human_adult_brain.list <- SplitObject(merged_obj,split.by = "pos")

human_adult_brain.anchors <- FindIntegrationAnchors(object.list = human_adult_brain.list,anchor.features = 2000,dims = 1:30)
human_adult_brain.integrated <- IntegrateData(anchorset = human_adult_brain.anchors, dims = 1:30,features.to.integrate = rownames(human_adult_brain.list[[1]]))
head(human_adult_brain.integrated@meta.data)
DefaultAssay(human_adult_brain.integrated) <- "integrated"

# scale and center features in the dataset
human_adult_brain.integrated <- ScaleData(human_adult_brain.integrated, features =rownames(human_adult_brain.integrated),verbose = FALSE)

# Perform linear dimensional reduction
human_adult_brain.integrated <- RunPCA(human_adult_brain.integrated, npcs = 50, verbose = FALSE)
human_adult_brain.integrated <- RunUMAP(human_adult_brain.integrated, reduction = "pca", dims = 1:30)   
pdf("/md01/nieyg/project/BBB/human_brain/plot/human_brain-Umap-sample.pdf.pdf")
p1 <- DimPlot(human_adult_brain.integrated, reduction = "umap", group.by = "ident")
p2 <- DimPlot(human_adult_brain.integrated, reduction = "umap", group.by = "pos", label = TRUE, repel = TRUE) + NoLegend()
p1 
p2
dev.off()
DefaultAssay(human_adult_brain.integrated) <- "RNA"
saveRDS(human_adult_brain.integrated,"human_adult_brain.integrated.rds")
# scale and center features in the dataset
all<-readRDS("human_adult_brain.integrated.rds")
fetalbrain <- FindVariableFeatures(fetalbrain, selection.method = "vst", nfeatures = 2000)

fetalbrain <- ScaleData(fetalbrain, features =rownames(fetalbrain),verbose = FALSE)
fetalbrain <-    RunPCA(fetalbrain, npcs = 50, verbose = FALSE)
fetalbrain <-   RunUMAP(fetalbrain, reduction = "pca", dims = 1:30)

pdf("/md01/nieyg/project/BBB/human_brain/plot/fetalbrain-Umap-sample.pdf.pdf")
p1 <- DimPlot(fetalbrain, reduction = "umap", group.by = "ident")
p2 <- DimPlot(fetalbrain, reduction = "umap", group.by = "type", label = TRUE, repel = TRUE) + NoLegend()
p3 <- DimPlot(fetalbrain, reduction = "umap", group.by = "Pos", label = TRUE, repel = TRUE) + NoLegend()

p1 
p2
p3
dev.off()

Decrease<-c("Tfrc","Mfsd2a","Tjp1","Slco1c1","Slco1a4","Slc29a1","Slco2b1")
Increase<-c("B2m","Tspo","Vwf","Ckb","Il18","Alpl")
features = c('Bsg','Abcg2','Abcb1a',"Slc16a4",'Slc30a1','Slc16a1','Slco1c1','Slc2a1',"Gpd2","Nt5c2",
  	'Ddc','Eogt','Isyna1','Ocln','Cgnl1','Sorbs2','Dnm3','Palm','Tfrc','Igf1r','Tsc22d1','Esyt2',"Palmd")
AgedDEG<-c("Ttr","Acer2","Depp1","Tmem252","Fmo2","Vwf","H2-K1","H2-Q7","Pglyrp1","Plat","Xdh","Aldoa","Ly6a",
	"Fxyd5","Gadd45g","Grrp1","Hbb-bs","Mt1","Apold1","Bnip3","Cp","H2-Q6","Kbtbd11","Litaf","Mt2","Scgb3a1","Selenop","Stra6")
YoungDEG<-c("Macf1","Gm42418","Ivns1abp","AY036118","Ndnf")
Decrease<-toupper(Decrease)
Increase<-toupper(Increase)
features<-toupper(features)
AgedDEG<-toupper(AgedDEG)
YoungDEG<-toupper(YoungDEG)
pdf("/md01/nieyg/project/BBB/human_brain/plot/Human_brain_Decrease_umap.pdf")
FeaturePlot(human_adult_brain.integrated, features = Decrease)
FeaturePlot(fetalbrain, features = Decrease)
dev.off()
pdf("/md01/nieyg/project/BBB/human_brain/plot/Human_brain_Increase_umap.pdf")
FeaturePlot(human_adult_brain.integrated, features = Increase)
FeaturePlot(fetalbrain, features = Increase)
dev.off()
pdf("/md01/nieyg/project/BBB/human_brain/plot/Human_brain_features_umap.pdf")
FeaturePlot(human_adult_brain.integrated, features = features)
FeaturePlot(fetalbrain, features = features)
dev.off()
pdf("/md01/nieyg/project/BBB/human_brain/plot/Human_brain_AgedDEG_umap.pdf")
FeaturePlot(human_adult_brain.integrated, features = AgedDEG)
FeaturePlot(fetalbrain, features = AgedDEG)
dev.off()
pdf("/md01/nieyg/project/BBB/human_brain/plot/Human_brain_YoungDEG_umap.pdf")
FeaturePlot(human_adult_brain.integrated, features = YoungDEG)
FeaturePlot(fetalbrain, features = YoungDEG)
dev.off()
new.cluster.ids <- factor(fetalbrain$type)
head(new.cluster.ids)
fetalbrain@active.ident <- new.cluster.ids
fetalEC <- subset(fetalbrain,idents=c("Endo"))
adultEC <- subset(human_adult_brain.integrated,idents=c("End"))
pdf("/md01/nieyg/project/BBB/human_brain/plot/UmapEC-sample.pdf.pdf")
p1 <- DimPlot(fetalEC, reduction = "umap", group.by = "ident")
p2 <- DimPlot(adultEC, reduction = "umap", group.by = "ident")
pdf("/md01/nieyg/project/BBB/human_brain/plot/EC_YoungDEG_umap.pdf")
FeaturePlot(adultEC, features = YoungDEG)
FeaturePlot(fetalEC, features = YoungDEG)
dev.off()

adultEC<-readRDS("human_adult_brain.integrated.rds")
EC_cells <- subset(adultEC,idents=c("End"))


saveRDS(EC_cells,"human_adult_brain.EC_cells.rds")
EC_cells <- readRDS("human_adult_brain.EC_cells.rds")
#######gene Umap###########
pdf("/md01/nieyg/project/BBB/human_brain/plot/t.pdf")
FeaturePlot(EC_cells, features = "A1BG")


#count<-as.matrix(EC_cells@assays$RNA)
#count<-t(count)
gene<-as.data.frame(EC_cells@assays$RNA["A1BG",])
gene<-t(gene)
gene<-as.data.frame(gene)

p<-ggplot(gene, aes(x = gene[,1]))
p + geom_density(color = “black”, fill = “gray”)+theme_bw() + theme(panel.grid=element_blank())

ggplot(dat,aes(x=Num))+
  geom_density(aes(fill=as.character(dat$Sample),color=as.character(dat$Sample)),alpha = 0.5,size=1,linetype="solid")

EC_cells <- readRDS("/md01/nieyg/project/BBB/human_brain/human_adult_brain.EC_cells.rds")
    tabPanel("Adulet Human Brain EC",
             fluidRow(
                 column(10,
                        wellPanel(
                            textInput("human_gene", label = h3("Gene symbols"), value = "Sox17"),
                            plotOutput("human_featureplot",height="400px",width = "400px"),
                            plotOutput("human_density",height="400px",width = "600px")))         
                 
             )
    ),

    output$human_density <- renderPlot({
        gene<-as.data.frame(EC_cells@assays$RNA[input$human_gene,])
        gene<-t(gene)
        gene<-as.data.frame(gene)       
        p<-ggplot(gene, aes(x = gene[,1]))+labs(x="log2(RPKM+1)",y="Density",title=input$human_gene)
        p + geom_density(color = "black")+theme_bw() + theme(panel.grid=element_blank())
    })
    output$human_featureplot <- renderPlot({
        FeaturePlot(EC_cells, features = input$human_gene)
    })


