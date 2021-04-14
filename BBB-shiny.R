#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

options(encoding = "UTF-8")

# Load required packages
library(shiny)
library(Seurat)
library(ggplot2)
library(magrittr)
# Load data
BBB_all<-readRDS("/md01/nieyg/project/BBB/BBB_YMO.rds")
BBB_all$celltype <- factor(BBB_all$celltype,levels=c("C_V_1","C_V_2","Capillary_1","Capillary_2","Capillary_3","C_A_1","C_A_2","Arterial_1","Arterial_2","Venous","Interferon",
                                                     "Erythroid","Choroid_plexus","Pericyte","Microglia_1","Microglia_2","Lymphatics","SMC"))
BBB_all$orig.ident <- factor(BBB_all$orig.ident,levels=c("Young","Middle1","Old"))

# Load data
YMO<-readRDS("/md01/nieyg/project/BBB/YMO_results/YMO_onlyEC_integrated.rds")
YMO$celltype <- factor(YMO$celltype,levels=c("Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4","C_V_1","C_V_2","Venous","Interferon","Choroid_plexus"))
YMO$orig.ident <- factor(YMO$orig.ident,levels=c("Young","Middle1","Old"))

####load bulkdata####
df2<-read.table("/md01/nieyg/BBB_RNAseq_barplot_RPKM_df2.txt",header=T)
df3<-read.table("/md01/nieyg/BBB_Proteomics_barplot_RPKM_df3.txt",header=T)
#####load humandata####
EC_cells <- readRDS("/md01/nieyg/project/BBB/human_brain/human_adult_brain.EC_cells.rds")

# Define UI for application 
ui <- navbarPage(
    title = "BBB scRNA-seq",  
    tabPanel("BBB RNAseq and Proteomics",
             fluidRow(
                 column(5,
                        wellPanel(
                            textInput("RNA_gene", label = h3("RNA seq Gene symbols"), value = "Tfrc"),
                            plotOutput("RNA_barplot",height="400px",width = "600px"))) ,
                 column(5,
                        wellPanel(
                            textInput("PRO_gene", label = h3("Proteomics Gene symbols"), value = "Tfrc"),
                            plotOutput("PRO_barplot",height="400px",width = "600px"))) 

                 
                  )
    ),
    tabPanel("Mouse Brain Celltypes",
             fluidRow(
                 column(4,
                        wellPanel(
                            h3(titlePanel("Umap for all celltypes")),
                            h3(titlePanel(".")),
                            plotOutput("all_distPlot",height="400px",width = "450px"))),
                 column(4,
                        wellPanel(
                            selectizeInput(inputId ="all_color_by",label = h3("Color by"),choices =c("Young","Aged","Middle"),selected ="Young"),
                            plotOutput("all_distPlot_age",height="400px",width = "400px"))),
                 column(4,
                        wellPanel(
                            textInput("all_gene", label = h3("Gene symbols"), value = "Sox17"),
                            plotOutput("all_FeaturePlot",height="400px",width = "400px"))),
                 column(12,
                        wellPanel(
                            #plotOutput("VlnPlot",height="450px",width = "850px"),
                            selectizeInput(inputId ="all_group_type",label = h3("Group by"),choices =c("sample_type","Cell type"),selected ="sample_type"),
                            selectizeInput(inputId ="all_add_point",label = h3("add point?"),choices =c("No","Yes"),selected ="Yes"),
                            br(),
                            plotOutput("all_violinplot",height="400px",width = "1200px"),
                            br(),
                            plotOutput("all_violinplot_split",height="400px",width = "1200px")))
                 
                 
             )
    ),
    tabPanel("onlyEC subtypes",
             fluidRow(
                 column(4,
                        wellPanel(
                            h3(titlePanel("Umap for EC subtypes")),
                            h3(titlePanel(".")),
                            plotOutput("distPlot",height="400px",width = "450px"))),
                 column(4,
                        wellPanel(
                            selectizeInput(inputId ="color_by",label = h3("Color by"),choices =c("Young","Aged","Middle"),selected ="Young"),
                            plotOutput("distPlot_age",height="400px",width = "400px"))),
                 column(4,
                        wellPanel(
                            textInput("whole_gene", label = h3("Gene symbols"), value = "Sox17"),
                            plotOutput("FeaturePlot",height="400px",width = "400px"))),
                 column(12,
                        wellPanel(
                            #plotOutput("VlnPlot",height="450px",width = "850px"),
                            selectizeInput(inputId ="group_type",label = h3("Group by"),choices =c("sample_type","Cell type"),selected ="sample_type"),
                            selectizeInput(inputId ="add_point",label = h3("add point?"),choices =c("No","Yes"),selected ="Yes"),
                            br(),
                            plotOutput("violinplot",height="400px",width = "1200px"),
                            br(),
                            plotOutput("violinplot_split",height="400px",width = "1200px")))
                 
                 
             )
    ),
    tabPanel("Adult Human Brain EC",
             fluidRow(
                 column(4,
                        wellPanel(
                            textInput("human_gene", label = h3("Please enter the Gene symbol in Uppercase:"), value = "A2M"),
                            plotOutput("human_featureplot",height="400px",width = "400px"))),      

    column(6,
           wellPanel(
               h3(titlePanel("Density distribution of gene expression")),
               h3(titlePanel("in human brain's Endothelial cells")),
               plotOutput("human_density",height="400px",width = "600px"))
           )         
    
)
)
)
# Define server logic required to draw a histogram

server<-function(input, output) {
    output$distPlot <- renderPlot({
        DimPlot(YMO,  label = TRUE, pt.size = 0.5,cols=c("Arterial"="#206A5D","C_A"="#81B264","Capillary_1"="#FFCC29","Capillary_2"="#FFCC29","Capillary_3"="#FFCC29","Capillary_4"="#FFCC29",
                                                         "C_V_1"="#F58634","C_V_2"="#F58634","Venous"="#BE0000","Choroid_plexus"="#31326F","Interferon"="#93329E")) + labs(title = "Brain EC subtypes")
    })
    output$distPlot_age <- renderPlot({
        Umap = YMO@reductions$umap@cell.embeddings %>%
            as.data.frame() %>% cbind(tx =YMO@meta.data$orig.ident)
        if (input$color_by=="Young"){
            p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
                geom_point(size = 0.2, alpha = 1) + 
                scale_color_manual(values=c("Young" = "#00A087", "Middle1"="grey","Old" = "grey"))   
            p+labs(title="Young cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
        }else if (input$color_by=="Aged"){
            p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
                geom_point(size = 0.2, alpha = 1) + 
                scale_color_manual(values=c("Young" ="grey" , "Middle1"="grey","Old" ="#E64B35" ))   
            p+labs(title="Aged cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
        }else if (input$color_by=="Middle"){
            p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
                geom_point(size = 0.2, alpha = 1) + 
                scale_color_manual(values=c("Young" = "grey", "Middle1"="#FFA040","Old" = "grey"))   
            p+labs(title="Middle cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
        }   
    })
    
    output$VlnPlot<-renderPlot({
        VlnPlot(YMO, features = input$whole_gene, pt.size = 0.2, ncol = 1)
    })
    
    output$FeaturePlot<-renderPlot({
        FeaturePlot(YMO, features = input$whole_gene, pt.size = 0.2, ncol = 1)
    })     
    output$violinplot <- renderPlot({
        if (input$group_type=="sample_type"){
            if(input$add_point=="No"){
                VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="orig.ident",pt.size = 0,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis()+geom_boxplot(width=0.2)  
            }
            else{
                VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="orig.ident",pt.size = 0.5,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis() 
            }
        }else if (input$group_type=="Cell type"){
            if(input$add_point=="No"){
                VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="celltype",pt.size = 0)+RotatedAxis()+geom_boxplot(width=0.2)  
            }
            else{
                VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="celltype",pt.size = 0.5)+RotatedAxis() 
            }        }
    })
    output$violinplot_split <- renderPlot({
        if(input$add_point=="No"){
            VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="celltype",split.by="orig.ident",pt.size = 0,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis()+geom_boxplot(width=0.2,outlier.size=0,position=position_dodge(width=0.9))
        }
        else{
            VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="celltype",pt.size = 0.5,split.by="orig.ident",col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis() 
        }
    })
    #############BBB all celltype 
    
    output$all_distPlot <- renderPlot({
        DimPlot(BBB_all,  label = TRUE, pt.size = 0.5) + labs(title = "Brain all celltypes")
    })
    output$all_distPlot_age <- renderPlot({
        Umap = BBB_all@reductions$umap@cell.embeddings %>%
            as.data.frame() %>% cbind(tx =BBB_all@meta.data$orig.ident)
        if (input$all_color_by=="Young"){
            p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
                geom_point(size = 0.2, alpha = 1) + 
                scale_color_manual(values=c("Young" = "#00A087", "Middle1"="grey","Old" = "grey"))   
            p+labs(title="Young cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
        }else if (input$all_color_by=="Aged"){
            p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
                geom_point(size = 0.2, alpha = 1) + 
                scale_color_manual(values=c("Young" ="grey" , "Middle1"="grey","Old" ="#E64B35" ))   
            p+labs(title="Aged cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
        }else if (input$all_color_by=="Middle"){
            p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
                geom_point(size = 0.2, alpha = 1) + 
                scale_color_manual(values=c("Young" = "grey", "Middle1"="#FFA040","Old" = "grey"))   
            p+labs(title="Middle cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
        }   
    })
    
    output$all_VlnPlot<-renderPlot({
        VlnPlot(BBB_all, features = input$all_gene, pt.size = 0.2, ncol = 1)
    })
    
    output$all_FeaturePlot<-renderPlot({
        FeaturePlot(BBB_all, features = input$all_gene, pt.size = 0.2, ncol = 1)
    })     
    output$all_violinplot <- renderPlot({
        if (input$all_group_type=="sample_type"){
            if(input$all_add_point=="No"){
                VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="orig.ident",pt.size = 0,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis()+geom_boxplot(width=0.2)  
            }
            else{
                VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="orig.ident",pt.size = 0.5,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis() 
            }
        }else if (input$all_group_type=="Cell type"){
            if(input$all_add_point=="No"){
                VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="celltype",pt.size = 0)+RotatedAxis()+geom_boxplot(width=0.2)  
            }
            else{
                VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="celltype",pt.size = 0.5)+RotatedAxis() 
            }        }
    })
    output$all_violinplot_split <- renderPlot({
        if(input$all_add_point=="No"){
            VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="celltype",split.by="orig.ident",pt.size = 0,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis()+geom_boxplot(width=0.2,outlier.size=0,position=position_dodge(width=0.9))
        }
        else{
            VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="celltype",pt.size = 0.5,split.by="orig.ident",col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis() 
        }
    })
    output$RNA_barplot <- renderPlot({
        p1 <- ggplot(df2[which(df2$gene==input$RNA_gene),], aes(x = tissue, y = count,fill=tissue)) +
            geom_bar(stat = "identity", color = "black", position = position_dodge()) +
            geom_errorbar(aes(ymin = count - sd, ymax = count + sd), 
                          width = 0.2, position = position_dodge(0.9))+
            labs(x="",
                 y="RNAseq-RPKM",
                 title=input$RNA_gene)  +
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5), 
                  legend.position="right", 
                  legend.title = element_blank())
        
        p1+ theme_bw() + theme(panel.grid=element_blank())
    })
    output$PRO_barplot <- renderPlot({
        p1 <- ggplot(df3[which(df3$gene==input$PRO_gene),], aes(x = tissue, y = count,fill=tissue)) +
            geom_bar(stat = "identity", color = "black", position = position_dodge()) +
            geom_errorbar(aes(ymin = count - sd, ymax = count + sd), 
                          width = 0.2, position = position_dodge(0.9))+
            labs(x="",
                 y="Proteomics signal",
                 title=input$PRO_gene)  +
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5), 
                  legend.position="right", 
                  legend.title = element_blank())
        
        p1+ theme_bw() + theme(panel.grid=element_blank())
    })
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
}

# Run the application 
shinyApp(ui = ui, server = server)
#rsconnect::setAccountInfo(name='bbb-v1', token='79DD2E0EB821DE578BC349068377A766', secret='0z3BCaAMeGdgDmqYXQXVKVFhnCgNvQ5EflMZuHv7')

