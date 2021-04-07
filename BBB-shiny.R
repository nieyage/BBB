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

# Load data
YMO<-readRDS("/md01/nieyg/project/BBB/YMO_results/YMO_onlyEC_integrated.rds")
YMO$celltype <- factor(YMO$celltype,levels=c("Arterial","C_A","Capillary_1","Capillary_2","Capillary_3","Capillary_4","C_V_1","C_V_2","Venous","Interferon","Choroid_plexus"))
YMO$orig.ident <- factor(YMO$orig.ident,levels=c("Young","Middle1","Old"))

# Define UI for application 
ui <- navbarPage(
    title = "BBB scRNA-seq",  
    tabPanel("BBB EC and other celltype",
             fluidRow(
                 column(4,
                        wellPanel(
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
    )
)
# Define server logic required to draw a histogram

server<-function(input, output) {
    output$distPlot <- renderPlot({
        DimPlot(YMO,  label = TRUE, pt.size = 0.5) + labs(title = "Brain EC subtypes")
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
                VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="celltype",pt.size = 0,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis()+geom_boxplot(width=0.2)  
            }
            else{
                VlnPlot(YMO,features =input$whole_gene,assay="RNA",group.by="celltype",pt.size = 0.5,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis() 
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
                VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="celltype",pt.size = 0,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis()+geom_boxplot(width=0.2)  
            }
            else{
                VlnPlot(BBB_all,features =input$all_gene,assay="RNA",group.by="celltype",pt.size = 0.5,col=c("Young" = "#00A087", "Middle1"="#FFA040","Old" = "#E64B35"))+RotatedAxis() 
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
}

# Run the application 
shinyApp(ui = ui, server = server)
