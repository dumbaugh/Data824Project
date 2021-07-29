#### Shiny Project -- David Umbaugh

# Load packages ----
library(shiny)
library(Seurat)
library(stringi)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(tidyverse)


## read in single cell RDS object

seuratObject <- readRDS("Data/AllCombined_1400kdownsampled_V4_3.RDS")

#change active identity to celltype2
Idents(seuratObject) <- 'CellType2'

## extract cell barcodes based on cell identity 
neutrophils <- WhichCells(seuratObject,idents = 'Neutrophils')
immuneCells <- WhichCells(seuratObject,idents='Immune cells')
HSC <- WhichCells(seuratObject,idents='Hepatic stellate cells')
EndoCells <- WhichCells(seuratObject,idents="Endothelial cells")
PCheps <- WhichCells(seuratObject,idents="Pericentral hepatocytes")
PPheps <- WhichCells(seuratObject,idents="Periportal hepatocytes")



# user interface

ui <- fluidPage(
  titlePanel(strong("Exploring Acetaminophen Hepatotoxicity with scRNAseq")),
  navbarPage( 
    NULL,  
    tabPanel(strong("Overview: Data Visualization and Acquisition"),
             sidebarLayout(
      sidebarPanel(h1("APAP Single Cell Project"),
                   h3("About the Data"),
                    p("Acetaminophen overdose is the most common cause of liver failure in the United States, resulting in over 10,000 hospitalizations annually.
To understand how different cell types within the liver respond to acetaminophen hepatotoxicity, single cell RNA sequencing experiments have been
performed. The data in the application comes from three different single cell RNA sequencing experiments:"),
                   h5((a("Spatial reconstruction of liver after APAP using sc-RNAseq",
                         href="https://pubmed.ncbi.nlm.nih.gov/33983442/"))),
                    h5((a("Functional compensation precedes recovery following acute liver injury",
                          href="https://pubmed.ncbi.nlm.nih.gov/33214549/"))),
                    h5((a("Acute liver failure is regulated by MYC- and microbiome dependent programs",
                          href="https://pubmed.ncbi.nlm.nih.gov/33106666/"))),
                     p("Each experiment explored different aspects of acetaminophen hepatotoxicty, therefore, these datasets are
                     complementary and when integrated allow for a more complete picture of the transcriptomic landscape of the liver
                     following drug-induced toxicity. The APAP toxicity timecourse can be divided into an injury phase and recovery phase.
                     The injury phase spans the first 24h after APAP exposure, followed by a recovery phase which lasts 48h to 72h.By combining
                     all publically available scRNAseq datasets that used APAP, the entire timeline can be observed. Additionally, one of the papers
                     used a higher dose allowing for an analysis of the dose-dependent response of cell types.")),
      mainPanel(titlePanel(strong("Transcriptomic landscape of the liver following Acetaminophen Overdose")),
                img(src="dimplott_default.png",height=800,width=1000)
       
      
    )
    )
    ),
    tabPanel(strong("Visualization by UMAP and VlnPlot"),br(),
  
  fluidRow(
    
    column(3,
           helpText("Create UMAP and VlnPlot with APAP scRNAseq data."),

                              textInput("var",label = "Type a gene here",
                                        value="Cyp2e1",placeholder="Cyp2e1"),
           selectInput("select",("Choose an identity class"), 
                        choices = list("BiologicalReplicate","CellType2"),selected="BiologicalReplicate"),
           selectInput("split",("Choose an identity class to split by"), 
                choices = list("None","BiologicalReplicate","CellType2"),selected="None"),
           checkboxGroupInput("cell",("Optional: Choose a cell subpopulation(s) to examine, only effects UMAP"), 
                              choices = list("Pericentral hepatocytes","Periportal hepatocytes","Endothelial cells","Immune cells",
                                             "Hepatic stellate cells"),selected=NULL),
           ),
                                                    
    
    column(9,
           plotOutput("umap"),
    
    
  )),
  
  fluidRow(
    
    column(3,
           helpText("Violin plot to visualize distributions of gene expression within groups.")),
    
    column(9,
           plotOutput("vln")),
    
    
  )
  
  ),
  tabPanel(strong("Dotplot comparing gene expression across groups"),br(),
           
           fluidRow(
    
    column(3,
           helpText("Create a dotplot to compare scaled gene expression levels"),
           textInput("listvar",label = "Type a list of genes separated by a comma (e.g. Cyp2e1,Cyp2f2)",
                     value="Cyp2e1",placeholder="Cyp2e1"),
           selectInput("select2",("Choose an identity class"), 
                       choices = list("BiologicalReplicate","CellType2"),selected="BiologicalReplicate"),
           ),
    
    column(9,
           plotOutput("dot"),
           
  )
),
),
  tabPanel(strong("Feature scatter comparing two different genes"),br(),        
          fluidRow(
          column(3,
                 helpText("Create a correlation plot to compare the relationship between two different genes."),
                 textInput("feature1",label = "Type your first gene (e.g. Cyp2e1)",
                           value="Cyp2e1",placeholder="Cyp2e1"),
                 textInput("feature2",label = "Type your second gene (e.g. Cyp2f2)",
                           value="Cyp2f2",placeholder="Cyp2f2"),
                 selectInput("select3",("Choose an identity class"), 
                             choices = list("BiologicalReplicate","CellType2"),selected="BiologicalReplicate"),
                 checkboxGroupInput("cell2",("Optional: Choose a cell subpopulation(s) to examine"), 
                                  choices = list("Pericentral hepatocytes","Periportal hepatocytes","Endothelial cells","Immune cells",
                                                 "Hepatic stellate cells"),selected=NULL)
          ),
          
          column(9,
                 plotOutput("featurescatter"),
        )
),
)
)
)




# Define server logic ----
server <- function(input, output) {
 #tab1
  output$umap <- renderPlot({
    i <- 1
    cellVector <- c()
    while(i <= length(input$cell)){
      testing <- switch(input$cell[i],"Pericentral hepatocytes"= PCheps,
                        "Periportal hepatocytes"=PPheps, 
                        "Endothelial cells"=EndoCells,
                        "Immune cells" = immuneCells,
                        "Hepatic stellate cells"= HSC)
      cellVector <- c(cellVector,testing)
      i=i+1
    }
    Idents(seuratObject) <- input$select
    x <- str_to_title(input$var)
    if(length(cellVector)==0){
    FeaturePlot(seuratObject,features=x,label=TRUE,pt.size = 1,repel = TRUE,label.size = 6)+FontSize(x.title=20,y.title=20)}else{
      FeaturePlot(seuratObject,features=x,label=TRUE,pt.size = 1,repel = TRUE,cells=cellVector,label.size=6)+FontSize(x.title=20,y.title=20)
    }
  })
  
  
  output$vln <- renderPlot({
    Idents(seuratObject) <- input$select

    x <- str_to_title(input$var)
    if(input$split =='None'){
      VlnPlot(seuratObject,features=x)+theme(legend.position='none')+FontSize(x.title=20,y.title=20)} 
    else{
      VlnPlot(seuratObject,features=x,split.by=input$split)+FontSize(x.title=20,y.title=20)
      }
    
  })
  
  ##tab 2
  output$dot <- renderPlot({
   
    
    Idents(seuratObject) <- input$select2
    y <- unlist(strsplit(input$listvar,split=",")) %>% 
      stri_trim() %>% 
      str_to_title()
      DotPlot(seuratObject,features=c(y))+coord_flip()})
  
  ##tab 3
  output$featurescatter <- renderPlot({
    j <- 1
    cellVector2 <- c()
    while(j <= length(input$cell2)){
      track <- switch(input$cell2[j],
                      "Pericentral hepatocytes"= PCheps,
                      "Periportal hepatocytes"=PPheps,
                      "Endothelial cells"=EndoCells,
                        "Immune cells" = immuneCells,
                        "Hepatic stellate cells"= HSC)
      cellVector2 <- c(cellVector2,track)
      j=j+1
    }
    Idents(seuratObject) <- input$select3
    q <- str_to_title(input$feature1)
    z <- str_to_title(input$feature2)
    if(length(cellVector2)==0){
      FeatureScatter(seuratObject,feature1=q,feature2=z)+FontSize(x.title=20,y.title=20)}
    else{
        FeatureScatter(seuratObject,feature1=q,feature2=z,cells=cellVector2)+FontSize(x.title=20,y.title=20)
      }
    
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)