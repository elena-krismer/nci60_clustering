library(shiny)
library(tidyverse)
library(plyr)
library(vroom)
library(msigdbr)
library(tidyverse)
library(msigdbr)
library(shiny)
library(plyr)
library(data.table)
library(vroom)
library(qdapTools)
library(corpcor)
library(qvalue)
library(splines)
library(ggfortify)
library(ggplot2)
library(janitor)
library(dplyr)
library(ComplexHeatmap)

source("code/gene_set.R")

active_pathway_options <- read.csv("data/flux_results_models_publication.txt", sep = "\t") %>%
  dplyr::select(Description) %>%
  as.list() %>%
  unlist() %>%unname()

gene_set_database_options <- c("Hallmark", "C5_gobp", "C5_gomf", "C5_gocc",
                               "C2_cgp", "combined_gs")

ui <- shinyUI(fluidPage(

  titlePanel("Clustering"),
  
  # Sidebar with a button
  sidebarLayout(
    sidebarPanel(
      # pca
      tags$form(
        selectInput(
        inputId = "selected_pathway", label = "Select active pathway",
        choices = active_pathway_options, selected = "CMP-N-acetylneuraminate synthesis"),
        selectInput(
          inputId = "selected_db_pca", label = "Select Gene Database",
          choices = gene_set_database_options, selected = "Hallmark"),
        textInput(
        inputId = "selected_subset", label = "enter geneset", placeholder = "ANGIOGENESIS"),
        actionButton("button_pca", "Submit")),
      # cluster
      tags$form(
        selectInput(
          inputId = "selected_db_heatmap", label = "Select Gene Database",
          choices = gene_set_database_options, selected = "Hallmark"),
        textInput(
          inputId = "search_term_heatmap", label = "Enter geneset", placeholder = "ANGIOGENESIS"),
        actionButton("button_heatmap", "Submit")
        
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("PCA",   plotOutput("pca_plot")),
        tabPanel("heatmap",   plotOutput("heatmap"))
      )
   )
)))


server <- shinyServer(function(input, output) {
  
  
  geneexpression <-vroom("data/nci_60_expression.txt", delim = "\t") %>%
    # remove ME MDA N as no expression data is available
    dplyr::select(-"ME_MDA_N", -starts_with("X"), -starts_with("..")) %>%
    # filter ensembl ids
    dplyr::filter(grepl("ENSG", Identifier)) %>%
    plyr::ddply("Identifier", numcolwise(mean))
  
  Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  C5_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  C5_gomf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
  C5_gocc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
  C2_cgp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
  combined_gs <- get_combined_gs()
  
  
  observeEvent(input$button_pca, {
    
    selected_pathway <- input$selected_pathway
    search_term <- input$selected_subset
    selected_db <- req(input$selected_db_pca)
    
    pathway_activity <- read.csv("data/flux_results_models_publication.txt",
                                 sep = "\t",
                                 colClasses = c(rep("character", 3), rep("double", 59))) %>%
      dplyr::select(!c("System", "Subsystem")) %>%
      column_to_rownames("Description") %>%
      mutate_all(function(x) as.numeric(as.character(x))) %>%
      apply(1, function(x)(x-min(x))/(max(x)-min(x))) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Description") %>%
      dplyr::filter(Description == selected_pathway) %>%
      pivot_longer(!Description, names_to = "cellline", values_to = "phenotype") %>%
      column_to_rownames("cellline")
    
    if (selected_db == "Hallmark") {
      gene_setdb <- Hallmark
    } else if (selected_db == "C5_gobp") {
      gene_setdb <- C5_gobp
    } else if (selected_db == "C5_gomf") {
      gene_setdb <- C5_gomf
    } else if (selected_db == "combined_gs") {
      gene_setdb <- combined_gs
    } else if (selected_db == "C2_cgp") {
      gene_setdb <- C2_cgp
    } else {
      gene_setdb <- C5_gocc
    }
    
    geneset_genes <- gene_setdb %>% 
      as.data.frame() %>%
      dplyr::filter(str_detect(gs_name, search_term)) %>%
      dplyr::select(ensembl_gene) %>% 
      as.vector() %>% 
      as.list() %>%
      unlist()
  
  pca_matrix <- geneexpression %>%
      filter(Identifier %in% geneset_genes) %>%
      column_to_rownames("Identifier") %>%
      as.matrix() %>%
      t()
  
  pca_results <- prcomp(pca_matrix)  
  output$pca_plot <- renderPlot({
    autoplot(pca_results, data = pathway_activity, colour = "phenotype")

    })
  })
  
  observeEvent(input$button_heatmap, {
    
    genedb_selected <- req(input$selected_genedb_heatmap)
    search_term <- input$search_term_heatmap
    
    if (genedb_selected  == "Hallmark") {
      gene_setdb <- Hallmark
    } else if (genedb_selected  == "C5_gobp") {
      gene_setdb <- C5_gobp
    } else if (genedb_selected  == "C5_gomf") {
      gene_setdb <- C5_gomf
    } else if (genedb_selected  == "combined_gs") {
      gene_setdb <- combined_gs
    } else if (genedb_selected  == "C2_cgp") {
      gene_setdb <- C2_cgp
    } else {
      gene_setdb <- C5_gocc
    }
    
    geneset_genes <- genedb %>% 
      as.data.frame() %>%
      dplyr::filter(str_detect(gs_name, search_term)) %>%
      dplyr::select(ensembl_gene) %>% 
      as.vector() %>% 
      as.list() %>%
      unlist()
    
    pathway_activity <- read.csv("data/flux_results_models_publication.txt",
                                 sep = "\t",
                                 colClasses = c(rep("character", 3), rep("double", 59))) %>%
      dplyr::select(!c("System", "Subsystem")) %>%
      column_to_rownames("Description") %>%
      mutate_all(function(x) as.numeric(as.character(x))) %>%
      apply(1, function(x)(x-min(x))/(max(x)-min(x))) %>%
      as.data.frame() 
    
    
    gen_heatmap_data <- geneexpression %>%
      subset(Identifier %in% geneset_genes) %>%
      dplyr::select(!Identifier)
    
    gen_heatmap_data <- as.matrix(sapply(gen_heatmap_data, as.numeric))
    
    extra_label <- HeatmapAnnotation(cmp = pathway_activity$`CMP-N-acetylneuraminate synthesis`,
                                     anchor = pathway_activity$`Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein`,
                                     glucoverbrosie = pathway_activity$`Synthesis of glucocerebroside`,
                                     galactosyl = pathway_activity$`Synthesis of galactosyl glucosyl ceramide (link with ganglioside metabolism)`
    )
    
    output$heatmap <- renderPlot({
      Heatmap(gen_heatmap_data, 
              top_annotation = extra_label)
    
    })
  })
  
})

shinyApp(ui = ui, server = server)




    

