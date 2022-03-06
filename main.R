# setwd("~/Documents/GitHub/nci60_gsea/")

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

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Shiny
# ----------------------------------------------------------------------------

gene_set_database_options <- c("Hallmark", "C5_gobp", "C5_gomf", "C5_gocc",
                               "C2_cgp", "combined_gs")

path_task_analysis_results <- "data/flux_results_models_publication.txt"
# input variables
active_pathway_options <- vroom(path_task_analysis_results,
                                delim = "\t"
) %>%
  mutate(
    sum =
      dplyr::select(., is_numeric) %>%
      rowSums(na.rm = TRUE)
  ) %>%
  dplyr::filter(sum != 61000) %>%
  dplyr::filter(sum != 0) %>%
  dplyr::select(Description) %>%
  as.list() %>%
  unlist() %>%
  unname()

C2_cgp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
#combined_gs <- get_combined_gs()


ui <- fluidPage(
  titlePanel("Clustering"),
  sidebarLayout(
    sidebarPanel(
        selectInput(
          inputId = "selected_pathway", label = "Select active pathway",
          choices = active_pathway_options, selected = "CMP-N-acetylneuraminate synthesis"),
        selectInput(
          inputId = "selected_genesetdb", label = "Select geneset database",
          choices = gene_set_database_options, selected = "Hallmark"),
        textInput(
          inputId = "selected_subset", label = "enter geneset", placeholder = "ANGIOGENESIS"), 
        actionButton("button", "Submit")
    ),
    mainPanel(
      tabsetPanel(
        plotOutput("pca_plot"),
        plotOutput("heatmap")
      ))
    )
)


server <- function(input, output) {
  
  
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
  
  
  observeEvent(input$button, {
    search_term <- input$selected_subset
    gene_set_1 <- input$selected_genesetdb
    
    pathway_activity <- read.csv("data/flux_results_models_publication.txt",
                                 sep = "\t",
                                 colClasses = c(rep("character", 3), rep("double", 59))) %>%
      dplyr::select(!c("System", "Subsystem")) %>%
      dplyr::filter(Description == input$selected_pathway) %>%
      pivot_longer(!Description, names_to = "cellline", values_to = "phenotype") %>%
      column_to_rownames("cellline")
    
    
    # find better way
    if (gene_set_1 == "Hallmark") {
      gene_setdb <- Hallmark
    } else if (gene_set_1 == "C5_gobp") {
      gene_setdb <- C5_gobp
    } else if (gene_set_1 == "C5_gomf") {
      gene_setdb <- C5_gomf
    } else if (gene_set_1 == "combined_gs") {
      gene_setdb <- combined_gs
    } else if (gene_set_1 == "C2_cgp") {
      gene_setdb <- C2_cgp
    } else {
      gene_setdb <- C5_gocc
    }
    
    #search_term <- "ANGIOGENESIS"
    #gene_setdb <- Hallmark
    geneset_genes <-gene_setdb %>% 
      as.data.frame() %>%
      #dplyr::filter(gs_name == geneset) %>% 
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
   # pca_results <- prcomp(pca_matrix)
    
  output$pca_plot <- renderPlot({
    autoplot(pca_results, data = pathway_activity, colour = "phenotype")
    #returN
    #print(hist(rnorm(100)))
  })
  })
  observeEvent(input$button_heatmap, {
    search_term <- input$selected_subset_heatmap
    gene_set_1 <- input$selected_genesetdb_heatmap
    
    
    
    pathway_activity <- read.csv("data/flux_results_models_publication.txt",
                                 sep = "\t",
                                 colClasses = c(rep("character", 3), rep("double", 59))) %>%
      dplyr::select(!c("System", "Subsystem")) %>%
      dplyr::filter(Description == input$selected_pathway) %>%
      pivot_longer(!Description, names_to = "cellline", values_to = "phenotype") %>%
      column_to_rownames("cellline")
    
    
    # find better way
    if (gene_set_1 == "Hallmark") {
      gene_setdb <- Hallmark
    } else if (gene_set_1 == "C5_gobp") {
      gene_setdb <- C5_gobp
    } else if (gene_set_1 == "C5_gomf") {
      gene_setdb <- C5_gomf
    } else if (gene_set_1 == "combined_gs") {
      gene_setdb <- combined_gs
    } else if (gene_set_1 == "C2_cgp") {
      gene_setdb <- C2_cgp
    } else {
      gene_setdb <- C5_gocc
    }
    
    #search_term <- "ANGIOGENESIS"
    #gene_setdb <- Hallmark
    geneset_genes <-gene_setdb %>% 
      as.data.frame() %>%
      #dplyr::filter(gs_name == geneset) %>% 
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
    # pca_results <- prcomp(pca_matrix)
    
    output$pca_plot <- renderPlot({
      autoplot(pca_results, data = pathway_activity, colour = "phenotype")
      #returN
      #print(hist(rnorm(100)))
    })
  })
}


shinyApp(ui, server)
