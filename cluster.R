library(mixOmics)
library(janitor)
library(bioDist)
library(ComplexHeatmap)
library(tidyverse)

geneexpression <-vroom("data/nci_60_expression.txt", delim = "\t") %>%
  # remove ME MDA N as no expression data is available
  dplyr::select(-"ME_MDA_N", -starts_with("X"), -starts_with("..")) %>%
  # filter ensembl ids
  dplyr::filter(grepl("ENSG", Identifier)) %>%
  plyr::ddply("Identifier", numcolwise(mean))


search_term <- "EPITHEL"


geneset_genes <-Hallmark %>% 
  as.data.frame() %>%
  #dplyr::filter(gs_name == geneset) %>% 
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

gen_heatmap_data <-   as.matrix(sapply(gen_heatmap_data, as.numeric))

extra_label <- HeatmapAnnotation(CMP= pathway_activity$`CMP-N-acetylneuraminate synthesis`,
                                "Anchor Protein"= pathway_activity$`Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein`,
                                Glucocerebroside = pathway_activity$`Synthesis of glucocerebroside`,
                                # "Galactosyl glucosyl ceramide" = pathway_activity$`Synthesis of galactosyl glucosyl ceramide (link with ganglioside metabolism)`,
                              "N-Acetylglucosamine synthesis" = pathway_activity$`N-Acetylglucosamine synthesis`,
                            "Glucuronate synthesis" = pathway_activity$`Glucuronate synthesis (via inositol)`,
                          "Ceramide synthesis" = pathway_activity$`Ceramide synthesis`,
                       "Glucosaminyl-acylphosphatidylinositol"= pathway_activity$`Phosphatidyl-inositol to glucosaminyl-acylphosphatidylinositol`
                    )

Heatmap(gen_heatmap_data, name = "Hallmark - Angiogenesis",
        top_annotation = extra_label)



plsda_for_geneset <- function(geneexpression, gene_database, geneset){
 geneset_genes <- gene_database %>% dplyr::filter(gs_name == geneset) %>% 
    dplyr::select(ensembl_gene) %>% as.list() %>% unlist() %>% unname()
  geneexpression <- geneexpression %>% as.data.frame()
  rownames(geneexpression) <- geneexpression$Identifier
  geneexpression$Identifier <- NULL
  
  geneexpression <- geneexpression %>% 
    t() %>% 
    as.data.frame() %>%
    # first row has gene ids - move them into header
    row_to_names(row_number = 1) %>%
    dplyr::select(intersect(geneset_genes, colnames(geneexpression)))
  
  cellline <- rownames(geneexpression) %>% as.vector()
  rownames(geneexpression) <- NULL
  pls_result <- mixOmics::splsda(X = geneexpression, Y = cellline)
  return(pls_result)
}


plotIndiv(pls_result, 
         centroid =  TRUE,
legend = TRUE)

d.sample.euc <- euc(as.matrix(gen))
heatmap(as.matrix(d.sample.euc),sym=TRUE,
        main='Between-sample distances (Euclidean)',
        xlab='cell type')

d.sample.euc <- euc(t(gen))
heatmap(as.matrix(d.sample.euc),sym=TRUE,col=hmcol,
        main='Between-sample distances (Euclidean)',
        xlab='cell type')




Heatmap(gen, clustering_distance_rows = )
pheatmap(gen_heatmap_data, 
         annotation_col = phenotype)

# PCA

pathway_activity <- read.csv("data/flux_table.txt",
                             sep = "\t",
                             colClasses = c(rep("character", 3), rep("double", 59))) %>%
  dplyr::select(!c("System", "Subsystem")) %>%
  dplyr::filter(Description == "CMP-N-acetylneuraminate synthesis") %>%
  pivot_longer(!Description, names_to = "cellline", values_to = "phenotype") %>%
  column_to_rownames("cellline")

pca_matrix <- geneexpression %>%
  filter(Identifier %in% geneset_genes) %>%
  column_to_rownames("Identifier") %>%
  as.matrix() %>%
  t()
 
pca_results <- prcomp(pca_matrix)
autoplot(pca_results, data = pathway_activity, colour = "phenotype")


