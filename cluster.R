library(mixOmics)
library(janitor)
library(bioDist)
library(ComplexHeatmap)
library(tidyverse)

#-------------------------------------------------------------------------------
# cluster geneexpression add activity state of reactions from tasks analysis
#-------------------------------------------------------------------------------

geneexpression <- vroom("data/nci_60_expression.txt", delim = "\t") %>%
  # remove ME MDA N as no expression data is available
  dplyr::select(-"ME_MDA_N", -starts_with("X"), -starts_with("..")) %>%
  # filter ensembl ids
  dplyr::filter(grepl("ENSG", Identifier)) %>%
  plyr::ddply("Identifier", numcolwise(mean))


search_term <- "EPITHEL"
geneset_genes <- Hallmark %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(gs_name, search_term)) %>%
  dplyr::select(ensembl_gene) %>%
  as.vector() %>%
  as.list() %>%
  unlist()

pathway_activity <- read.csv("data/flux_results_models_publication.txt",
  sep = "\t",
  colClasses = c(rep("character", 3), rep("double", 59))
) %>%
  dplyr::select(!c("System", "Subsystem")) %>%
  column_to_rownames("Description") %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  apply(1, function(x) (x - min(x)) / (max(x) - min(x))) %>%
  as.data.frame()


gen_heatmap_data <- geneexpression %>%
  subset(Identifier %in% geneset_genes) %>%
  dplyr::select(!Identifier)

gen_heatmap_data <- as.matrix(sapply(gen_heatmap_data, as.numeric))

extra_label <- HeatmapAnnotation(
  "CMP" = pathway_activity$`CMP-N-acetylneuraminate synthesis`,
  "Anchor Protein" = pathway_activity$`Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein`,
  "Glucocerebroside" = pathway_activity$`Synthesis of glucocerebroside`,
  "N-Acetylglucosamine synthesis" = pathway_activity$`N-Acetylglucosamine synthesis`,
  "Glucuronate synthesis" = pathway_activity$`Glucuronate synthesis (via inositol)`,
  "Ceramide synthesis" = pathway_activity$`Ceramide synthesis`,
  "Glucosaminyl-acylphosphatidylinositol" = pathway_activity$`Phosphatidyl-inositol to glucosaminyl-acylphosphatidylinositol`
)

Heatmap(gen_heatmap_data,
  name = "Hallmark - EMT",
  top_annotation = extra_label
)
