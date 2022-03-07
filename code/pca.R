# PCA
library(tidyverse)
library(janitor)
library(ggfortify)
library(ComplexHeatmap)
library(plyr)

`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# get mean flux CMP-N-Acetylneuroaminate reaction from previous sampling
# analysis combine to one dataframe for clustering and pca
#-------------------------------------------------------------------------------

cmp_reactions <- c(
  "CACNAO", "CACNAO2", "ACNML", "ACNAMtr", "CMPACNAtg", "NS26T2g",
  "NS26Tg", "S23T2g", "S23T3g", "S23T4g", "S23Tg", "S26Tg", "SIAT4Bg",
  "SIAT9g", "ST3GAL21g", "ST3GAL22g", "ST3GAL23g", "ST3GAL31g",
  "ST3GAL61g", "ST3GAL62g",
  "ST6GALNAC21", "ST6GALNAC22", "ST6GALNAC23", "ST6GALNAC24",
  "ST6GALNAC25", "ST6GALNAC26", "ST6GALNAC27", "ST6GALNAC28",
  "ST6GALNAC31", "ST6GALNAC61", "ST6GALNAC62", "ST8SIA11",
  "ST8SIA12", "ST8SIA51g", "ST8SIA53g", "ST8SIA54g", "ST8SIA55g",
  "ST8SIA56g", "EX_cmpacna_e",
  "CMPACNAtn", "CMPSAS", "CMPACNAtg"
)

models <- c(
  "786-0", "A498", "A549", "ACHN", "ADR-RES", "BT-549", "CAKI-1",
  "CCRF-CEM", "COLO-205", "DU-145", "EKVX", "HCC2998", "HCT-15", "HCT-116",
  "HL-60", "HOP-62", "HOP-92", "HS-578-T", "HT29", "IGROV-1", "K-562", "KM12",
  "LOXIMVI", "M14", "MALME-3M", "MCF7", "MDA-MB-231", "MDA-MB-435", "MOLT-4",
  "NCI-H23", "NCI-H226", "OVCAR-3", "NCI-H322M", "NCI-H460", "NCI-H522",
  "OVCAR-4", "OVCAR-5", "OVCAR-8", "PC-3", "RPMI-8226", "RXF393", "SF268",
  "SF295", "SF539", "SK-MEL-2", "SK-MEL-5", "SK-MEL-28", "SK-OV-3",
  "SN12C", "SNB-19", "SNB-75", "SR", "SW620", "T47D", "TK10", "U251",
  "UACC-62", "UACC-257", "UO-31"
)

cmp_sampling_df <- data.frame(matrix(
  ncol = length(cmp_reactions), nrow = 0,
  dimnames = list(NULL, cmp_reactions)
))


for (model in models) {
  path_to_sampling <- paste0("~/Documents/GitHub/nci_60/data/sampling/", model, "_sampling.csv")
  sampling <- read.csv(path_to_sampling)
  model_sampling <- sampling %>%
    select_if(names(.) %in% cmp_reactions) %>%
    colMeans() %>%
    as.data.frame() %>%
    t()
  rownames(model_sampling) <- model
  cmp_sampling_df <- rbind(cmp_sampling_df, model_sampling, fill = TRUE)
}

cmp_sampling_df <- cmp_sampling_df %>%
  filter(!grepl("fill", rownames(cmp_sampling_df)))


# clustering
cmp_matrix <- cmp_sampling_df %>%
  as.matrix() %>%
  t() %>%
  scale()
cmp_heatmap <- Heatmap(cmp_matrix,
  name = "CMP-N-Acetylneuraminate", cluster_rows = TRUE,
  cluster_columns = TRUE, cluster_column_slices = TRUE
)

#-------------------------------------------------------------------------------
# PCA Geneexpression of Angiogenesis genes
# Label - CMP-N-Acetylneuroaminate reactions grouped in 4
#-------------------------------------------------------------------------------

geneexpression <- vroom("data/nci_60_expression.txt", delim = "\t") %>%
  # remove ME MDA N as no expression data is available
  dplyr::select(-"ME_MDA_N", -starts_with("X"), -starts_with("..")) %>%
  # filter ensembl ids
  dplyr::filter(grepl("ENSG", Identifier)) %>%
  plyr::ddply("Identifier", numcolwise(mean))

search_term <- "C5_GOBP"

geneset_genes <- combined_gs %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(gs_name, search_term)) %>%
  dplyr::select(ensembl_gene) %>%
  as.vector() %>%
  as.list() %>%
  unlist()

# perform PCA
gen_pca_data <- geneexpression %>%
  subset(Identifier %in% geneset_genes) %>%
  remove_rownames() %>%
  column_to_rownames("Identifier") %>%
  t() %>%
  as_tibble() %>%
  prcomp(center = TRUE, scale. = TRUE)

# group celllines according to CMP-N-Acetylneuraminate reaction flux
# clustering
cmp_matrix <- cmp_sampling_df %>%
  as.matrix() %>%
  t() %>%
  scale()

h <- heatmap(cmp_matrix %>% t(), keep.dendro = TRUE)
row.clusters <- as.hclust(h$Rowv)
cmp_cluster_4 <- cutree(row.clusters, k = 4) %>%
  as.data.frame() %>%
  mutate("CMP_cluster" = cmp_cluster_4$.)

autoplot(gen_pca_data,
  data = cmp_cluster_4, colour = factor(cmp_cluster_4$CMP_cluster),
  fill = factor(cmp_cluster_4$CMP_cluster),
  frame = TRUE,
  frame.type = "norm"
) +
  # autoplot does not take the color for what reason ever
  scale_color_manual(values = c("#8691C7", "#CF98C8", "#7FC9BA", "#F3828C")) +
  scale_fill_manual(values = c("#8691C7", "#CF98C8", "#7FC9BA", "#F3828C")) +
  labs(color = "Pathway Cluster", fill = "Pathway Cluster", title = "PCA - C5_GOBP Angiogenesis")

#
# EPITHELIAL MESENCHYMAL TRANSITION

search_term <- "EPITH"

geneset_genes <- Hallmark %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(gs_name, search_term)) %>%
  dplyr::select(ensembl_gene) %>%
  as.vector() %>%
  as.list() %>%
  unlist()

# perform PCA
gen_pca_data <- geneexpression %>%
  subset(Identifier %in% geneset_genes) %>%
  remove_rownames() %>%
  column_to_rownames("Identifier") %>%
  t() %>%
  as_tibble() %>%
  prcomp(center = TRUE, scale. = TRUE)

# group celllines according to CMP-N-Acetylneuraminate reaction flux
# clustering
cmp_matrix <- cmp_sampling_df %>%
  as.matrix() %>%
  t() %>%
  scale()

h <- heatmap(cmp_matrix %>% t(), keep.dendro = TRUE)
row.clusters <- as.hclust(h$Rowv)
cmp_cluster_4 <- cutree(row.clusters, k = 4) %>%
  as.data.frame() %>%
  mutate("CMP_cluster" = cmp_cluster_4$.)

autoplot(gen_pca_data,
  data = cmp_cluster_4, colour = factor(cmp_cluster_4$CMP_cluster),
  fill = factor(cmp_cluster_4$CMP_cluster),
  frame = TRUE,
  frame.type = "norm"
) +
  # autoplot does not take the color for what reason ever
  scale_color_manual(values = c("#8691C7", "#CF98C8", "#7FC9BA", "#F3828C")) +
  scale_fill_manual(values = c("#8691C7", "#CF98C8", "#7FC9BA", "#F3828C")) +
  labs(color = "Pathway Cluster", fill = "Pathway Cluster", title = "PCA - Epithelial Mesenchymal Transition") +
  theme_minimal()

#-------------------------------------------------------------------------------
# PCA CMP-N-Acetylneuroaminate flux from sampling analysis
# Label - Hallmark geneexpression clustering
#-------------------------------------------------------------------------------


pca_results <- cmp_sampling_df %>%
  prcomp(center = TRUE, scale. = TRUE)

# without major outliers
# remove CAKI-0 and NCI-H522
pca_cmp_no_outlier <- cmp_sampling_df %>%
  filter(rownames(cmp_sampling_df) %notin% c("CAKI-1", "NCI-H522")) %>%
  prcomp(center = TRUE, scale. = TRUE)

search_term <- "EPITH"

geneset_genes <- Hallmark %>%
  as.data.frame() %>%
  dplyr::filter(str_detect(gs_name, search_term)) %>%
  dplyr::select(ensembl_gene) %>%
  as.vector() %>%
  as.list() %>%
  unlist()

gen_heatmap_data <- geneexpression %>%
  subset(Identifier %in% geneset_genes) %>%
  dplyr::select(!Identifier) %>%
  as.matrix()

gen_heatmap_data <- as.matrix(sapply(gen_heatmap_data, as.numeric))

h <- heatmap(gen_heatmap_data %>% t(), keep.dendro = TRUE)
row.clusters <- as.hclust(h$Rowv)
angiogenesis_cluster_4 <- cutree(row.clusters, k = 4) %>%
  as.data.frame() %>%
  filter(rownames(angiogenesis_cluster_4) %notin% c("CAKI-1", "NCI-H522"))
colnames(angiogenesis_cluster_4) <- "Angiogenesis"


autoplot(pca_cmp_no_outlier,
  data = angiogenesis_cluster_4, colour = factor(angiogenesis_cluster_4$Angiogenesis),
  fill = factor(angiogenesis_cluster_4$Angiogenesis),
  frame = TRUE,
  frame.type = "norm"
) +
  # autoplot does not take the color for what reason ever
  scale_color_manual(values = c("#8691C7", "#CF98C8", "#7FC9BA", "#F3828C")) +
  scale_fill_manual(values = c("#8691C7", "#CF98C8", "#7FC9BA", "#F3828C")) +
  labs(
    color = "EMT Expression Cluster \nHallmark",
    fill = "EMT Expression Cluster \nHallmark",
    title = "PCA - CMP-N-Acetylneuroaminate reactions sampling"
  ) +
  theme_minimal()
