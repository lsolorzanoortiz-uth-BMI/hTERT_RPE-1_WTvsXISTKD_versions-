# Set CRAN mirror to avoid prompts during installation #####
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Install WGCNA and dependencies from CRAN #####
install.packages("WGCNA", dependencies = TRUE)

install.packages(c(
  "BiocManager",      # Manages the installation of Bioconductor 
  "ggplot2",          # For vizualization
  "gplots",           # For heatmap.2
  "reshape2",         # For data manipulation
  "igraph",           # For network analysis and visualization
  "pheatmap",         # For heatmap visualization
  "RColorBrewer",     # For color palettes
  "corrplot",         # For correlation plots
  "ggrepel" ,         # For non-overlapping labels in plots
  "doParallel",       # For parallel processing
  "dplyr",            # Data manipulation and transformation
  'GOplot',           # Enrichment Analysis
  "gggenes",          # Plotting 
  "circlize",
  "viridisLite",      # Color
  "ggridges",         # For visualizing the GSEA 
  'VennDiagram',      # Plotting sig results
  "purrr",
  "tidyestimate",     # For immune scoring 
  "remotes",          # Dependency for immune scoring with omnideconv/immunedeconv
  "data.table",       # File writing (.gz)
  "broom",            # Stats to tidy output
  "tidyr",            # Data wrangling 
  "stringr",          # String manipulation
  "fastDummies"      # One hot encoding/Dummy variables 
))

# Install Bioconductor packages #####
BiocManager::install(c(
  "limma",            # For normalization and QC
  "clusterProfiler",  # For functional enrichment analysis
  "org.Hs.eg.db",     # Human gene annotations
  "STRINGdb",         # For protein-protein interaction networks
  "impute",           # Missing value imputation
  "preprocessCore",   # Low-level Preprocessing and Normalization
  "GO.db",            # Gene Ontology (GO) Database
  "AnnotationDbi",    # Database Interface/Query Tool
  "GEOquery",         # Downloads Data from NCBI Geodataset
  "biomaRt",          # To retrieve ESEMBL IDs
  "EnhancedVolcano",  # To plot and visualize DEG results 
  "ComplexHeatmap",   # To visualize clusters in count data 
  "ImmuneSigR",       # To get immune markers present in the sig results enriched in WT
  "omnideconv/immunedeconv" # To get immune scoring 
))

# Load required libraries #####
library(WGCNA)
library(ggplot2)
library(gplots)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(readr)
library(doParallel)
library(GEOquery)
library(biomaRt)
library(dplyr) 
library(limma)
library(GO.db)
library(GOplot)
library(karyoploteR)
library(Gviz)
library(gggenes)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(ggridges)
library(VennDiagram)
library(rtracklayer)
library(DOSE)
library(purrr)
library(ImmuneSigR)
library(tidyestimate)
library(remotes)
library(immunedeconv)
library(data.table)
library(tidyr)
library(broom)
library(stringr)
library(fastDummies)

# Increasing timeout cutoff 
options(timeout = 1200)

# Set working directory #####
setwd("/Users/luzsmac/Desktop/version_control/hTERT_RPE_1_GSE305810/hTERT_RPE_1_WTvsXISTKD_versions ")

if(!dir.exists("/Users/luzsmac/Desktop/version_control/hTERT_RPE_1_GSE305810/hTERT_RPE_1_WTvsXISTKD_versions /GSE305810")){
# Loading Data #####
getGEOSuppFiles("GSE305810")

}

# Extract sample information (metadata) #####
sample_info <- getGEO("GSE305810", GSEMatrix = TRUE)
sample_info <- pData(sample_info[[1]])
  
# List all files in the downloaded folder #####
files <- list.files("GSE305810", full.names = TRUE)

# Read the count data file (In RPKM) #####
count_data <- read.delim(files[1])

# Number of columns #####
last_col <- ncol(count_data)

# Filter rows where the sum of columns 5 through the end is > 0 #####
count_data <- count_data[rowSums(count_data[, 5:last_col], na.rm = TRUE ) >0, ] 

# Naming the rows #####                                                               
rownames(count_data) <- count_data$Gene

# Genes in count data#####
write.table(count_data$Gene, "Data/count_data_genes.txt")

# Cleaning sample info #####
sample_info_cleaned <- data.frame(
  sample_id = sample_info$title,
  treatment = gsub("^[^,]+,([^,]+),.+", "\\1" ,sample_info$title)
)
rownames(sample_info) = sample_info$sample_id

######### LIMMA########

# Subsetting to get the RNA-seq reads only #####
matrix <- count_data[ ,c("WT1", "WT2", "WT3","XISTAC1A5","XISTAC1B3","XISTAC2B2", "XISTAC2B4")]  

# Metadata #####
gene_metadata <- count_data[,c("Gene", "Chromosome", "Start", "End")]

group_names = factor(c(
  rep("WT",3),
  rep("XIST_KD",4)
))

# Log Transform #####
log_matrix <- as.matrix(log2(matrix+1))

# Construct Matrix #####
design <- model.matrix(~0 + group_names)
colnames(design) <- levels(group_names)

# Fitting #####
fit <- lmFit(log_matrix, design, genes = gene_metadata)

# Contrast #####
contrast_matrix <- makeContrasts(
  WT_vs_XISTKD = WT-XIST_KD,
  levels = design
)

# Fitting with contrasts #####
fit2 <- contrasts.fit(fit, contrast_matrix)

# Fitting Ebayes #####
fit3 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# Summary of results #####
results <- decideTests(fit3)

# Getting chrom metadata #####
chrom_metadata <- count_data$Chromosome

# Results for Limma #####
all_results <- topTable(
  fit3,
  coef = "WT_vs_XISTKD",
  sort.by = "none",
  n = Inf,
  adjust.method = "BH"
)
#Saving the limma results 
write.table(all_results, file="Data/all_results_limma_WT_vs_XISTKD.txt", sep ="/t")

# Plotting the differential analysis results 
all_results$status <- "Not significant"

all_results$status[all_results$adj.P.Val < 0.05 & all_results$logFC > 0] <- "Upregulated"
all_results$status[all_results$adj.P.Val < 0.05 & all_results$logFC < 0] <- "Downregulated"

counts <- table(all_results$status)

png("plots/DEG_counts_text.png",
    width = 1200, height = 800, res = 150)

# blank canvas with nicer limits
plot(1, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "",
     main = "Differential Expression Summary GSE305810 WT vs XIST KD hTERT RPE-1")

# optional subtitle-style total
total_genes <- sum(counts)

text(0.5, 0.85,
     paste0("Total genes analyzed: ", total_genes),
     cex = 1.4, col = "black")

# Upregulated
text(0.5, 0.65,
     paste0("Upregulated: ", counts["Upregulated"]),
     col = "#D62728", cex = 2.2, font = 2)

# Downregulated
text(0.5, 0.45,
     paste0("Downregulated: ", counts["Downregulated"]),
     col = "#1F77B4", cex = 2.2, font = 2)

# Not significant
text(0.5, 0.25,
     paste0("Not significant: ", counts["Not significant"]),
     col = "gray40", cex = 2.0, font = 2)

dev.off()


# Transforming into X, Y or Other #####
custom_chrom <- function(x){
  if(x == "X"){
    return("X")
  }
  else if(x =="Y"){
    return("Y")
  }
  else{
    return("Other")
  }
}

all_results['Chrom'] <- sapply(chrom_metadata, custom_chrom)

# Plotting DEGS####
png("plots/Differentially_Enhaced_Volcano_GSE305810.png", width = 1000, height =500, res =100 )
EnhancedVolcano(
  all_results,
  x = 'logFC',
  y = 'adj.P.Val',
  lab = rownames(all_results),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 2,
  labSize = 3,
  title = 'Differentially Expressed Genes in GSE305810',
  subtitle = 'WT vs XISTKD using limma'
)
dev.off()

# Color coding by chrom ###
keyvals.col <- ifelse(all_results$Chrom == "X", "forestgreen","lightcoral")

names(keyvals.col) <- all_results$Chrom


# Odering ddf by X so that the volcano plot will show the X color on top 
all_results <- all_results[order(all_results$Chrom == "X"),]
keyvals.col <- sort(keyvals.col, decreasing= TRUE)

# Plotting, chrom by color ####                                
png("plots/Differentially_Enhaced_Volcano_Colored_by_Chrom_GSE305810.png", width = 1000, height =500, res =100)
p <-EnhancedVolcano(
  all_results,
  x = 'logFC',
  y = 'adj.P.Val',
  lab = rownames(all_results),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 2,
  labSize = 3,
  colCustom = keyvals.col,
  shape = 6,
  colAlpha = 1,
  title = 'Differentially Expressed Genes in GSE305810 Colored by Chromosome',
  subtitle = 'WT vs XISTKD using limma'
)
p + theme(plot.title = element_text(hjust = 0))
dev.off()

# Plotting the differential analysis results 
all_results_x_chrom <- subset(all_results, Chrom == "X")

all_results_x_chrom$status <- "Not significant"

# Upregulated (significant + positive logFC)
all_results_x_chrom$status[
  all_results_x_chrom$adj.P.Val < 0.05 & all_results_x_chrom$logFC > 0
] <- "Upregulated (FDR + logFC)"

# Downregulated (significant + negative logFC)
all_results_x_chrom$status[
  all_results_x_chrom$adj.P.Val < 0.05 & all_results_x_chrom$logFC < 0
] <- "Downregulated (FDR + logFC)"

counts <- table(all_results_x_chrom$status)
total_genes <- sum(counts)

png("plots/DEG_counts_text_X_chom.png",
    width = 1200, height = 800, res = 150)
plot(1, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "",
     main = "Differential Expression GSE305810 WT vs XIST KD hTERT RPE-1 (X Chrom)")

text(0.5, 0.85,
     paste0("Total genes analysed: ", total_genes),
     cex = 1.4, col = "black")

# Upregulated
text(0.5, 0.65,
     paste0("Upregulated (FDR + logFC): ", counts["Upregulated (FDR + logFC)"]),
     col = "#D62728", cex = 2.0, font = 2)

# Downregulated
text(0.5, 0.45,
     paste0("Downregulated (FDR + logFC): ", counts["Downregulated (FDR + logFC)"]),
     col = "#1F77B4", cex = 2.0, font = 2)

# Not significant
text(0.5, 0.25,
     paste0("Not significant: ", counts["Not significant"]),
     col = "gray40", cex = 2.0, font = 2)

dev.off()

# Getting significant results #####
sig_results <- all_results[which(all_results$adj.P.Val<0.05), ]

# Subsetting for sig on count data #####
count_data_sig <- count_data[rownames(sig_results), ]

# Saving sig count_data #####
write.table(count_data_sig, "Data/count_data_sig.txt", sep ="\t")

# Matrix sig #####
matrix_sig <- matrix[rownames(sig_results), ]

matrix_sig_scaled <- t(scale(t(matrix_sig)))

# Create annotation dataframe for samples #####
annotation_col <- data.frame(
  Group = sample_info_cleaned$treatment,
  row.names = sample_info_cleaned$sample_id
)

# Hierarchical tree of the samples #####
hc_samples <- hclust(as.dist(1-cor(matrix_sig_scaled ,
                                   method = "pearson")), method = "complete")
# Getting dendrogram
sampleTree = as.dendrogram(hc_samples, method = average)

# Plotting Samples #####
png("plots/Sample_clustering_GSE305810.png", width = 800, height = 600, res = 100)

plot(sampleTree,
     main = "Sample Clustering in GSE305810 dataset ")

dev.off()

# Hierarchical tree of genes #####
hc_genes <- hclust(as.dist(1-cor(t(matrix_sig_scaled) ,
                                 method ="pearson")), method = "complete")
# Getting dendrogram #####
geneTree = as.dendrogram(hc_genes, method = average)

# Plotting Samples #####
png("plots/Gene_clustering_GSE305810.png", width = 800, height = 600, res = 100)

plot(geneTree,
     main = "Gene Clustering in GSE305810 dataset ",
     leaflab = "none",
     ylab = "Height")

dev.off()

# Heatmap #####
png("plots/Heatmap_GSE305810.png", width = 800, height = 600, res = 100)

heatmap.2(as.matrix(matrix_sig),
          Rowv = as.dendrogram(hc_genes),
          Colv = as.dendrogram(hc_samples),
          col = redgreen(100),
          scale ="row",
          margins = c(7,7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap of GSE305810 dataset",
          trace = "none"
)

dev.off()

# Clustering #####
n_clusters <- 2

gene_clusters <- hc_genes |>
  cutree( k = n_clusters) |>
  sort()

# Cluster dataframe #####
cluster_annotation<- data.frame(
  SYMBOL = names(gene_clusters),
  Cluster = gene_clusters
)

# Reorder genes by cluster assignment #####
deg_scaled_ordered<- matrix_sig_scaled %>% 
                                  data.frame() %>%  
                                  filter(rownames(.) %in% cluster_annotation$SYMBOL)%>% 
                                  mutate(
                                    WT_mean= rowMeans(pick("WT1", "WT2", "WT3"),na.rm = TRUE),
                                    KD_mean = rowMeans(pick("XISTAC1A5","XISTAC1B3","XISTAC2B2", "XISTAC2B4"),na.rm = TRUE)) 

deg_scaled_ordered <-data.frame(deg_scaled_ordered[cluster_annotation$SYMBOL,], cluster_annotation$Cluster)


# Calculate gap positions (where clusters change) #####
gap_positions <- which(diff(cluster_annotation$Cluster) != 0)

columns_to_exclude_in_matrix<- "WT_mean|KD_mean|cluster_annotation.Cluster"
num_cols_to_exclude_matrix <-grep(columns_to_exclude_in_matrix, colnames(deg_scaled_ordered))
# Create annotated heatmap with gaps between clusters #####
png("plots/Heatmap_GSE305810_gapped.png", width = 800, height = 600, res = 100)
pheatmap(
  deg_scaled_ordered[,-num_cols_to_exclude_matrix],
  cluster_rows = FALSE,  
  cluster_cols = hc_samples,
  annotation_row = NA,
  annotation_col = NA,
  gaps_row = gap_positions,  
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "DEG Hierarchical Clusters (Separated by Gaps)",
  fontsize = 10,
  filename = "hierarchical_clustering_DEGs_gapped.png",
  width = 10,
  height = 8
)
dev.off()

# Enrichment analysis
# Convert gene symbols to Entrez IDs 
gene_entrez <- bitr(cluster_annotation$SYMBOL, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Setting the row names to the symbols #####
row.names(gene_entrez) <- gene_entrez$SYMBOL

# Saving entrez of clusters 
write.table(gene_entrez, file="Data/All_Clusters_entrez.txt", sep="\t")

#Getting the missing symbols #####
missing_genes <- setdiff(cluster_annotation$SYMBOL, gene_entrez$SYMBOL)

# Saving Symbols not mapped for documentation purposes #####
write.table(missing_genes, "Data/Genes_not_mapped_to_ENTREZID.txt")

# subset for cluster 1 #####
cluster_1 <- subset(cluster_annotation, Cluster == 1)
cluster_2 <- subset(cluster_annotation, Cluster == 2)

gene_entrez_cluster_1 <- gene_entrez[cluster_1$SYMBOL, ]
gene_entrez_cluster_2 <- gene_entrez[cluster_2$SYMBOL, ]

# Perform GO enrichment for WT-associated module #####
go_enrichment_cluster_1 <- enrichGO(
  gene = gene_entrez_cluster_1$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",              
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

png("plots/kegg_cluster1.png", width = 800, height = 600, res = 100)
go_enrichment_kegg_cluster1 <- enrichKEGG(
  gene = gene_entrez_cluster_1$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.02,
  qvalueCutoff = 0.05
)
dev.off()

go_enrichment_cluster_2 <- enrichGO(
  gene = gene_entrez_cluster_2$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",              
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.02,
  qvalueCutoff = 0.05,
  readable = TRUE
)
count_data[,"SYMBOL"] <- count_data[,"Gene"]

#GO Enrichment Barplot #####
png("plots/GO_enrichment_barplot_cluster_1_GSE305810.png",
    width = 2000, height = 4000, res = 100)
barplot(arrange(go_enrichment_cluster_1, qvalue), 
        showCategory = 100,
        title = "GO Enrichment in Cluster 1 WT hTERT RPE-1 cells (GSE305810)")
dev.off()


png("plots/GO_enrichment_barplot_cluster_2_GSE305810.png",
    width = 2000, height = 4000, res = 100)
barplot(arrange(go_enrichment_cluster_2, qvalue), 
        showCategory = 100,
        title = "GO Enrichment in Cluster 2 XISTKD hTERT RPE-1 cells (GSE305810) ")
dev.off()

# Entrez ID for count_data 
count_data_gene_entrez <- bitr(count_data$SYMBOL, 
                               fromType = "SYMBOL",
                               toType = "ENTREZID", 
                               OrgDb = org.Hs.eg.db)

# Check for NA and duplicates 0.01% not mapped 
count_data_gene_entrez<- count_data_gene_entrez[!is.na(count_data_gene_entrez$ENTREZID),]
count_data_gene_entrez<- count_data_gene_entrez[!duplicated(count_data_gene_entrez$ENTREZID),]
count_data_gene_entrez<- count_data_gene_entrez[!duplicated(count_data_gene_entrez$SYMBOL),]

# Rename rows on gee entrez count data
rownames(count_data_gene_entrez) <- count_data_gene_entrez$SYMBOL

# Prepare for enrichment analysis 
count_data_gene_entrez <- cbind(count_data_gene_entrez[rownames(count_data_gene_entrez),], all_results[rownames(count_data_gene_entrez), ]) 
gene_entrez_t <- count_data_gene_entrez$t
names(gene_entrez_t) <- count_data_gene_entrez$ENTREZID

# Gse enrichment 
gsea_do <- gseDO(
  geneList = sort(gene_entrez_t, decreasing= TRUE),
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# Filtering by sig 
sig_gsea <- filter(gsea_do, p.adjust < 0.05)

# Plotting the gsea
png("plots/gsea_do_enrichment.png", width = 2000,
    height = 4000,res = 100)
ridgeplot(sig_gsea, showCategory = 35, orderBy= "NES", fill= "p.adjust")+ 
  scale_fill_continuous(low="#FF0000CC", high="#3182bdCC") +
  ggtitle("GSEA Results in GSE305810 WT vs XIST KD hTERT RPE-1") +
  xlab("NES")+
  theme(text = element_text(size = 26)) +
  theme(axis.text.y = element_text(size = 26))

dev.off()


# Df of differentially expressed genes and chroms #####
sig_genes_df <- data.frame(
  Symbol = rownames(count_data_sig),
  Chrom = count_data_sig$Chromosome
)

possibly(function(){mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   host = "may2021.archive.ensembl.org")

# Retrieve coordinates #####
sig_genes_coordinates <- getBM(attributes = c("external_gene_name", "chromosome_name", "strand",
                                              "start_position", "end_position"),
                               filters = c("external_gene_name", "chromosome_name"),
                               values = list(sig_genes_df$Symbol, c(1:22, "X", "Y", "MT")),
                               mart = mart)

# Formatting to convert to granges #####
sig_genes_coordinates$chromosome_name <- paste0("chr", sig_genes_coordinates$chromosome_name)

# Formating strands #####
formatting_strands <- function(x){
  if(x==1){
    return("+")
  }
  else if (x==-1){
    return("-")
  }
  else{
    return("*")
  }
}

# Calling formatting_strands #####
sig_genes_coordinates$strand<-lapply(sig_genes_coordinates$strand, formatting_strands)
# Formatting 
sig_genes_coordinates$strand<-as.character(sig_genes_coordinates$strand)
# Saving coordinates of sig genes 
write.table(sig_genes_coordinates, file="Data/sig_genes_coordinates.txt", sep="\t")
return(sig_genes_coordinates)
},
otherwise=(
  sig_genes_coordinates<-read.table("Data/sig_genes_coordinates.txt")
  ))

# Df of unique genes #####                        
sig_genes_coordinates_unique <- sig_genes_coordinates[!duplicated(sig_genes_coordinates$external_gene_name), ]

# Indices of duplicated genes #####
duplicated_indices <- which(duplicated(sig_genes_coordinates$external_gene_name))

# Df of duplicated genes #####
df_duplicated <- sig_genes_coordinates[duplicated_indices, ]

# Define chrom order #####
chromosome_order <- factor(paste0("chr", c(1:22, "X", "Y")))

# Order df #####
sig_genes_coordinates_unique$chromosome_name <- factor(
  sig_genes_coordinates_unique$chromosome_name, 
  levels = chromosome_order
)
# Plotting genes and coordinates #####
png("plots/genes_coordinates.png", width = 4500, height = max(2000, 23 * 250), res = 300)
ggplot(sig_genes_coordinates_unique, aes(xmin = start_position, xmax = end_position,
                                         y= chromosome_name)) +
  geom_gene_arrow() +
  facet_wrap(~ chromosome_name, ncol = 1) +
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    strip.text.y = element_blank()
  )+
  ggtitle("Differentially Expressed Genes In WT vs XISTKD In GSE305810 By Chromosome")

dev.off()

# subsetting and saving data with X chrom #####
sig_genes_coordinates_unique_X_chrom <- filter(data.frame(sig_genes_coordinates_unique), chromosome_name == "chrX")

# Formatting 
sig_genes_coordinates_unique_X_chrom[,'strand' ] <- as.character(sig_genes_coordinates_unique_X_chrom[,'strand' ])

# Saving 
write.table(sig_genes_coordinates_unique_X_chrom , file = "Data/sig_genes_coordinates_unique_X_chrom.txt",  col.names=TRUE, row.names= FALSE)


# Subseting key genes 
sig_results <- subset(all_results, adj.P.Val< 0.05) 
vec_MIR <- grep("MIR", rownames(sig_results), value = TRUE)
vec_IL <- grep("IL", rownames(sig_results), value = TRUE)

possibly(function(){
# List of attributes available in mart
attributes_df <- listAttributes(mart)

# Retriving gene information 
gene_function1 <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name","description", "goslim_goa_description","definition_1006"),
  values = rownames(sig_results),
  mart = mart
)
data.table::fwrite(gene_function1,"Data/gene_function1.txt.gz", sep="\t")


gene_function2 <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name","namespace_1003", "biogrid", "reactome"),
  values = rownames(sig_results),
  mart = mart
)
write.table(gene_function2, file="Data/gene_function2.txt", sep="\t")
return(attributes_df, gene_function1, gene_function2)},
otherwise=(
  gene_function1 <- fread("Data/gene_function1.txt.gz", sep="\t")
  
))

# The following block of code will result in 

# Df of cluster 1 enrichment 
WT_enrichment <- as.data.frame(go_enrichment_cluster_1)
############## Delete?

# Subsetting the gene_function1 df with genes of interest from cluster2
gene_function1_subset_WT <- subset(gene_function1, external_gene_name %in% cluster_1$SYMBOL)

# Building the search term of functions 
search_terms_function <- "metabolism|immunity|autoimmunity|inflammation|pro-inflammatory|cytok|defense|bacteria|interleukin|immune"

# Retrieving the index that match the search terms 
index_description_of_interest_WT_1006 <- grep(search_terms_function, gene_function1_subset_WT$definition_1006) 
index_description_of_interest_WT_goa_description <- grep(search_terms_function, gene_function1_subset_WT$goslim_goa_description) 

# Finding unique indices in genes and descrition of interest in 1006 and goslim 
unique_index_of_interest <- unique(index_description_of_interest_WT_1006,index_description_of_interest_WT_goa_description)

# Gene function df of unique indices 
gene_function1_subset_description_of_interest<- gene_function1_subset_WT[unique_index_of_interest,  ]

# Retrieve sig results df for genes of interest based on 1006 function and goslim goa description of cluster 1 genes that match immune descriptions 
sig_results_interest_immune_related_description <- sig_results[unique(gene_function1_subset_description_of_interest$external_gene_name) , ]

# Initialize de gsea_do object
gsea_do_df_all_results <- as.data.frame(gsea_do)

# Retrieve markers 
all_markers <- Get_Markers()

all_markers_df <- data.frame(Marker = names(unlist(all_markers)),
                             Genes = unname(unlist(all_markers)))

# Subset Markers in cluster 1 (WT)
all_markers_df_cluster_WT <- subset(all_markers_df, Genes %in% cluster_1$SYMBOL)

# Retrieve stats for markers
all_markers_df_cluster_WT_limma_stats <- data.frame(all_markers_df_cluster_WT, sig_results[all_markers_df_cluster_WT$Genes, c(1,5)])

# Clean list of markers 
all_markers_df_cluster_WT_limma_stats$Markers_cleaned<-gsub("^(.+)_PMID.+","\\1",all_markers_df_cluster_WT$Marker)

all_markers_df_cluster_WT_limma_stats$Markers_cleaned<-gsub(" ","_",all_markers_df_cluster_WT_limma_stats$Markers_cleaned )
all_markers_df_cluster_WT_limma_stats$Markers_cleaned<-gsub("/|-|,","_",all_markers_df_cluster_WT_limma_stats$Markers_cleaned )
all_markers_df_cluster_WT_limma_stats$Markers_cleaned<-gsub("__","_",all_markers_df_cluster_WT_limma_stats$Markers_cleaned )


# Subset Markers in cluster 1 (XISTKD)
all_markers_df_cluster_XISTKD <- subset(all_markers_df, Genes %in% cluster_2$SYMBOL)

# Retrieve stats for markers
all_markers_df_cluster_XISTKD_limma_stats <- data.frame(all_markers_df_cluster_XISTKD, sig_results[all_markers_df_cluster_XISTKD$Genes, c(1,5)])

# Clean list of markers 
all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned<-gsub("^(.+)_PMID.+","\\1",all_markers_df_cluster_XISTKD$Marker)

all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned<-gsub(" ","_",all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned )
all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned<-gsub("/|-|,","_",all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned )
all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned<-gsub("__","_",all_markers_df_cluster_XISTKD_limma_stats$Markers_cleaned )

# Adding 
# Using immunedeconv to test for immune score 
methods<- c(
  "quantiseq",
  "mcp_counter",
  "xcell",
  "epic",
  "abis",
  "estimate")

immune_test_list <- list()
immune_test_names <- c()

for(i in seq_along(methods)){
  # If the method is not Epic evaluate 
  if(methods[i] != "epic"){
  immune_test_list[[i]] <-immunedeconv::deconvolute(matrix, method= methods[i], arrays= FALSE)
  immune_test_names[i]<- methods[i]}
  
  # Epic has a specific hyper-parameter so it must be passed in 
  else if(methods[i] == "epic"){
    immune_test_list[[i]] <-immunedeconv::deconvolute(matrix, method= methods[i], arrays= FALSE,tumor=FALSE)
    immune_test_names[i]<- methods[i]}}

# Preparing for statistical testing return a list of all the dfs with the test info
all_immune_test_list_with_metadata<- list()
retrieving_stats<-function(df, name){
  df %>%
  data.frame()%>%
    mutate(Test= rep(c(name), length.out = nrow(df)))}

# Calling the stats function 
all_immune_test_list_with_metadata<-Map(retrieving_stats,immune_test_list,immune_test_names)

# Binding all the results for each immune scoring test 
all_immune_test_df_with_metadata<-c()
for(i in seq_along(all_immune_test_list_with_metadata)){
  all_immune_test_df_with_metadata<- rbind(all_immune_test_df_with_metadata,
                                           all_immune_test_list_with_metadata[[i]])}
# Creating rownames for matrix/df
col1<- all_immune_test_df_with_metadata["cell_type"]
test<- all_immune_test_df_with_metadata["Test"]

# Creating rownames with cell_type and test info
rownames_immune_tests<-unname(unlist(Map(paste0, col1,"_",test)))
rownames_immune_tests<- gsub(" ", "_",rownames_immune_tests)
rownames( all_immune_test_df_with_metadata)<- rownames_immune_tests
all_immune_test_df_with_metadata<- all_immune_test_df_with_metadata[rowSums(all_immune_test_df_with_metadata[, -c(1,9)], na.rm = TRUE ) >0, ] 

## Performing Wilcoxon Rank-Sum Test with greater 

# Retrieving cell_type + test names
col_1 <-rownames(all_immune_test_df_with_metadata)

# Retrieving sample names 
sample_names <- colnames(all_immune_test_df_with_metadata[ ,-c(1,9)])

# Running Wilcoxon Rank Sum, greater for each celltype+test WT vs XISTAC

all_immune_test_df_wilcox_results<-all_immune_test_df_with_metadata %>%
  select(-c(1,9)) %>%
  t() %>%
  data.frame() %>%
  rename_with(function(x) col_1) %>%
  mutate(group = c("WT", "WT", "WT", "XISTAC", "XISTAC", "XISTAC", "XISTAC"),
         sample = sample_names) %>%
  pivot_longer(cols = -c(group, sample), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  group_split() %>% # group and split by variable
  set_names(map_chr(., ~ unique(.x$variable))) %>%  # set names to retain varIable name 
  map(~ wilcox.test(value ~ group, data = ., alternative = "greater")) %>% # wilcoxon with greater 
  map_dfr(~ tidy(.), .id = "variable") %>%
  mutate(p.adjust=p.adjust(p.value, method = "BH")) %>% # FDR corection
  rename(w_statistic =statistic)


# Adding r effect and direction 
all_immune_test_fold_change_effect_size<-all_immune_test_df_wilcox_results %>%
  left_join(
    all_immune_test_df_with_metadata %>%
      select(cell_type, Test) %>%
      rename(variable = cell_type, method = Test),
    by = "variable"
  ) %>%
  mutate(
    direction = ifelse(w_statistic > 6, "WT > XISTAC", "XISTAC > WT"),
    Z        = (w_statistic - (3*4/2)) / sqrt(3*4*(3+4+1)/12),
    r_effect = abs(Z) / sqrt(3+4),  
    magnitude = case_when(
      abs(r_effect) >= 0.50 ~ "large",
      abs(r_effect) >= 0.30 ~ "medium",
      abs(r_effect) >= 0.10 ~ "small",
      TRUE                  ~ "negligible"))
  
# Adding cell_type to the Wilcoxon test with effect size 
# Reordering all_immune_test_df_with_metadata
all_immune_test_df_with_metadata<-all_immune_test_df_with_metadata[all_immune_test_fold_change_effect_size$variable, ]

# Testing for row order and adding cell info
if(all.equal(all_immune_test_fold_change_effect_size$variable, rownames(all_immune_test_df_with_metadata))){
  all_immune_test_fold_change_effect_size["cell_type"] <- all_immune_test_df_with_metadata$cell_type}

# Retrieving dfs, splitting by direction 
unique_cell_types_deconvolution <- unique(all_immune_test_fold_change_effect_size["cell_type"])

# Cells that are statistically more significant in WT 
cell_deconvolution_upregulated_in_WT<- subset(all_immune_test_fold_change_effect_size, direction == "WT > XISTAC")

# Cells that are statistically more significant in XIST KD
cell_deconvolution_upregulated_in_XISTKD<- subset(all_immune_test_fold_change_effect_size, direction == "XISTAC > WT")

# Retrieving unique cells in each df split above 
unique_cells_WT<-cell_deconvolution_upregulated_in_WT %>%
                              select(cell_type) %>%
                              unique() %>%
                              data.frame() 


unique_cells_XISTKD<-cell_deconvolution_upregulated_in_XISTKD %>%
                            select(cell_type) %>%
                            unique() %>%
                            data.frame() 




# Separate Wider to One hot encode 
pathway_genes_wider_WT <- WT_enrichment %>% separate_wider_delim(geneID, delim = "/", names_sep = "",
                                                       too_few="align_start", names_repair = "universal" )  

# Retrieving cols to one hot encode 
gene_cols_WT <- paste0("geneID", seq_along(1:47))

# Retrieving list of unique genes 

all_genes_pathway_vector_WT<- c()
removing_dups<-function(x){
  all_genes_pathway_vector_WT<-c(all_genes_pathway_vector_WT,x)
}
vector_of_all_genes_in_pathway_WT<-lapply(pathway_genes_wider_WT[,gene_cols_WT], function(x) removing_dups(x))
# Retrieve unique gene names 
vector_unique_genes_pathway_WT<-unique(unlist(vector_of_all_genes_in_pathway_WT))
# Remove duplicates 
vector_unique_genes_pathway_WT<-sort(vector_unique_genes_pathway_WT[which(!is.na(vector_unique_genes_pathway_WT))])

# Creating dummies for each gene column 
pathway_genes_one_hot_encoded_WT <- dummy_cols(pathway_genes_wider_WT, select_columns = gene_cols_WT, remove_selected_columns = TRUE)

library(tibble)
##Setting rownames 
pathway_genes_one_hot_encoded_WT <- pathway_genes_one_hot_encoded_WT %>%
                                    mutate(Description_copy = Description) %>% 
                                    column_to_rownames(var = "Description_copy")

# Remove NA cols after encoding 
pathway_genes_one_hot_encoded_WT <- pathway_genes_one_hot_encoded_WT[ ,  -(grep("geneID\\d+_NA$",names(pathway_genes_one_hot_encoded_WT)))]



## Write an assertion test just to check for correct one hot encoding 
#test_1<- pathway_genes_one_hot_encoded[ "purine ribonucleotide metabolic process",c(1:11,which(pathway_genes_one_hot_encoded["purine ribonucleotide metabolic process", ]== 1)) ]
#test-2 <- pathway_genes_one_hot_encoded[]
#stopifnot(
#all.equal(test_1, test_2)
#)
 ###########################

# Gene cols one hot encoded 

# Renaming one hot encoded genes
pathway_genes_one_hot_encoded_cols_only_WT<- gsub("[geneID0-9]*_([A-Za-z0-9]+)", "\\1",names(pathway_genes_one_hot_encoded_WT)[12:ncol(pathway_genes_one_hot_encoded_WT)])
last_col<-dim(pathway_genes_one_hot_encoded_WT)[2]
names(pathway_genes_one_hot_encoded_WT)[12:last_col] <- pathway_genes_one_hot_encoded_cols_only_WT
# Dropping non gene cols 
genes_only_pathway_onehot_encoded <-pathway_genes_one_hot_encoded_WT[12: last_col]

#Gene info only 
non_gene_data <-pathway_genes_one_hot_encoded_WT[1:11]

# Order colnames in alphabetical order
genes_only_pathway_onehot_encoded <-genes_only_pathway_onehot_encoded[order(colnames(genes_only_pathway_onehot_encoded))]

## Index check!!!!
all.equal(rownames(genes_only_pathway_onehot_encoded),rownames(non_gene_data))
all.equal(rownames(pathway_genes_one_hot_encoded_WT), rownames(non_gene_data))

# Consolidate duplicated cols 
index_list_TRUE <- list() # To store all the indices that have val ==1
True_at_i<-c() # To store the indices per consolidated col that have val==1
target <-1
last_col_df<- ncol(genes_only_pathway_onehot_encoded)
unique_names_check <- vector_unique_genes_pathway_WT
length_unique_names <- length(unique_names_check) # To keep track of how many cols I retrieve 
unique_names <- list()
k=1 #Counter for unique cols
for(i in 1:last_col_df){
  
 if(i==1){
   index<- which(genes_only_pathway_onehot_encoded[ , i]== 1) 
   True_at_i <- c(True_at_i,index) # Adds col 1 to the per unique col vector 
 }

 if(i!=1){
   index_at_i<-i
   index_at_i_minus <-(index_at_i-1)
   # Evaluates consecutive colnames
   boolean_index<-(names(genes_only_pathway_onehot_encoded)[index_at_i]==names(genes_only_pathway_onehot_encoded)[index_at_i_minus])
   
   # If consecutive names are the same add to per unique col vector
   if( boolean_index==TRUE){
     True_at_i <-c(True_at_i, which(genes_only_pathway_onehot_encoded[ , index_at_i]== target))
   }
   # If consecutive names are not the same add to list of unique cols with index at which val ==1 
   if(boolean_index==FALSE){
     if(k !=length_unique_names){
      index_list_TRUE[[k]]<-True_at_i # Add first
      unique_names[[k]]<-unique_names_check[k] # Keeping track of the unique colname order
      True_at_i <-c() # Restart vector per col
      True_at_i <-c(True_at_i, which(genes_only_pathway_onehot_encoded[ , index_at_i]== target)) # Evaluate then add 
      k=k+1 # Add to the index in unique colnames
     }
     if(k ==length_unique_names){ # if last unique col add
       index_list_TRUE[[k]]<- True_at_i
       unique_names[[k]]<-unique_names_check[k] # Keeping track of the unique colname order
     }
   }
 }
}

# Retrieve results
unique_names<-unlist(unique_names)
# Check lengths of unique names found by the unique function and the iterative function
stopifnot(length(index_list_TRUE)==length_unique_names & (length_unique_names== length(unique_names)))

# Checking for order of names of unique cols that were merged 
stopifnot(all.equal(unique_names, vector_unique_genes_pathway_WT))

setdiff(unique_names,vector_unique_genes_pathway_WT)

# Checking for proper consolidation 
target<- 1
cbinded_vectors<- c()
num_rows <- nrow(pathway_genes_one_hot_encoded_WT)
vector<- rep(0,num_rows)

for(i in seq_along(1:length(index_list_TRUE))){
  current_cols_index_TRUE<-unlist(index_list_TRUE[i]) # Retrieve the indices that are TRUE for each gene
  vector[current_cols_index_TRUE]<- target # Adding the the target to the indices for each col
  cbinded_vectors<-cbind(cbinded_vectors,vector) # bind all the cols together 
}
df_not_dup_cols<-data.frame(cbinded_vectors) # Convert to df 
if(dim(df_not_dup_cols)[2]== length(unique_names)){ # Check for the correct # of cols 
names(df_not_dup_cols)<-unique_names}

# Creating random 20 random numbers from 1 to 329 for each column. 
random_nums_to_test_merge<-sample(ncol(df_not_dup_cols), size = 20, replace = FALSE) # Random sampling to check for proper consolidation 
test_merge<-function(ran){
all.equal(sort(which(df_not_dup_cols[ ,ran]==1)),sort(unlist(index_list_TRUE[ran])))
}

results<-lapply(random_nums_to_test_merge,function(x) test_merge(x)) # Retrieve results of the consolidation test 


stopifnot(length(which(results==FALSE))==0) # Stop if any tests came out not equal/FALSE 

# Index check 
all.equal(non_gene_data$Description, pathway_genes_one_hot_encoded_WT$Description)

# Merge df with pathway info
WT_cleaned_pathway_genes_df<-data.frame(non_gene_data, df_not_dup_cols)

# Making sure all the genes are present
setdiff( unique_names , cluster_1$SYMBOL)
setdiff(unique(all_markers_df_cluster_WT_limma_stats$Genes), cluster_1$SYMBOL)

#Renaming for joining 
all_results[ "Genes"]<- rownames(all_results)
names(count_data)[1]<- "Genes"

### Not able to see deubiquination pathway !!!!!!!!!!!!!!!!########

# Aggregating pathway,marker, limma stats, count data 
markers_pathway_WT <- WT_cleaned_pathway_genes_df %>%
  pivot_longer(
    cols= c(12:363),
    names_to= "Genes",
    values_to= "value") %>%
    filter(value ==1) %>%
    full_join(all_markers_df_cluster_WT_limma_stats[c(1,2,5)],by="Genes") %>%
    left_join(all_results[ , -(8)], by= "Genes") %>%
    left_join(count_data[ , 1:7], by = "Genes") 

# Quick look at the available pathways for WT   
unique_pathway_WT<-data.frame(unique(markers_pathway_WT$Description))

# Retrieve the pathways of interest
pathways_of_interest<- c()
search_order<-c("hypoxia", "protein", "deubiquitination","apopt", "membrane fusion", "phagy")
for(term in search_order){
  current_search<-grep(term, markers_pathway_WT$Description, value=TRUE)
  pathways_of_interest<-c(pathways_of_interest, current_search)
}

# Retrieving the mean expr of WT cells 
markers_pathway_WT<-markers_pathway_WT %>%
                        mutate(
                        WT_mean= rowMeans(pick("WT1", "WT2", "WT3"),na.rm = TRUE))

# Subset and order selected pathways 
markers_selected_pathways_WT<-subset(markers_pathway_WT, Description %in% pathways_of_interest)
markers_selected_pathways_WT <- markers_selected_pathways_WT %>% 
  mutate(Description = factor(Description, levels = unique(pathways_of_interest)))

# Plotting Pathway and gene info colored by avg exp in WT 
n_genes <- length(unique(markers_selected_pathways_WT$Genes))
n_pathways <- length(unique(markers_selected_pathways_WT$Description))

px_per_gene <- 90
px_per_pathway <- 40
res_val <- 300

plot_width  <- max(8000, n_genes * px_per_gene)
plot_height <- max(6000, n_pathways * px_per_pathway)

png("plots/pathway_genes.png", width = plot_width, height = plot_height, res = res_val)

ggplot(markers_selected_pathways_WT, aes(Genes, Description, fill = WT_mean)) + 
  geom_tile(aes(width = 0.9, height = 0.9)) + 
  scale_fill_stepsn(n.breaks = 15, colours = c("#FFB6C1", "#DA70D6", "#4B0082")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = "Genes", y = "Pathway", fill = "WT_mean", 
       title="Avg WT expression of each gene in apoptosis and pyroptosis related pathways")+
  theme(plot.title = element_text(hjust = 0.5, size = 18))

dev.off()

# Gene and markers 
## Plotting Pathway and gene info colored by avg exp in WT 
n_genes <- length(unique(markers_selected_pathways_WT$Genes))
n_markers<- length(unique(markers_selected_pathways_WT$Marker))

px_per_gene <- 90
px_per_markers <- 40
res_val <- 100

plot_width  <- max(8000, n_genes * px_per_gene)
plot_height <- max(4000, n_markers * px_per_markers)

png("plots/marker_genes.png", width = plot_width, height = plot_height, res = res_val)

ggplot(markers_selected_pathways_WT, aes(Marker, Genes, fill = WT_mean)) + 
  geom_tile(aes(width = 0.9, height = 0.9)) + 
  scale_fill_stepsn(n.breaks = 15, colours = c("#FFB6C1", "#DA70D6", "#4B0082")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = "Genes", y = "Marker", fill = "WT_mean", 
       title="Avg WT expression of each gene in apoptosis and pyroptosis related markers")+
  theme(plot.title = element_text(hjust = 0.5, size = 18))

dev.off()

markers_genes_expr_WT <-markers_selected_pathways_WT %>%
  select(2,12,14,15,29) %>%
  pivot_longer(
    cols=c(3),
    values_to= "Markers"
  ) %>%
  select(c(1,2,3,4,6))
 
markers_genes_expr_WT_split<- markers_genes_expr_WT %>%
  group_split(Description, Markers_cleaned) %>%
  bind_rows() %>%
  as.data.frame() %>%
  mutate(General_pathway = case_when(
    str_detect(Description, "(?i)hypoxia") ~ "hypoxia",
    str_detect(Description, "(?i)protein") ~ "protein remodeling",
    str_detect(Description, "(?i)deubiquitination") ~ "deubiquitination",
    str_detect(Description, "(?i)apopt") ~ "apoptosis",
    str_detect(Description, "(?i)fusion") ~ "membrane fusion",
    str_detect(Description, "(?i)phagy") ~ "paghy"))


png("plots/genes_markers_general_pathway.png", width = 3000, height = 2000, res = res_val)
ggplot(markers_genes_expr_WT_split, aes(x = Genes, y = Markers_cleaned, fill = WT_mean)) +
    geom_tile(aes(width = 0.9, height = 0.9))+
    facet_wrap(~General_pathway)
dev.off()

length(unique(markers_genes_expr_WT_split$Description))
lapply(markers_genes_expr_WT_split, function(x) plotting_by_description(x))
## Add cyber scores 
# Any sig in the immune scores? look at those markers and find the genes, descriptions and pathways, disease
# Look at which genes are contributing to disease and sig enrichment pathways 
# Add anything to search terms? lipoprotein 
# Look which sig immune genes are in the same TAD 
# Map by coordinate 
# Description for each biomarker found to be sig 
# Look at info on reactive species, lactate pathways?
# APOE? lipids APOC1??
# Look at hormone marker? 
# Look at react to me 
# Co-expression network
# Cytokines 
# Back everything up
# reproduce for KD cluster!!! def and 

