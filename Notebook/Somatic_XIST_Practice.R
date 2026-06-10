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
  "tidyr"             # data wrangling 
))
remotes::install_github("omnideconv/immunedeconv")

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

# Increasing timeout cutoff 
options(timeout = 1200)

# Set working directory #####
setwd("/Users/luzsmac/Desktop/version_control/hTERT_RPE_1_GSE305810/hTERT_RPE_1_WTvsXISTKD_versions ")

if(!dir.exists("/Users/luzsmac/Desktop/version_control/hTERT_RPE_1_GSE305810/hTERT_RPE_1_WTvsXISTKD_versions /GSE305810")){
# Loading Data #####
getGEOSuppFiles("GSE305810")
  
# Extract sample information (metadata) #####
sample_info <- getGEO("GSE305810", GSEMatrix = TRUE)
sample_info <- pData(sample_info[[1]])
}
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
summary(results)
test <- as.data.frame(fit3)
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

# Retriving key genes 
phagy_index <- grep("phag", WT_enrichment$Description)
kinase_index <- grep("kinase", WT_enrichment$Description)
oxy_index <- grep("oxy", WT_enrichment$Description)

# Building the search term 
search_term <- "pag|kinase|oxygen|hypoxia|cata|substrate|ER|interleukin|stress|tyrosine|transferase|vacuole|steroid|bacteria "

# Retrieving genes that match the search term 
index_enrichment_retrieval <- grep(search_term, WT_enrichment$Description)
WT_enrichment_subset <- WT_enrichment[index_enrichment_retrieval, ]

# Retrieving gene names and cleaning 
WT_genes_of_interest <- WT_enrichment_subset$geneID
WT_genes_of_interest <- sapply(WT_genes_of_interest, function(x) gsub("/", ",", x))
WT_genes_of_interest <- strsplit(WT_genes_of_interest, split =",")
WT_genes_of_interest <- unlist(unname(WT_genes_of_interest))

# Subsetting the gene_function1 df with genes of interest from cluster2
gene_function1_subset <- subset(gene_function1, external_gene_name %in% WT_genes_of_interest)

# Building the search term of functions 
search_terms_function <- "metabolism|immunity|autoimmunity|inflammation|pro-inflammatory|cytok|defense|bacteria|interleukin|immune"

# Retrieving the index that match the search terms 
index_description_of_interest_WT_1006 <- grep(search_terms_function, gene_function1_subset$definition_1006) 
index_description_of_interest_WT_goa_description <- grep(search_terms_function, gene_function1_subset$goslim_goa_description) 

# Finding unique indices in genes and descrition of interest in 1006 and goslim 
unique_index_of_interest <- unique(index_description_of_interest_WT_1006,index_description_of_interest_WT_goa_description)

# Gene function df of unique indices 
gene_function1_subset_description_of_interest<- gene_function1_subset[unique_index_of_interest,  ]

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
all_markers_df_cluster_WT$Markers_cleaned<-gsub("^(.+)_PMID.+","\\1",all_markers_df_cluster_WT$Marker)

all_markers_df_cluster_WT$Markers_cleaned<-gsub(" ","_",all_markers_df_cluster_WT$Markers_cleaned )
all_markers_df_cluster_WT$Markers_cleaned<-gsub("/|-|,","_",all_markers_df_cluster_WT$Markers_cleaned )
all_markers_df_cluster_WT$Markers_cleaned<-gsub("__","_",all_markers_df_cluster_WT$Markers_cleaned )

# Using immunedeconv to test for immune score 
# 
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

#Extract each df 
for(i in seq_along(immune_test_list)){
  assign(paste0("immune_",immune_test_names[i],"_df"), data.frame(immune_test_list[i]))
}

# Preparing for statistical testing 
Group <- c("WT","WT","WT","XISTAC","XISTAC","XISTAC","XISTAC")
immune_test_list_stats<-list()
retrieving_stats<-function(df){
  col_1<-as.character(df[[1]])
  df %>%
  select(-1) %>%
  t() %>%
  data.frame() %>%
  rename_with(function(x) col_1) %>%
  mutate(Group = c("WT", "WT", "WT", "XISTAC", "XISTAC", "XISTAC", "XISTAC")) %>%
  pivot_longer(cols = -Group, names_to = "variable", values_to = "value") %>%
  
  # Grouping by variable
  group_by(variable) %>% 
  
  # T test between XISTKD and WT
    summarise(
      # Compares the groups against each other within each variable
      p_val = t.test(value ~ Group, data = pick(value, Group))$p.value, 
      .groups = "drop"
    )}
# Calling the stats function 
  immune_test_list_stats <-lapply(immune_test_list, function(x) retrieving_stats(x))

# Unpacking the results of the t test for each deconvolution method 
for(i in seq_along(immune_test_list_stats)){
  assign(paste0("immune_",immune_test_names[i],"_df_stats"), data.frame(immune_test_list_stats[i]))
}
# Normality test  
p_values<- c()
shapiro<-function(current_df){
  working_df<- current_df[ ,-1]
    p_values <-lapply(working_df, function(x) shapiro.test(x)$p.value)
  }
normality_test_list <-lapply(immune_test_list, function(x) shapiro(x))

# Iterate find test with normality , found test sets 2 and 6 to have normal cols 
for(i in seq_along(normality_test_list)){
  normal_found<-unlist(normality_test_list[[i]])>0.05
  if(any(normal_found==TRUE)){
    print(paste0(immune_test_names[[i]], " at index ", i, " in immune_test_names"))}
}
# Print which cols are normal in test sets 2 and 6
mcp_counter_normal_cols<- which(normality_test_list[[2]]>=0.05)
estimate_normal_cols <- which(normality_test_list[[6]]>=0.05)

install.packages("fitdistrplus")
library(fitdistrplus)

distrubutions<- c("norm", "lnorm", "exp", "pois", "cauchy", "gamma", "logis", "nbinom", "geom", "beta", "weibull")


descdist(df, discrete = FALSE) 
dictribution_results<-list()
retrieving_distribution_results<-function(df){
  df%>%
    select(-1)%>%
      summarise(across(all_of(my_columns), list(fitdist(across(all_of(my_columns))), na.rm = TRUE))
    )}
lapply(immune_test_list, function(x) retrieving_distribution_results(x))  


# Update t.test to something that desn't require normality 

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

