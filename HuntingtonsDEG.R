setwd("C:/R_projects/HDNEW.IISC")
getwd()


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install(c("limma", "Biobase"))
library(limma)
library(Biobase)

files <- list.files(
  pattern = "\\.txt$",
  full.names = TRUE
)

length(files)
head(readLines(files[1]), 40)
feat_line <- grep("^FEATURES", lines)
feat_line


dt <- fread(
  files[1],
  skip = feat_line + 1
)

dim(dt)
head(colnames(dt), 10)

lines <- readLines(files[1])

lines[feat_line:(feat_line + 5)]



#########################################################
#correct way to read aligned sample files
library(data.table)

read_agilent <- function(file) {
  
  lines <- readLines(file)
  feat_line <- grep("^FEATURES", lines)
  
  # extract header
  header <- strsplit(lines[feat_line], "\t")[[1]]
  
  # read data
  dt <- fread(
    file,
    skip = feat_line + 1,
    header = FALSE
  )
  
  # DROP first column ("DATA")
  dt <- dt[, -1]
  
  # assign column names (now lengths MATCH)
  colnames(dt) <- header[-1]
  
  # keep only biological probes
  dt <- dt[ControlType == 0]
  
  # return probe + expression
  dt[, .(ProbeName, gProcessedSignal)]
}

####################################################
#buid expression matrix 
expr_list <- lapply(files, read_agilent)

length(expr_list)       # should be 85
head(expr_list[[1]])


expr_list[[1]] #merege all samples in one list 
length(expr_list)
#now name each sample column properly
sample_names <- basename(files)

sample_names <- sub("\\.txt$", "", sample_names)
head(sample_names)
length(sample_names)

colnames(expr_list[[1]])


for(i in seq_along(expr_list)) {
  setnames(expr_list[[i]],
           old= "gProcessedSignal",
           new= sample_names[i])
}


head(expr_list[[1]])
#here we have duplicate probes so we will collapse first
#collapse duplicates inside each sample
expr_list <- lapply(expr_list, function(dt) {
  dt[, .(
    expression = mean(get(names(dt)[2]), na.rm = TRUE)
  ), by = ProbeName]
})
#rename the expression column properly 
for (i in seq_along(expr_list)) {
  setnames(expr_list[[i]],
           "expression",
           sample_names[i])
}
#merge
expr_merged <- Reduce(
  function(x, y) merge(x, y, by = "ProbeName", all = TRUE),
  expr_list
)

dim(expr_merged)
#sanity check
any(duplicated(expr_merged$ProbeName))
# should be FALSE

expr_mat <- as.matrix(expr_merged[, -1])
rownames(expr_mat) <- expr_merged$ProbeName

dim(expr_mat)

#Log2 transform
summary(expr_mat[,1])

expr_mat <- log2(expr_mat + 1)

summary(expr_mat[,1])
#we will use limma from here
library(limma)
expr_norm <- normalizeBetweenArrays(expr_mat, method = "quantile")

#check if normalization worked
boxplot(expr_norm[,1:10],
        outline = FALSE,
        main = "Post-normalization (first 10 samples)",
        las = 2)

dim(expr_norm)

sample_names <- colnames(expr_norm)

design <- model.matrix(~ group, data = samples)
colnames(design)
#metadata from scratch
#get sample names i.e columns of the matrix
sample_names <- colnames(expr_norm)
length(sample_names)
condition <- c( rep("Control", 42), 
                rep("Disease", 43) 
                )


meta <- data.frame(
  sample = sample_names,
  condition = factor(condition)
)

table(meta$condition)

#build a design matrix
design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)
design
#fit model for limma 
fit <- lmFit(expr_norm, design)
class(fit)
contrast.matrix <- makeContrasts(
  Disease_vs_Control = Disease - Control,
  levels = design
)
contrast.matrix
#apply contrast to limma model
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
#extract degs
deg <- topTable(
  fit2,
  coef = "Disease_vs_Control",
  number = Inf,
  adjust.method = "BH"
)

#sanity check imp
head(deg)
summary(deg$logFC)
table(deg$adj.P.Val < 0.05)

#count significant genes
#get stastical degs and strong effect degs
sum(deg$adj.P.Val < 0.05) #which genes are statistically confidently diff btwn disease and control
sum(abs(deg$logFC) > 1) #which gene change by at least  2-fold i.e strong effect degs


#this will give u strong effect genes only
deg_filt <- deg[
  abs(deg$logFC) >= 1 & deg$adj.P.Val < 0.05,
]
nrow(deg_filt)

up   <- deg_filt[deg_filt$logFC > 0, ]
down <- deg_filt[deg_filt$logFC < 0, ]


#primary set to be used for protein,ml
deg_strong <- deg[
  deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1,
]

#sec set for pathway,GSEA,KEGG,disease signature
deg_all_sig <- deg[deg$adj.P.Val < 0.05, ]



####################################################################################

#plot pca
graphics.off()
library(ggplot2)

pca <- prcomp(t(expr_norm), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  condition = meta$condition
)

ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: Expression profiles")


#######################################################################################
#heatmaps of all degs wrt to probe ids /unannoted
BiocManager::install("pheatmap")
library(pheatmap)
sample_dist <- dist(t(expr_norm))
sample_dist_mat <- as.matrix(sample_dist)
dim(sample_dist_mat)
rownames(sample_dist_mat) <- colnames(expr_norm)
colnames(sample_dist_mat) <- colnames(expr_norm)
#verify. must be TRUE
all(colnames(sample_dist_mat)==meta$sample)
#bcs i moved gpl files in the same folder as gsm samples, and it got stored in meta
#while re running these script so this is an extra step for that
meta <- meta[grepl("^GSM", rownames(meta)), ]
#fix rowm names in meta from numeric to ids
rownames(meta) <- meta$sample
head(rownames(meta))
annotation_col <- meta[, "condition", drop = FALSE]

ann_colors <- list(
  condition = c(
    Control = "#E76F51",
    Disease = "#2A9D8F"
  )
)



pheatmap(
  sample_dist_mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "complete"
)
###################################################################################

#volcano plots unannoted
library(ggplot2)

volcano_df <- deg
volcano_df$gene <- rownames(volcano_df)

volcano_df$significance <- "Not significant"
volcano_df$significance[
  volcano_df$adj.P.Val < 0.05 & volcano_df$logFC > 1
] <- "Upregulated"

volcano_df$significance[
  volcano_df$adj.P.Val < 0.05 & volcano_df$logFC < -1
] <- "Downregulated"

ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c(
    "Upregulated" = "#D62828",
    "Downregulated" = "#1D3557",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot: Disease vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal(base_size = 14)

###########################################################################

#safety check ...prepare deg table
head(rownames(deg))
#select top genes by adjusted p value
deg_strong <- deg[
  deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1,
]
#sort
deg_strong <- deg_strong[
  order(deg_strong$adj.P.Val),
]
#choose top 50 genes
top_n <- 50
top_genes <- rownames(deg_strong)[1:top_n]

#extract expression matrix
expr_top<- expr_norm[top_genes, ]
dim(expr_top)
# Z score scaling for standardisation..crucial
#it makes each gene to have mean of 0 and sd of 1 
expr_top_scaled<- t(scale(t(expr_top)))
#fix metadata
meta_heat<- meta
rownames(meta_heat)<- meta_heat$sample
meta_heat$sample<- NULL


#plot heat maps of top 50 gened unannoted just wrt to their probe ids 
library(pheatmap)

pheatmap(
  expr_top_scaled,
  annotation_col = meta_heat,
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  scale = "none",
  fontsize_row = 7,
  main = "Top 50 Differentially Expressed Genes"
)

#######################################################################

#for annotation in volcano plot
#load the GPL platform id for mouse 
library(dplyr)
#read annotation
annot_file <- "GPL13912_old_annotations.txt"
readLines(annot_file, n = 30)
annot <- read.delim(
  annot_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  skip = 19 # header line - 1
)
#check
colnames(annot)
head(annot$NAME)
head(annot$GENE_SYMBOL)
#clean annotation
annot_clean <- annot %>%
  dplyr::filter(CONTROL_TYPE == "FALSE") %>%
  dplyr::filter(GENE_SYMBOL != "" & !is.na(GENE_SYMBOL)) %>%
  dplyr::select(NAME, GENE_SYMBOL)
#check
head(annot_clean$GENE_SYMBOL)
head(annot_clean$NAME)
#verify
head(annot_clean)
#done with mapping probe  to gene names
# convert deg to dataframe
deg_df <- as.data.frame(deg)
deg_df$PROBE_ID <- rownames(deg_df)
# merge with annotation
deg_annot <- merge(
  deg_df,
  annot_clean,
  by.x = "PROBE_ID",
  by.y = "NAME",
  all.x = TRUE
)

# check
head(deg_annot[, c("PROBE_ID", "GENE_SYMBOL", "logFC", "adj.P.Val")])
colnames(annot_clean)
sum(is.na(deg_annot$GENE_SYMBOL))
nrow(deg_annot)
#plot

library(ggplot2)
library(ggrepel)
library(dplyr)
deg_annot_clean <- deg_annot %>%
  filter(!is.na(GENE_SYMBOL)) %>%
  mutate(negLogAdjP = -log10(adj.P.Val))

label_genes <- deg_annot_clean %>%
  filter(adj.P.Val < 0.01, abs(logFC) > 1.5) %>%
  arrange(adj.P.Val) %>%
  distinct(GENE_SYMBOL, .keep_all = TRUE) %>%
  head(15)

ggplot(deg_annot_clean, aes(x = logFC, y = negLogAdjP)) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1),
             alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = label_genes,
    aes(label = GENE_SYMBOL),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Disease vs Control Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Significant"
  )
#the above code will work but if it doesnt maually force it to open in png format
#by running the code below 


png("volcano_test.png", width = 1200, height = 1000, res = 150)

ggplot(deg_annot_clean, aes(x = logFC, y = negLogAdjP)) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1),
             alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = label_genes,
    aes(label = GENE_SYMBOL),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal(base_size = 14)

dev.off()
getwd()

###########################################################################
#annoted heatmaps of top  50 genes
# remove rows without gene symbols
deg_annot_clean <- deg_annot[!is.na(deg_annot$GENE_SYMBOL), ]
# Get top 50 genes by adjusted p-value
top50_genes <- deg_annot_clean %>%
  dplyr::distinct(GENE_SYMBOL, .keep_all = TRUE) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(50) %>%                     # ← base R head
  dplyr::pull(GENE_SYMBOL)

# extract gene names
top50_gene_names <- top50_genes$GENE_SYMBOL
# merge expression with annotation
expr_df <- as.data.frame(expr_norm)
expr_df$PROBE_ID <- rownames(expr_df)

expr_annot <- merge(
  expr_df,
  annot_clean,
  by.x = "PROBE_ID",
  by.y = "NAME"
)

# collapse probes → genes (mean expression)
library(dplyr)

expr_gene <- expr_annot %>%
  select(-PROBE_ID) %>%
  group_by(GENE_SYMBOL) %>%
  summarise(across(where(is.numeric), mean))

# convert to matrix
expr_gene <- as.data.frame(expr_gene)
rownames(expr_gene) <- expr_gene$GENE_SYMBOL
expr_gene$GENE_SYMBOL <- NULL

#subset expression matrix to top 50 geens
expr_top50<- expr_gene[rownames(expr_gene)%in% top50_gene_names, ]

#scale expression row
expr_top50_scaled <- t(scale(t(expr_top50)))
#build column annotation 
annotation_col <- data.frame(
  Condition = meta$condition
)
rownames(annotation_col) <- meta$sample
#now align it with matrix
annotation_col <- annotation_col[colnames(expr_top50_scaled), , drop = FALSE]

#final check
all(colnames(expr_top50_scaled) == rownames(annotation_col))
#plot annoted hetamap of top 50 genes
library(pheatmap)

pheatmap(
  expr_top50_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  scale = "none",
  main = "Top 50 Differentially Expressed Genes (Disease vs Control)"
)
ls()
head(expr_norm)
save.image("HD_project_workspace.RData")
load("HD_project_workspace.RData")
save(
  # expression
  expr_norm,
  expr_gene,
  
  # metadata
  meta,
  
  # DEG results
  deg,
  deg_annot,
  deg_annot_clean,
  deg_strong,
  
  # gene lists
  top50_genes,
  top50_gene_names,
  top_genes,
  
  # annotation
  annot_clean,
  
  file = "HD_core_ML_protein_ready.RData"
)
load("HD_core_ML_protein_ready.RData")

