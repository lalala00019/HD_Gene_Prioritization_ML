#ML-based gene prioritisation
#huntington's disease
setwd("C:/R_projects/HDNEW.IISC/outputs")
getwd()
list.files()
#well use deg_strong matrix that has probe level deg and then later convert to gene level
head(rownames(deg_strong)) 
#extract probe ids from rows 
deg_probes <- rownames(deg_strong)
length(deg_probes)
#subset expression matrix
expr_deg <- expr_norm[rownames(expr_norm) %in% deg_probes, ]
dim(expr_deg)
#transpose for ml 
ml_data <- t(expr_deg)
dim(ml_data)
#add class label
ml_data <- as.data.frame(ml_data)
ml_data$Condition <- as.factor(meta$condition)#error cuz meta has 83 samples and ml_data has 85
table(ml_data$Condition)
rownames(ml_data)[1:10]  # first 10 sample names in expression
head(meta$sample)     # assuming meta has sample IDs
#remove non-sample rows from ml_data, remove the 2 gpl ids 
# keep only rows that are in meta$sample
ml_data <- ml_data[rownames(ml_data) %in% meta$sample, ]
dim(ml_data)
#align meta to expr rows
meta <- meta[match(rownames(ml_data), meta$sample), ]
length(meta$sample)  # should be 83
#add condition column
ml_data$Condition <- as.factor(meta$condition)
table(ml_data$Condition)
#confirm
nrow(ml_data)
ncol(ml_data) # 237 + 1 (Condition)
head(ml_data$Condition)


#train random forest on ml_data
install.packages("randomForest")
library(randomForest)
#train Rf
rf_model<- randomForest(Condition~ ., data=ml_data,
                        ntree=500, important=TRUE)
#prediction
pred_class<- predict(rf_model)
table(pred_class, ml_data$Condition) #confusion matrix on training data
#feature importance
importance_mat <- importance(rf_model)
head(importance_mat)
#build a clean importance table
importance_df <- data.frame(
  Probe = rownames(importance_mat),
  MeanDecreaseGini = importance_mat[, "MeanDecreaseGini"]
)

importance_df <- importance_df[
  order(importance_df$MeanDecreaseGini, decreasing = TRUE),
]

head(importance_df, 10)
#take top probes
top_probes <- importance_df$Probe[1:20]
top_probes
#map probes to gene symbols
top_genes <- annot_clean$GENE_SYMBOL[
  match(top_probes, annot_clean$NAME)
]

top_genes
#clean the NA's...final gene prioritization
top_genes <- unique(na.omit(top_genes))
top_genes
length(top_genes)

#ROC curve
#get class probabilities
install.packages("pROC")
library(pROC)
#get class probabilities
pred_prob <- predict(rf_model, type = "prob")[, "Disease"]
#compute ROC
roc_obj <- roc(ml_data$Condition, pred_prob)
#plot
plot(roc_obj, main = "ROC Curve – Random Forest Classifier")
auc(roc_obj)

#create a ranked table 
# Map probes to gene symbols in the importance table
importance_df$GENE_SYMBOL <- annot_clean$GENE_SYMBOL[
  match(importance_df$Probe, annot_clean$NAME)
]

# Keep only rows with valid gene symbols
importance_df_clean <- importance_df[!is.na(importance_df$GENE_SYMBOL), ]

# Take top 14
top14_df <- importance_df_clean[1:14, ]

top14_df


write.csv(top14_df,
          file = "Top14_ML_prioritized_genes_RandomForest.csv",
          row.names = FALSE)

#here we trained and tested on the same data which is only for gene prioritization
#now we'll try train-test split
set.seed(123)
train_idx<- sample(
  1:nrow(ml_data),
  size=0.7*nrow(ml_data)
)
train_data<- ml_data[train_idx, ]
test_data<- ml_data[-train_idx, ]
#train RF only on training data
rf_model <- randomForest(
  Condition ~ .,
  data = train_data,
  ntree = 500,
  importance = TRUE
)
#predict on test data 
test_pred <- predict(rf_model, test_data)
table(test_pred, test_data$Condition)
#permutation sanity check, test by shuffling the labels 
set.seed(123)
shuffled <- test_data
shuffled$Condition <- sample(shuffled$Condition)
shuf_pred <- predict(rf_model, shuffled)
table(shuf_pred, shuffled$Condition)
#OOB error test 
rf_model$err.rate[rf_model$ntree, ]
