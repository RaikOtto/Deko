library("caret")
library("stringr")
set.seed(1)

meta_info = read.table("~/Deco//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

# scRNA section

path_transcriptome_file = "~/Deco/Results/Cell_fraction_predictions/RepSet_Cibersort_Baron.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace_all(colnames(expr_raw), pattern = "\\.", "_")
colnames(expr_raw) = str_replace_all(colnames(expr_raw), pattern = "-", "_")
colnames(expr_raw) = make.names(colnames(expr_raw))
rownames(expr_raw) = str_replace(rownames(expr_raw), pattern = "^X", "")

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
rownames(expr_raw) = str_replace(rownames(expr_raw), pattern = "^X", "")
meta_data = meta_info[rownames(expr_raw),]

expr_raw = expr_raw[rownames(meta_data),]
rownames(expr_raw) = make.names(rownames(expr_raw))

expr_raw = expr_raw[, colnames(expr_raw) != "Subtype"  ]
col_var = apply( expr_raw, FUN = var, MARGIN = 2)
expr_raw = expr_raw[, col_var > 0.1]

# supervised section

path_transcriptome_file = "~/Deco/Data/Bench_data/MAPTor_NET.S57.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

marker_genes = c()
nr_marker_genes = 50
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_alpha_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_beta_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_gamma_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_delta_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_acinar_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_ductal_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deco/Results/Cell_fraction_predictions/Dif_exp_hisc_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = unique(marker_genes)
marker_genes = marker_genes[ marker_genes %in% rownames(expr_raw)]
table(marker_genes %in% rownames(expr_raw))

expr_raw = expr_raw[marker_genes,]
dim(expr_raw)

row_var = apply( expr_raw, FUN = var, MARGIN = 1)
summary(row_var)
expr_raw = expr_raw[ row_var > 10, ]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)
expr_raw = t(expr_raw)

###

table(meta_data$Grading)

truth_vec = meta_data$Grading
#truth_vec[truth_vec %in% c("G1","G2")] = "G1_G2"
#truth_vec[truth_vec %in% c("G3")] = "G3"
dim(expr_raw)
length(truth_vec)

bootControl = trainControl(
    method = "repeatedcv",
    classProbs=T,
    savePredictions = T,
    number = 10
    #summaryFunction = twoClassSummary
)

model = caret::train(
    x = expr_raw,
    y = truth_vec,
    method = "multinom",
    trControl = bootControl,
    scaled = FALSE,
    tuneLength = 10,
    preProcess = c("center", "scale")
)
summary(model)
model$finalModel
plot(model)
plot(model, metric = "Kappa")
resampleHist(model)

predValues = extractPrediction(
    models = list(model)#,
    #testX = t(expr_raw[good_genes,]),
    #testY = truth_vec#,
#    type = "prob"
)

plotClassProbs(predValues$model)

preds = predict(
    model,
    newdata = expr_raw
)
confusionMatrix(
    #as.factor(model$pred$pred),
    #as.factor(model$pred$obs)
    as.factor(preds),
    reference = as.factor(truth_vec)
)
svmROC <- pROC::roc(
    predValues$pred,
    predValues$obs
)

library(pROC)
# Select a parameter setting
selectedIndices <- predValues$pred$mtry == 2
# Plot:
plot.roc(
    predValues$pred$obs,
    predValues$pred$M
)

gbmImp <- varImp(model, scale = FALSE)
plot(varImp(model), top = 20)
