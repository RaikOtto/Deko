library(pROC)
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

# supervised RNA-seq section

path_transcriptome_file = "~/Deco/Data/Bench_data/MAPTor_NET.S57.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

row_var = apply( expr_raw, FUN = var, MARGIN = 1)
summary(row_var)
expr_raw = expr_raw[ row_var > quantile(row_var)[2], ]
dim(expr_raw)

descrCor <- cor(t(expr_raw))
summary(descrCor[upper.tri(descrCor)])

highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
expr_raw <- expr_raw[-highlyCorDescr,] ###
dim(expr_raw)
descrCor2 <- cor(t(expr_raw))
summary(descrCor2[upper.tri(descrCor2)])

expr_raw_save = expr_raw
#comboInfo <- findLinearCombos(t(expr_raw))
#filteredDescr = t(expr_raw)[-comboInfo$remove,]
#dim(expr_raw)

marker_genes = c()
nr_marker_genes = 100
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

meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)
expr_raw = t(expr_raw)

###

table(meta_data$Grading)

truth_vec = meta_data$Grading
truth_vec[truth_vec %in% c("G1","G2")] = "G1_G2"
truth_vec[truth_vec %in% c("G3")] = "G3"
dim(expr_raw)
length(truth_vec)

bootControl = trainControl(
    method = "LOOCV",
    classProbs=T,
    savePredictions = T#,
    #number = 10
    #summaryFunction = twoClassSummary
)

model = caret::train(
    x = expr_raw,
    y = truth_vec,
    method = "multinom",
    trControl = bootControl,
    scaled = FALSE,
    tuneLength = 1,# model$finalModel$decay 0.0005623413
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

predictions = data.frame(
    "predictions" = model$pred$pred,
    "observations" = model$pred$obs
)
confusionMatrix(
    as.factor(predictions$predictions),
    as.factor(predictions$observations),
    positive = "G3"
)
svmROC <- pROC::roc(
    predValues$pred,
    predValues$obs
)

# Select a parameter setting
selectedIndices <- predValues$pred$mtry == 2
# Plot:
plot.roc(
    predValues$pred$obs,
    predValues$pred$M
)

gbmImp <- varImp(model, scale = FALSE)
plot(varImp(model), top = 20)

indices = list(
    "1" = seq(1,1+5),
    "2" = seq(7,7+5),
    "3" = seq(13,13+5),
    "4" = seq(19,19+5),
    "5" = seq(25,25+5),
    "6" = seq(31,31+5),
    "7" = seq(37,37+4),
    "8" = seq(42,42+4),
    "9" = seq(47,47+4),
    "10" = seq(52,57)
)

sensitivity <- c()
specificity <- c()
accuracy <- c()
PPV <- c()
#Kappa <<- c()

for (i in 1:10){
    predictions = as.character(model$pred$pred)[
        as.integer(unlist(indices[i]))
    ]
    observations = as.character(model$pred$obs)[
        as.integer(unlist(indices[i]))
    ] 
    predictions = data.frame(
        "predictions" = predictions,
        "observations" = observations
    )
    conf_mat = confusionMatrix(
        as.factor(predictions$predictions),
        as.factor(predictions$observations),
        positive = "G3"
    )
    sensitivity = c(sensitivity, as.double(conf_mat$byClass["Sensitivity"]))
    specificity = c(specificity, as.double(conf_mat$byClass["Specificity"]))
    accuracy = c(accuracy, as.double((conf_mat$byClass["Balanced Accuracy"])))
    PPV = c(PPV, as.double((conf_mat$byClass["Pos Pred Value"])))
    Kappa = c(Kappa, as.double((conf_mat$byClass["Pos Pred Value"])))
}

sd(as.double(as.character(unlist(sensitivity))))
sd(as.double(as.character(unlist(specificity))))
sd(as.double(as.character(unlist(accuracy))))
sd(as.double(as.character(unlist(PPV)))[-9])

