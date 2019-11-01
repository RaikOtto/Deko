library("caret")
library("stringr")
set.seed(1)

path_transcriptome_file = "~/Deco/Data/Bench_data/MAPTor_NET.S57.tsv"
path_transcriptome_file = "~/Deco/Data/Bench_data/"

meta_info = read.table("~/Deco//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
meta_data = meta_info[colnames(expr_raw),]
table(meta_data$Grading)
expr_raw = expr_raw[,rownames(meta_data)]

truth_vec = meta_data$Grading
truth_vec[truth_vec %in% c("G1","G2")] = 0
truth_vec[truth_vec %in% c("G3")] = 1

row_var = apply( expr_raw, MARGIN = 1, FUN = var)
percentiles_lowest_10 = quantile(col_var, probs = 1:100/100)[1]
dim(expr_raw)
expr_raw = expr_raw[row_var >= percentiles_lowest_10,]
dim(expr_raw)

bootControl = trainControl(
    method = "cv",
    #number = 200,
    classProbs=T,
    savePredictions = T
)
dim(expr_raw)
length(truth_vec)

model = caret::train(
    t(expr_raw),
    truth_vec,
    method = "multinom",
    trControl = bootControl,
    scaled = FALSE
)

model$finalModel
plot(model)
plot(model, metric = "Kappa")
resampleHist(model)

predict(model, newdata = t(expr_raw))[1:5]

predValues = extractPrediction(
    models = list(model),
    testX = t(expr_raw),
    testY = truth_vec
)

plotClassProbs(testProbs)

svmPred <- subset(
    predValues,
    model == "multinom"
)

#svmProb <- subset(testProbs, model == "svmRadial")

confusionMatrix(predValues$pred, predValues$obs)
svmROC <- roc(predValues$pred, predValues$obs)

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
