library("stringr")
library("grid")
library("caret")
library("ranger")
library("tidyverse")
library("e1071")

opt.cut = function( perf, pred ){
    cut.ind = mapply( FUN = function( x, y, p ){
        d = (x - 0)^2 + ( y - 1 )^2
        ind = which(d == min(d))
        c(
            sensitivity = y[[ ind ]],
            specificity = 1 - x[[ ind ]], 
            cutoff = p[[ ind ]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}

# split into training and testing
set.seed(23489)

# fit a random forest model (using ranger)

res_mat = read.table("~/Dropbox/testproject/Results/Prediction_results_baseline.tsv", header = TRUE)

truth_vec = res_mat$Grading_Binary
truth_vec[truth_vec == "G3"] = 0
truth_vec[truth_vec == "G1_G2"] = 1
truth_vec = as.double(truth_vec)
table(truth_vec)
prediction_vec = res_mat$Prediction_Grading_binary
prediction_vec[prediction_vec == "G3"] = 0
prediction_vec[prediction_vec == "G1_G2"] = 1
prediction_vec = as.double(prediction_vec)
table(prediction_vec)

pred_obj = ROCR::prediction(
    labels = truth_vec,
    predictions = prediction_vec
)
ROCR::performance(pred_obj,"sens")

rocr_auc = round(as.double(unclass(ROCR::performance(pred_obj,"auc"))@"y.values"),2 )

print( paste0( c( "Sensitivity:", sensitivity,", Specificity:", specificity, ", F1:", F1_score, ", ROC:",rocr_auc) ), collapse = " " )
perf_vec = ROCR::performance(
    prediction.obj = pred_obj,
    measure = "tpr",
    x.measure = "fpr"
);

perf_vec <<- c(perf_vec,perf)
res_vec = c(dataset_query, subtype,sensitivity,specificity,F1_score,rocr_auc)
res_mat =  rbind(res_mat,res_vec)

predictor_mat = cbind(train_mat,meta_data$Grading)
#write.table(predictor_mat,"~/Downloads/ml_data.tsv",quot =F, sep ="\t")
