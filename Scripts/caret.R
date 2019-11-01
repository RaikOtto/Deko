library(caret)
library(ranger)
library(tidyverse)
library(e1071)
library(readr) # Useful for importing data
library(knitr) # Useful for creating nice tables
library(data.table)
library(tm) #Text Manager

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

# load in abalone dataset
  <- read.table("~/Downloads/abalone.tsv", sep = ",")
# load in column names
colnames(abalone_data) <- c("sex", "length", "diameter", "height", 
                            "whole.weight", "shucked.weight", 
                            "viscera.weight", "shell.weight", "age")
# add a logical variable for "old" (age > 10)
abalone_data <- abalone_data %>%
    mutate(old = age > 10) %>%
    # remove the "age" variable
    select(-age)

# split into training and testing
set.seed(23489)
train_index <- sample(1:nrow(abalone_data), 0.9 * nrow(abalone_data))
abalone_train <- abalone_data[train_index, ]
abalone_test <- abalone_data[-train_index, ]
# remove the original dataset
rm(abalone_data)
# view the first 6 rows of the training data
head(abalone_train)

dim(abalone_train)

# fit a random forest model (using ranger)
rf_fit <- train(as.factor(old) ~ ., 
                data = abalone_train, 
                method = "ranger")

rf_fit


###
library("stringr")
library("grid")

path_transcriptome_file = "~/Deco/Data/Bench_data/MAPTor_NET.S57.tsv"

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
    return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/Deco//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_data = meta_info[colnames(expr_raw),]

truth_vec = meta_data$Grading
truth_vec[truth_vec %in% c("G1","G2")] = 0
truth_vec[truth_vec %in% c("G3")] = 1

train_mat = t(expr_raw)[,1:100]

t_data = data.frame(
    cbind(
        Grading = as.double(truth_vec),
        train_mat
    )
)

rf_fit <- glm(
    Grading ~ ., data = t_data, family=binomial(link="logit")
    #method = "glm"#,
    #lambda = 0
)
prediction_vector <- plogis(predict(rf_fit, as.data.frame(train_mat)))  # predicted scores
optCutOff = InformationValue::optimalCutoff(as.double(t_data$Grading), predicted)[1] 
sensitivity = round(InformationValue::sensitivity(actuals = as.double(truth_vec),predicted, threshold = optCutOff),2)
specificity = round(InformationValue::specificity(actuals = as.double(truth_vec),predicted, threshold = optCutOff),2)

F1_score =  round(2*(sensitivity*specificity/(sensitivity+specificity)),2)

pred_obj = ROCR::prediction(
    predictions = as.double( prediction_vector  ),
    labels = truth_vec
)

rocr_auc = as.character(
    round(
        as.double(
            unclass(
                ROCR::performance(
                    pred_obj,
                    "auc"
                )
            )@"y.values"
        ),
        2 )
)

print( paste0( c( "Sensitivity:", sensitivity,", Specificity:", specificity, ", F1:", F1_score, ", ROC:",rocr_auc) ), collapse = " " )
perf = ROCR::performance(
    prediction.obj = pred_obj,
    measure = "tpr",
    x.measure = "fpr"
);

perf_vec <<- c(perf_vec,perf)
res_vec = c(dataset_query, subtype,sensitivity,specificity,F1_score,rocr_auc)
res_mat =  rbind(res_mat,res_vec)