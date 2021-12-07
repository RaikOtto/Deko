devtools::install_github("ERGlass/lrcde.dev", build_vignettes=TRUE)
library(lrcde) # Load the lrcde package

path_visualization_files = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.vis.tsv"

scRNA_file_path_query = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.tsv"
scRNA_file_path_train = "~/Deko/Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.tsv"
model_name = str_replace_all(scRNA_file_path,pattern = "\\.tsv","")
model_name = tail(str_split(model_name,pattern = "/")[[1]],1)

t = read.table(scRNA_file_path_train, sep ="\t", header = T, row.names = 1, nrows = 1)
colnames(t) = str_replace_all(colnames(t),pattern ="\\.","_")
colnames(t) = str_replace_all(colnames(t),pattern ="^X","")
dim(t)
meta_data = meta_info[colnames(t),]
subtype_vector = str_to_lower(meta_data$Subtype)
table(subtype_vector)
transcriptome_data = read.table(scRNA_file_path,sep="\t",header  = T)



###

timestamp <- Sys.time()
library(caret)
library(plyr)
library(recipes)
library(dplyr)

model <- "nnls"

#########################################################################

set.seed(251)
ctrl <- trainControl(method = "repeatedcv", repeats = 10)

rctrl1 <- trainControl(method = "cv", number = 10, returnResamp = "all")

test_reg_cv_model <- train(
    trainX,
    trainY,
    method = "nnls",
    trControl = rctrl1
)


test_reg_pred <- predict(test_reg_cv_model, testX)
test_reg_cv_form <- train(y ~ ., data = training, method = "nnls", trControl = rctrl1)
test_reg_pred_form <- predict(test_reg_cv_form, testX)

test_reg_loo_model <- train(trainX, trainY, method = "nnls", trControl = rctrl2)

Marker_Gene_List <<- list()
for( subtype in unique(subtype_vector) ){
    Marker_Gene_List[[subtype]] = identify_marker_genes(
        expression_training_mat = expression_training_mat,
        subtype_vector = subtype_vector,
        subtype = subtype,
        nr_marker_genes = training_nr_marker_genes
    )
}

# Prepare bseq training
training_mat_bseq = new(
    "ExpressionSet",
    exprs = as.matrix(expression_training_mat)
)

test_reg_none_model <- train(
    transcriptome_data,
    #subtype_vector, 
    method = "nnls", 
    #trControl = rctrl3,
    tuneLength = 1,
    preProc = c("center", "scale")
)
test_reg_none_pred <- predict(test_reg_none_model, testX)

test_reg_predictors1 <- predictors(test_reg_cv_model)

test_reg_imp <- varImp(test_reg_cv_model)
