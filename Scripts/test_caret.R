library("dplyr")
library(pROC)
library("caret")
library("stringr")
library("e1071")
set.seed(1)

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

# scRNA section

path_transcriptome_file = "~/Deko_Projekt/Results/Bseq_results_fractions_p_values.tsv"
cell_type_predictions = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F)
colnames(cell_type_predictions) = str_replace(colnames(cell_type_predictions), pattern = "^X", "")
colnames(cell_type_predictions) = str_replace_all(colnames(cell_type_predictions), pattern = "\\.", "_")
colnames(cell_type_predictions) = str_replace_all(colnames(cell_type_predictions), pattern = "-", "_")
colnames(cell_type_predictions) = make.names(colnames(cell_type_predictions))
rownames(cell_type_predictions) = str_replace(rownames(cell_type_predictions), pattern = "^X", "")

cell_type_predictions = cell_type_predictions %>% filter(model %in% "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")
#cell_type_predictions = cell_type_predictions %>% filter(model %in% "Alpha_Beta_Gamma_Delta_HISC_Baron")
expr_raw = cell_type_predictions %>% filter( (Dataset %in% "RepSet"))
expr_raw = expr_raw %>% filter(grading %in% c("G1","G2","G3"))
rownames(expr_raw) = expr_raw$sample_id
expr_raw = expr_raw[,c("alpha","beta","gamma","delta","acinar","ductal","p_value")]

rownames(expr_raw) = str_replace( rownames(expr_raw), pattern = "^X", "")
dim(expr_raw)

# discern NECs from NETs

expr_raw = t(expr_raw)
meta_data = meta_info[rownames(expr_raw),]
truth_vec = as.character(meta_data$NEC_NET)
expr_raw[1:5,1:5]

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
    metric = "Accuracy",
    trControl = bootControl,
    scaled = FALSE,
    tuneLength = 5,# model$finalModel$decay 0.0005623413
    preProcess = c("center", "scale")
)

prediction_ml = model$pred$pred
observations_ml = model$pred$obs
predictions = data.frame(
    "predictions" = prediction_ml,
    "observations" = observations_ml
)
m = confusionMatrix(
    as.factor(predictions$predictions),
    as.factor(predictions$observations),
    positive = "NEC"
)

# supervised RNA-seq section

path_transcriptome_file = "~/Deko_Projekt/Data/Bench_data/Riemer_Scarpa.S69.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
rownames(expr_raw)

row_var = apply( expr_raw, FUN = var, MARGIN = 1)
summary(row_var)
expr_raw = expr_raw[ row_var > quantile(row_var)[2], ]
dim(expr_raw)

expr_raw_save = expr_raw

marker_genes = c()
nr_marker_genes = 100
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_alpha_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_beta_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_gamma_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_delta_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_acinar_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_ductal_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_hisc_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = unique(marker_genes)
marker_genes = marker_genes[ marker_genes %in% rownames(expr_raw)]
table(marker_genes %in% rownames(expr_raw))

expr_raw = expr_raw[marker_genes,]
dim(expr_raw)

meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)
expr_raw = t(expr_raw)

###

truth_vec = meta_info[rownames(expr_raw),"Grading"]
truth_vec[truth_vec %in% c("G1","G2")] = "G1_G2"
truth_vec[truth_vec %in% c("G3")] = "G3"
dim(expr_raw)

length(truth_vec)
#truth_vec = factor(truth_vec, levels = c("G3","G2","G1"))

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
    metric = "Accuracy",
    trControl = bootControl,
    scaled = FALSE,
    tuneLength = 5,# model$finalModel$decay 0.0005623413
    preProcess = c("center", "scale")
)
summary(model)
#model$finalModel
#plot(model, metric = "Accuracy")
#plot(model, metric = "Kappa")
#resampleHist(model)

#predValues = extractPrediction(
#    models = list(model)#,
    #testX = t(expr_raw[good_genes,]),
    #testY = truth_vec#,
#    type = "prob"
#)

#plotClassProbs(predValues$obs)


prediction_ml = model$pred$pred
observations_ml = model$pred$obs
#prediction_data = cell_type_predictions %>% filter(cell_type_predictions$Dataset %in% "RepSet") # Missiaglia RepSet Sadanandam Scarpa Wiedenmann
#observations_ml = as.character(meta_info[as.character(prediction_data$sample_id),"Grading"])
#observations_ml[observations_ml %in% c("G1","G2")] = "G1_G2"
#prediction_ml = as.character(predict(model,newdata = prediction_data[,c("alpha","beta","gamma","delta","acinar","ductal","p_value")]))

predictions = data.frame(
    "predictions" = prediction_ml,
    "observations" = observations_ml
)
m = confusionMatrix(
    as.factor(predictions$predictions),
    as.factor(predictions$observations),
    positive = "G3"
)
colMeans(m$byClass)
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
    #Kappa = c(Kappa, as.double((conf_mat$byClass["Pos Pred Value"])))
}

sd(as.double(as.character(unlist(sensitivity))))
sd(as.double(as.character(unlist(specificity))))
sd(as.double(as.character(unlist(accuracy))))
sd(as.double(as.character(unlist(PPV)))[-9])

### visualization

bench_data = read.table("~/Deko_Projekt/Results/Classification_Performance.tsv",sep="\t",stringsAsFactors = T,header = T)
bench_data$Dataset = factor(bench_data$Dataset,levels = c("Missiaglia","Scarpa","RepSet","Riemer","Sadanandam","Average"))
bench_data$Type = factor(bench_data$Type, levels = c("Unsupervised","Supervised"))

p = ggplot(
    data = bench_data,
    aes(
        x = Dataset,
        y = Accuracy,
        fill = Type
) )
p = p + geom_bar(stat="identity", position=position_dodge())
p = p + scale_fill_manual(values = c("Blue","red"))
p = p + theme(legend.position="top") + xlab(label = "NEN data set")
p = p + guides(fill=guide_legend(title="Classification model")) 

p
