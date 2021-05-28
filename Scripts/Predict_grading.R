library("dplyr")
library(pROC)
library("caret")
library("stringr")
library("e1071")
set.seed(1)

# scRNA section

path_transcriptome_file = "~/Deko_Projekt/Results/All.S200.CIBERSORT.tsv"

cell_type_predictions = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = T, row.names = 1)
cell_type_predictions[1:5,1:5]

###

fitControl <- trainControl(
    method = "repeatedcv",
    number = 10,
    sampling = "down",
    savePred=T
)

pred_grading = function(
    train_mat,
    truth_vec,
    model = ""
){
    if ( model == ""){
        grading_model = caret::train(
            x = train_mat,
            y = truth_vec,
            method = "rf",
            norm.votes=T,
            #predict.all=FALSE,
            #preProcess = c("scale","center"),
            type = "Classification",
            metric= "Accuracy",
            ntree = 500,
            trControl = fitControl)
    }
    return(grading_model)
}

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S96.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
#colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
#expr_raw[1:5,1:5]

matcher = match(rownames(cell_type_predictions),meta_info$Sample,nomatch = 0)
rownames(cell_type_predictions)[matcher == 0]
rownames(cell_type_predictions)[matcher == 0] = paste("X",rownames(cell_type_predictions)[matcher == 0],sep ="")
matcher = match(rownames(cell_type_predictions),meta_info$Sample,nomatch = 0)
matcher == 0
meta_data = meta_info[rownames(cell_type_predictions),]
#

selector = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","RMSE","Correlation","P_value")
train_mat = cell_type_predictions[,selector]

randomizer_test = sample(1:nrow(train_mat),size = nrow(train_mat)*.80)
randomizer_hold_out = 1:200
randomizer_hold_out = randomizer_hold_out[-randomizer_test]

train_mat_test = train_mat[randomizer_test,]
train_mat_hold_out = train_mat[randomizer_hold_out,]

truth_vec = meta_info[rownames(train_mat),"Grading"]
truth_vec_test = meta_info[rownames(train_mat[randomizer_test,]),"Grading"]
truth_vec_hold_out = meta_info[rownames(train_mat[randomizer_hold_out,]),"Grading"]

grading_model = pred_grading(
    train_mat = train_mat,
    truth_vec = truth_vec
)

predictions = predict(grading_model, train_mat)
predictions_hold_out = predict(grading_model, train_mat_hold_out)
truth_hold_out = factor(truth_vec_hold_out, levels = c("G1","G2","G3"))

con_mat = confusionMatrix(
    predictions_hold_out,
    truth_hold_out,
    #factor(truth_vec,levels = c("G1","G2","G3")),
    positive = "G3"
)

con_mat$byClass

gbmImp <- varImp(grading_model, scale = FALSE)
plot(varImp(grading_model), top = 8)

#meta_info[rownames(train_mat),"Predicted_Grading"] = as.character(prediction_ml)

# discern NECs from NETs

# supervised RNA-seq section

path_transcriptome_file = "~/Deko_Projekt/Data/Bench_data/Riemer_Scarpa.S69.tsv"
path_transcriptome_file = "~/Deko_Projekt/Data/Bench_data/Scarpa.S29.tsv"
path_transcriptome_file = "~/Deko_Projekt/Data/Bench_data/Riemer.S40.tsv"
path_transcriptome_file = "~/Deko_Projekt/Data/Bench_data/Sadanandam.S29.tsv"
path_transcriptome_file = "~/Deko_Projekt/Data/Bench_data/Missaglia.S75.tsv"
train_mat = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(train_mat) = str_replace(colnames(train_mat), pattern = "^X", "")

row_var = apply( train_mat, FUN = var, MARGIN = 1)
train_mat = train_mat[ row_var > quantile(row_var)[2], ]
marker_genes = c( "MKI67")

nr_marker_genes = 100
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_alpha_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_beta_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_gamma_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_delta_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_acinar_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_ductal_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = c( marker_genes, read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/Dif_exp_hisc_400.tsv", stringsAsFactors = F,header = T)[1:nr_marker_genes,1] )
marker_genes = unique(marker_genes)
marker_genes = marker_genes[ marker_genes %in% rownames(train_mat)]
which(!(marker_genes %in% rownames(train_mat)))

train_mat = train_mat[marker_genes,]
"MKI67" %in% rownames(train_mat)

truth_vec = meta_info[colnames(train_mat),"Grading"]

res = pred_data(
    train_mat = t(train_mat),
    truth_vec = truth_vec
)
d = res$byClass

###

svmROC <- pROC::multiclass.roc(
    as.factor(predictions$predictions),
    as.factor(predictions$observations)
    #predValues$pred,
    #predValues$obs
)
gbmImp <- varImp(model_exp)
plot(varImp(model_exp, scale = FALSE), top = 20)

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

bench_data = read.table("~/Deko_Projekt/Results/Classification_Performance.tsv",sep="\t",stringsAsFactors = F,header = T)
bench_data = bench_data %>% filter(Dataset %in% c("RepSet"))
bench_data$Model[bench_data$Model == "Expression"] = "Ki-67"

# accuracy
# remotes::install_github("coolbutuseless/ggpattern")
#library("ggpattern")

p_accuracy_plot = ggplot(
    data = bench_data,
    aes(
        x = Grading,
        y = Balanced.Accuracy*100,
        fill = Model
) )+ theme(legend.position="none") + theme(legend.text=element_text(size=14))
p_accuracy_plot = p_accuracy_plot + geom_bar(stat="identity", position=position_dodge())
p_accuracy_plot = p_accuracy_plot + scale_fill_manual(values = c("Blue","red"))
p_accuracy_plot = p_accuracy_plot + xlab(label = "Accuracy") + ylab(label = "%")
p_accuracy_plot = p_accuracy_plot + theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14))

p_sensitivity_plot = ggplot(
    data = bench_data,
    aes(
        x = Grading,
        y = Sensitivity*100,
        fill = Model
    ) )+ theme(legend.position="none") + theme(legend.text=element_text(size=14))
p_sensitivity_plot = p_sensitivity_plot + geom_bar(stat="identity", position=position_dodge())
p_sensitivity_plot = p_sensitivity_plot + scale_fill_manual(values = c("Blue","red"))
p_sensitivity_plot = p_sensitivity_plot + xlab(label = "Sensitivity") + ylab(label = "%")
p_sensitivity_plot = p_sensitivity_plot + theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14))

p_specificity_plot = ggplot(
    data = bench_data,
    aes(
        x = Grading,
        y = Specificity * 100,
        fill = Model
    ) ) + theme(legend.position="none") + theme(legend.text=element_text(size=14))
p_specificity_plot = p_specificity_plot + geom_bar(stat="identity", position=position_dodge())
p_specificity_plot = p_specificity_plot + scale_fill_manual(values = c("Blue","red"))
p_specificity_plot = p_specificity_plot  + xlab(label = "Specificity") + ylab(label = "%")
p_specificity_plot = p_specificity_plot + theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14))

joint_plot = ggpubr::ggarrange(
    p_accuracy_plot,
    p_sensitivity_plot,
    p_specificity_plot,
    labels = c("", "", ""),
    ncol = 3,
    nrow = 1,
    common.legend = TRUE
    #legend.grob = get_legend(p_exo)
)

#svg(filename = "~/Deko_Projekt/Results/Images/Figure_4_G1_vs_G2_vs_G3_Grading_classification.svg", width = 10, height = 10)
joint_plot
dev.off()
