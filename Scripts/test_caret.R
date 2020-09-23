library("dplyr")
library(pROC)
library("caret")
library("stringr")
library("e1071")
set.seed(1)

fitControl <- trainControl(
    method = "cv",
    number = 5,
    sampling = "down",
    savePred=T
)

pred_data = function(
    train_mat,
    truth_vec,
    model = ""
){
    if ( model == ""){
        model = caret::train(
            x = train_mat,
            y = truth_vec,
            method = "rf",
            norm.votes=T,
            #predict.all=FALSE,
            type = "Classification",
            metric= "Accuracy",
            ntree = 500,
            trControl = fitControl)
    }

    truth_vec = factor(truth_vec, levels = c("G1","G2","G3"))
    prediction_ml = predict(model, train_mat)
    con_mat = confusionMatrix(
        prediction_ml,
        truth_vec,
        positive = "G3"
    )

    return(con_mat)
}

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

# scRNA section

path_transcriptome_file = "~/Deko_Projekt/Results/Cell_fraction_predictions/Riemer_Scarpa.S69.Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.bseqsc..dec_res.tsv"
path_transcriptome_file = "~/Deko_Projekt/Results/Cell_fraction_predictions/Riemer.S40.Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.bseqsc..dec_res.tsv"
path_transcriptome_file = "~/Deko_Projekt/Results/Cell_fraction_predictions/Scarpa.S29.Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.bseqsc..dec_res.tsv"
path_transcriptome_file = "~/Deko_Projekt/Results/Cell_fraction_predictions/Sadanandam.S29.Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.bseqsc..dec_res.tsv"
path_transcriptome_file = "~/Deko_Projekt/Results/Cell_fraction_predictions/Missaglia.S75.Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.bseqsc..dec_res.tsv"

cell_type_predictions = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = T)

colnames(cell_type_predictions) = str_replace(colnames(cell_type_predictions), pattern = "^X", "")
colnames(cell_type_predictions) = str_replace_all(colnames(cell_type_predictions), pattern = "\\.", "_")
colnames(cell_type_predictions) = str_replace_all(colnames(cell_type_predictions), pattern = "-", "_")
colnames(cell_type_predictions) = str_to_lower(colnames(cell_type_predictions))
rownames(cell_type_predictions) = 1:nrow(cell_type_predictions)
cell_type_predictions = cell_type_predictions %>% filter(model %in% "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")
rownames(cell_type_predictions) = str_replace(cell_type_predictions$sample, pattern = "^X", "")

remove(train_mat)
train_mat = cell_type_predictions[,
   match(c("alpha","beta","gamma","delta","acinar","ductal","p_value","correlation","rmse"), colnames(cell_type_predictions), nomatch = 0)
]
truth_vec = meta_info[rownames(train_mat),"Grading"]

res = pred_data(
    train_mat = train_mat,
    truth_vec = truth_vec
)
d = res$byClass

gbmImp <- varImp(model_deco, scale = FALSE)
plot(varImp(model_deco), top = 20)

# discern NECs from NETs

meta_data = meta_info[rownames(train_mat),]
train_mat = t(train_mat)
truth_vec = as.character(meta_data$NEC_NET)
train_mat[1:5,1:5]

model = caret::train(
    x = train_mat,
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

bench_data = read.table("~/Deko_Projekt/Results/Classification_Performance.tsv",sep="\t",stringsAsFactors = T,header = T)
bench_data = bench_data %>% filter(Dataset %in% c("RepSet"))

# accuracy
# remotes::install_github("coolbutuseless/ggpattern")
library("ggpattern")

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

svg(filename = "~/Deko_Projekt/Results/Images/Figure_4_G1_vs_G2_vs_G3_Grading_classification.svg", width = 10, height = 10)
joint_plot
dev.off()
