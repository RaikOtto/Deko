###

deconvolution_results$NEC_NET = meta_data[rownames(deconvolution_results),"NEC_NET"]
deconvolution_results_selection = subset(deconvolution_results, Grading %in% c("G3"))
vis_mat = deconvolution_results_selection[,c("NEC_NET","ductal","MKI67","hisc")]
vis_mat = reshape2::melt(vis_mat)
colnames(vis_mat) = c("NEC_NET","Type","Value")

p = ggplot( data = vis_mat,aes( x = NEC_NET, y = Value, fill =Type) )
p = p + geom_bar(stat="identity", position=position_dodge())
p

anova_mki67 = aov(deconvolution_results$MKI67 ~ as.factor(as.character(deconvolution_results$NEC_NET)) )
TukeyHSD(anova_mki67)$`as.factor(as.character(deconvolution_results$NEC_NET))`

anova_ductal = aov(deconvolution_results$ductal ~ as.factor(as.character(deconvolution_results$NEC_NET)) )
TukeyHSD(anova_ductal)$`as.factor(as.character(deconvolution_results$NEC_NET))`

anova_hisc = aov(deconvolution_results$hisc ~ as.factor(as.character(deconvolution_results$NEC_NET)) )
TukeyHSD(anova_hisc)$`as.factor(as.character(deconvolution_results$NEC_NET))`

###

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


### ROCR
library("ROCR")
library("InformationValue")

perf_vec <<- c()
subtypes = c( "HISC", "ductal", "mki67")
res_mat <<-  matrix(as.character(),ncol = 6)

for( i in seq(7,24,by=3)){
    
    dataset_query = tail(as.character(unlist(str_split(transcriptome_files[i],pattern = "/"))),1)
    dataset_query = str_replace_all(dataset_query,".tsv","")
    
    if (type == "ductal") {
        models = models_ductal
    } else if  (type == "hisc") {
        models = models_hisc
    }
    
    dataset_training = as.character(unlist(models[((i-1) %% 3) + 1]))
    
    transcriptome_file = transcriptome_files[i]
    visualization_file = visualization_files[i]
    
    print(i)
    
    model_path = paste0(
        c(
            "~/Deko/Data/Bench_data/Models/",
            dataset_query,
            dataset_training[2],
            algorithm,
            "RDS"
        ),
        collapse = "."
    )
    model_path = str_replace_all(model_path, pattern = "/\\.","/")
    deconvolution_results = readRDS(model_path)
    
    meta_data = meta_info[rownames(deconvolution_results),]
    meta_data = meta_data[meta_data$Grading!="",]
    deconvolution_results = deconvolution_results[rownames(meta_data),]
    deconvolution_results$Grading = meta_data$Grading
    
    target_vector = deconvolution_results$Grading
    #inclusion_vector = which(target_vector %in% c("G1","G2","G3"))
    #target_vector = target_vector[inclusion_vector]
    target_vector[target_vector %in% c("G1","G2")] = 0
    target_vector[target_vector != 0] = 1
    
    HISC = deconvolution_results$hisc#[inclusion_vector]
    HISC = HISC / 100
    ductal = deconvolution_results$ductal#[inclusion_vector]
    ductal = ductal / 100
    
    transcriptome_data = read.table(transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
    ki_index = which(rownames(transcriptome_data) == "MKI67")
    colnames(transcriptome_data) = str_replace_all(colnames(transcriptome_data),pattern="^X","")
    deconvolution_results$MKI67 = rep("",nrow(deconvolution_results))
    deconvolution_results$MKI67 = as.double(transcriptome_data[ki_index,rownames(deconvolution_results)])
    
    mki67 = deconvolution_results$MKI67#[inclusion_vector]
    mki67 = mki67 / 100
    
    ###
    
    print ( dataset_query)

    for (subtype in subtypes ){
        
        prediction_vector = eval(
            parse(
                text = subtype
            )
        )
        
        t_data = data.frame(
            Grading = as.factor(target_vector),
            Value =  as.double( prediction_vector  )
        )
        
        rf_fit <- glm(
            Grading ~ Value, data = t_data, family=binomial(link="logit")
            #method = "glm"#,
            #lambda = 0
        )
        predicted <- plogis(predict(rf_fit, t_data))  # predicted scores
        optCutOff <- optimalCutoff(t_data$Grading, predicted)[1] 
        sensitivity = round(InformationValue::sensitivity(actuals = as.double(target_vector),predicted, threshold = optCutOff),2)
        specificity = round(InformationValue::specificity(actuals = as.double(target_vector),predicted, threshold = optCutOff),2)

        F1_score =  round(2*(sensitivity*specificity/(sensitivity+specificity)),2)
        
        pred_obj = ROCR::prediction(
            predictions = as.double( prediction_vector  ),
            labels = target_vector
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

        print( paste0( c( subtype,"Sensitivity:", sensitivity,", Specificity:", specificity, ", F1:", F1_score, ", ROC:",rocr_auc) ), collapse = " " )
        perf = ROCR::performance(
            prediction.obj = pred_obj,
            measure = "tpr",
            x.measure = "fpr"
        );
        
        perf_vec <<- c(perf_vec,perf)
        res_vec = c(dataset_query, subtype,sensitivity,specificity,F1_score,rocr_auc)
        res_mat =  rbind(res_mat,res_vec)
    }
}
colnames(res_mat) = c("Dataset","Proportion","Sensitivity","Specificity","F1","AUC")
write.table(res_mat,"~/Deco/Results/Figure_4_Classification_G1_&_G2_versus_G3_Performance.tsv",sep="\t",row.names = F)
# graphics

perf_1 = unlist(perf_vec)[[1]]
perf_2 = unlist(perf_vec)[[2]]
perf_3 = unlist(perf_vec)[[3]]
perf_4 = unlist(perf_vec)[[4]]

x_val_1 = ( as.double(unlist(perf_1@x.values ) ) ) * 100
x_val_2 = ( as.double(unlist(perf_2@x.values ) ) ) * 100
x_val_3 = ( as.double(unlist(perf_3@x.values ) ) ) * 100
x_val_4 = ( as.double(unlist(perf_4@x.values ) ) ) * 100

y_val_1 = ( as.double( unlist(perf_1@y.values) ) ) * 100
y_val_2 = ( as.double( unlist(perf_2@y.values) ) ) * 100
y_val_3 = ( as.double( unlist(perf_3@y.values) ) ) * 100
y_val_4 = ( as.double( unlist(perf_4@y.values) ) ) * 100

sub_classifier = as.factor( 
    c( 
        rep(subtypes[1],length(x_val_1)),
        rep(subtypes[2],length(x_val_2)),
        rep(subtypes[3],length(x_val_3)),
        rep(subtypes[4],length(x_val_4))
    )
)

plots_frame = as.data.frame( 
    cbind( 
        c(
            x_val_1,
            x_val_2, 
            x_val_3,
            x_val_4
        ),
        c(
            y_val_1,
            y_val_2,
            y_val_3,
            y_val_4
        )
    )
)
plots_frame$sub_classifier = sub_classifier
colnames(plots_frame) = colnames = c("X","Y","Subtype_classifier")

q_bird = ggplot(
    data = plots_frame, 
    aes( x = X, y = Y )
) + geom_line( size = 2, aes( linetype = Subtype_classifier, color = Subtype_classifier) )
q_bird  = q_bird + xlab( "False Positive Rate" ) + ylab( "Sensitivity" ) + geom_abline( intercept = 0, slope =1)
q_bird = q_bird + theme( 
    panel.background = element_blank(),
    legend.position  = "top",
    panel.border = element_rect(colour = "black", fill=NA, size = 2 )
)
q_bird
