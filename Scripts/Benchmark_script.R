### benchmark runs

library(devtools)
load_all("~/artdeco")
library(stringr)
library("MuSiC")
library("xbioc")

path_transcriptome_files = c(
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.Primary.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.Pancreas.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.tsv",6)
)

path_visualization_files = c(
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.Primary.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.Pancreas.vis.tsv",6),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.vis.tsv",6)
)

models_exokrine = c(
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_progenitor_stanescu_hisc_haber")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Progenitor_Stanescu_HISC_Haber")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe_progenitor_stanescu_hisc_haber")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe","Progenitor_Stanescu_HISC_Haber")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor_progenitor_stanescu_hisc_haber")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor","Progenitor_Stanescu_HISC_Haber"))
)
models_endocrine = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Baron_progenitor_stanescu_hisc_haber")),
    list(c("Alpha_Beta_Gamma_Delta_Baron","Progenitor_Stanescu_HISC_Haber")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Segerstolpe_progenitor_stanescu_hisc_haber")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Progenitor_Stanescu_HISC_Haber")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Lawlor_progenitor_stanescu_hisc_haber")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Progenitor_Stanescu_HISC_Haber"))
)

benchmark_results_t_ori = data.frame(
    "Dataset_query" = as.character(),
    "Dataset_training" = as.character(),
    "P_val_cor_ratio_numerical_mki67" = as.double(),
    "P_val_cor_ductal_numerical_mki67" = as.double(),
    "P_val_cor_hisc_numerical_mki67" = as.double(),
    "P_val_cor_ratio_categorical_mki67" = as.double(),
    "P_val_cor_ductal_categorical_mki67" = as.double(),
    "P_val_cor_hisc_categorical_mki67" = as.double(),
    "P_val_chisq_ratio_categorical_grading_categorical" = as.double(),
    "P_val_chisq_ductal_categorical_grading_categorical" = as.double(),
    "P_val_chisq_hisc_categorical_grading_categorical" = as.double(),
    "Anova_G1_G2_Ratio" = as.double(),
    "Anova_G1_G3_Ratio" = as.double(),
    "Anova_G2_G3_Ratio" = as.double(),
    "Anova_G1_G2_Ductal" = as.double(),
    "Anova_G1_G3_Ductal" = as.double(),
    "Anova_G2_G3_Ductal" = as.double(),
    "Anova_G1_G2_hisc" = as.double(),
    "Anova_G1_G3_hisc" = as.double(),
    "Anova_G2_G3_hisc" = as.double()
)

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko/Scripts/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
meta_data = meta_info[colnames(expr_raw),]
#meta_data = meta_data[which(meta_data$Location == "Primary"),]
#meta_data = meta_data[which(meta_data$Grading != ""),]
#meta_data = meta_data[which(meta_data$KI67 > 0),]
#meta_data = meta_data[which(meta_data$Location == "pancreas"),]
expr_raw = expr_raw[,rownames(meta_data)]
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]

run_benchmark = function(
    dataset_query,
    dataset_training,
    algorithm,
    path_transcriptome_file,
    path_visualization_file,
    path_benchmark_files
){

    # prep
    
    transcriptome_data = read.table(path_transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
    colnames(transcriptome_data) = str_replace_all(colnames(transcriptome_data),pattern="^X","")
    meta_data = meta_info[colnames(transcriptome_data),]
    if( !nrow(meta_data) == ncol(transcriptome_data))
        stop("Inconclusive meta data and transcriptome file dimensions")
    row_names = as.character(rownames(transcriptome_data))
    col_names = as.character(colnames(transcriptome_data))
    transcriptome_data = matrix(as.integer(as.character(unlist(transcriptome_data))),ncol = length(col_names))
    rownames(transcriptome_data) = row_names
    colnames(transcriptome_data) = col_names
    
    visualization_data = read.table(path_visualization_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
    colnames(visualization_data) = str_replace_all(colnames(visualization_data),pattern="^X","")

    deconvolution_results = Determine_differentiation_stage(
        transcriptome_data = transcriptome_data,
        deconvolution_algorithm = str_to_lower(algorithm),
        models = dataset_training,
        nr_permutations = 1000,
        output_file = ""
    )
    
    if (!(""%in% meta_data[rownames(deconvolution_results),"Grading"]))
        deconvolution_results[,"Grading"] = meta_data[rownames(deconvolution_results),"Grading"]
    
    #deconvolution_results$Confidence_score_dif = log(deconvolution_results$Confidence_score_dif+1)
    
    ki_index = which(rownames(transcriptome_data) == "MKI67")
    if( length(ki_index) != 0 ){

        deconvolution_results[,"MKI67"] = rep(0,nrow(deconvolution_results))
        deconvolution_results[,"MKI67"] = log(as.double(transcriptome_data[ki_index[1],])+1)
        
    } else {
        
        deconvolution_results[,"MKI67"] = as.double(meta_data$KI67)
    }
    
    deconvolution_results$Strength_de_differentiation[which(is.infinite(as.double(deconvolution_results$Strength_de_differentiation)))] = -4
    deconvolution_results$Confidence_score_dif[which(is.infinite(as.double(deconvolution_results$Confidence_score_dif)))] = 0

    ### results parsing
    
    vis_mat = create_heatmap_differentiation_stages(
        visualization_data,
        deconvolution_results,
        #high_threshold = 10,
        confidence_threshold = .9,
        show_colnames = F,
        aggregate_differentiated_stages = F
    )
    
    cor_1 = cor.test((deconvolution_results[rownames(vis_mat),"MKI67"]),(vis_mat$Ratio_numeric))
    cor_1_p_value = cor_1$p.value
    if (length(deconvolution_results$ductal) > 0){
        cor_2 = cor.test((deconvolution_results[,"MKI67"]),(deconvolution_results$ductal))
        cor_2_p_value = cor_2$p.value
    } else {cor_2_p_value = 1.0}
    cor_3 = cor.test((deconvolution_results[,"MKI67"]),(deconvolution_results$hisc))
    cor_3_p_value = cor_3$p.value
    cor_4 = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$MKI67)),as.factor(as.character(vis_mat$Ratio))))
    cor_4_p_value = cor_4$p.value
    if (length(vis_mat$ductal) > 0){
        cor_5 = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$MKI67)),as.factor(as.character(vis_mat$ductal))))
        cor_5_p_value = cor_5$p.value
    } else {
        cor_5_p_value = 1.0
    }
    cor_6 = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$MKI67)),as.factor(as.character(vis_mat$hisc))))
    cor_6_p_value = cor_6$p.value
    
    if ( length(vis_mat$Grading) > 0){
    
            cor_7 = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$Grading)),as.factor(as.character(vis_mat$Ratio))))
            cor_7_p_value = cor_7$p.value
            
            if (length(vis_mat$ductal) > 0){
                cor_8 = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$Grading)),as.factor(as.character(vis_mat$ductal))))
                cor_8_p_value = cor_8$p.value
            } else {
                cor_8_p_value = 1.0
            }
            cor_9 = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$Grading)),as.factor(as.character(vis_mat$hisc))))
            cor_9_p_value = cor_9$p.value
            
            cor_10 = aov(as.double(vis_mat$Ratio_numeric) ~ as.factor(as.character(vis_mat$Grading)) )
            cor_10_p_value = as.double(TukeyHSD(cor_10)$`as.factor(as.character(vis_mat$Grading))`[,4])
            
            if (length(deconvolution_results[rownames(vis_mat),"ductal"]) > 0){
                cor_11 = aov(as.double(deconvolution_results[rownames(vis_mat),"ductal"]) ~ as.factor(as.character(vis_mat$Grading)) )
                cor_11_p_value = as.double(TukeyHSD(cor_11)$`as.factor(as.character(vis_mat$Grading))`[,4])
            } else {
                cor_11_p_value = c(1.0,1.0,1.0)
            }
            cor_12 = aov(as.double(deconvolution_results[rownames(vis_mat),"hisc"]) ~ as.factor(as.character(vis_mat$Grading)) )
            cor_12_p_value = as.double(TukeyHSD(cor_12)$`as.factor(as.character(vis_mat$Grading))`[,4])
    } else {
        cor_7_p_value = cor_8_p_value = cor_9_p_value = 1
        cor_10_p_value = cor_11_p_value = cor_12_p_value = rep(1,3)
    }
    
    identifier = tail(as.character(unlist(str_split(dataset_training[1],pattern = "_"))),1)
    if (str_detect(dataset_training[2], "HISC")){
        dataset_training_label = paste(identifier,"hisc",sep="_")
    } else {
        dataset_training_label = identifier
    }
    
    results_vec = c(
        dataset_query,
        dataset_training_label,
        cor_1_p_value,
        cor_2_p_value,
        cor_3_p_value,
        cor_4_p_value,
        cor_5_p_value,
        cor_6_p_value,
        cor_7_p_value,
        cor_8_p_value,
        cor_9_p_value,
        cor_10_p_value,
        cor_11_p_value,
        cor_12_p_value
    )
    
    benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
    benchmark_results_t = rbind(benchmark_results_t,results_vec)
    colnames(benchmark_results_t) = colnames(benchmark_results_t_ori)
    
    write.table(benchmark_results_t, path_benchmark_files,sep="\t",quote=F,row.names= F)
    
}

algorithm = "bseqsc"
type = "endocrine"
benchmark_results_t = benchmark_results_t_ori
path_benchmark_files = paste0(c("~/Deko/Results/Benchmark_results",algorithm,type,"tsv"),collapse = ".")

for( i in 1:length(path_transcriptome_files)){
    
    if (!file.exists(path_benchmark_files))
        write.table(benchmark_results_t,path_benchmark_files,sep="\t",quote=F,row.names= F)
    
    dataset_query = tail(as.character(unlist(str_split(path_transcriptome_files[i],pattern = "/"))),1)
    dataset_query = str_replace_all(dataset_query,".tsv","")
    
    if (type == "exocrine") {
        models = models_exokrine
    } else if  (type == "endocrine") {
        models = models_endocrine
    }
    
    dataset_training = as.character(unlist(models[((i-1) %% 6) + 1]))
    
    path_transcriptome_file = path_transcriptome_files[i]
    path_visualization_file = path_visualization_files[i]
    
    print(i)
    
    benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
    
    if (nrow(benchmark_results_t) >= i)
        next(paste0("Skipping ",i))

    run_benchmark(
        dataset_query = dataset_query,
        dataset_training = dataset_training,
        algorithm = algorithm,
        path_transcriptome_file = path_transcriptome_file,
        path_visualization_file = path_visualization_file,
        path_benchmark_files = path_benchmark_files
    )
}
