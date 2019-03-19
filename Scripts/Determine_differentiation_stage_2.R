library(devtools)
load_all("~/artdeco")
library(stringr)
library("MuSiC")
library("xbioc")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko/Scripts/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
meta_data = meta_info[colnames(expr_raw),]
#col_names = colnames(expr_raw)
#row_names = rownames(expr_raw)
#expr_raw = as.double(as.character(unlist(expr_raw)))
#expr_raw = matrix(expr_raw,ncol = length(col_names))
#colnames(expr_raw) = col_names
#rownames(expr_raw) = row_names
#meta_data = meta_data[which(meta_data$Location == "Primary"),] # GSE73338
#meta_data = meta_data[which(meta_data$Subtype == "Primary"),]# GSE98894
#meta_data = meta_data[which(meta_data$Location == "pancreas"),]# GSE98894
#meta_data = meta_data[which(meta_data$Grading != ""),] # GSE73338
#meta_data = meta_data[which(! is.na(meta_data$KI67) ),] # GSE73338
#meta_data = meta_data[which(meta_data$Location == "pancreas"),]
expr_raw = expr_raw[,rownames(meta_data)]
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]

#write.table(expr_raw,path_transcriptome_file,sep="\t",quote=F)
#write.table(expr,visualization_data_path,sep="\t",quote=F)

path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.41.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.24_samples.pancreatic_NEN.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Scarpa.29_samples.pancreatic_NEN.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.Pancreas.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.tsv"

visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.41.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.24_samples.pancreatic_NEN.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Scarpa.29_samples.pancreatic_NEN.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.ki67.Grading.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.Primary.Pancreas.vis.tsv"

visualization_data = read.table(visualization_data_path, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
transcriptome_data = read.table(path_transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
meta_data = meta_info[colnames(transcriptome_data),]
row_names = as.character(rownames(transcriptome_data))
col_names = as.character(colnames(transcriptome_data))
transcriptome_data = matrix(as.integer(as.character(unlist(transcriptome_data))),ncol = length(col_names))
rownames(transcriptome_data) = row_names
colnames(transcriptome_data) = col_names
colnames(visualization_data) = str_replace(colnames(visualization_data), pattern = "^X", "")
colnames(transcriptome_data) = str_replace(colnames(transcriptome_data), pattern = "^X", "")

selected_models = c(
    #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",

    #"Alpha_Beta_Gamma_Delta_Baron",

    #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe",

    #"Alpha_Beta_Gamma_Delta_Segerstolpe",

    #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor",

    "Alpha_Beta_Gamma_Delta_Lawlor",
    "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor"
)

deconvolution_results = Determine_differentiation_stage(
    transcriptome_data = transcriptome_data,
    deconvolution_algorithm = "music",
    models = selected_models,
    nr_permutations = 1000,
    output_file = ""
)

if (!(""%in% meta_info[rownames(deconvolution_results),"Grading"]))
    deconvolution_results[,"Grading"] = meta_info[rownames(deconvolution_results),"Grading"]

#deconvolution_results$Confidence_score_dif = log(deconvolution_results$Confidence_score_dif+1)

ki_index = which(rownames(transcriptome_file) == "MKI67")
if( length(ki_index) != 0 ){
    deconvolution_results[,"MKI67"] = rep(0,nrow(deconvolution_results))
    deconvolution_results[,"MKI67"] = log(as.double(transcriptome_file[ki_index[1],rownames(deconvolution_results)])+1)
    #deconvolution_results[,"MKI67"] = as.double(transcriptome_file[ki_index[1],])
} else {
    deconvolution_results[,"MKI67"] = as.double(meta_data[rownames(deconvolution_results),"KI67"])
}

deconvolution_results$Strength_de_differentiation[which(is.infinite(as.double(deconvolution_results$Strength_de_differentiation)))] = -4
deconvolution_results$Confidence_score_dif[which(is.infinite(as.double(deconvolution_results$Confidence_score_dif)))] = 0

vis_mat = create_heatmap_differentiation_stages(
    visualization_data,
    #transcriptome_file,
    deconvolution_results,
    #high_threshold = 10,
    confidence_threshold = .9,
    show_colnames = F,
    aggregate_differentiated_stages = F
)
visualization_data_path
selected_models

###

cor_t = cor.test((deconvolution_results[rownames(vis_mat),"MKI67"]),(vis_mat$Ratio_numeric));cor_t$p.value
cor_t = cor.test((deconvolution_results[,"MKI67"]),(deconvolution_results$ductal));cor_t$p.value
cor_t = cor.test((deconvolution_results[,"MKI67"]),(deconvolution_results$hisc));cor_t$p.value
chi_t = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$MKI67)),as.factor(as.character(vis_mat$Ratio))));chi_t$p.value
chi_t = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$MKI67)),as.factor(as.character(vis_mat$ductal))));chi_t$p.value
chi_t = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$MKI67)),as.factor(as.character(vis_mat$hisc))));chi_t$p.value
chi_t = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$Grading)),as.factor(as.character(vis_mat$Ratio))));chi_t$p.value
chi_t = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$Grading)),as.factor(as.character(vis_mat$ductal))));chi_t$p.value
chi_t = suppressWarnings(chisq.test(as.factor(as.character(vis_mat$Grading)),as.factor(as.character(vis_mat$hisc))));chi_t$p.value

anova_res = aov(as.double(vis_mat$Ratio_numeric) ~ as.factor(as.character(vis_mat$Grading)) );TukeyHSD(anova_res)
anova_res = aov(as.double(deconvolution_results[rownames(vis_mat),"ductal"]) ~ as.factor(as.character(vis_mat$Grading)) );TukeyHSD(anova_res)
anova_res = aov(as.double(deconvolution_results[rownames(vis_mat),"hisc"]) ~ as.factor(as.character(vis_mat$Grading)) );TukeyHSD(anova_res)

Grading_G1_G2 = vis_mat$Grading
Grading_G1_G2[Grading_G1_G2 %in% c("G1","G2")] = "G1_G2"
anova_res = aov(as.double(vis_mat$Ratio_numeric) ~ as.factor(as.character(Grading_G1_G2)) );TukeyHSD(anova_res)
aov(as.double(deconvolution_results$Ratio_numeric) ~ as.factor(as.character(vis_mat$Grading)) );TukeyHSD(anova_res)

#
car::leveneTest(as.double(vis_mat$Ratio_numeric)~grading)
aov_residuals <- residuals(object = anova_res )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )
stats::kruskal.test(as.double(vis_mat$Ratio_numeric)~grading)


### additional

mki67_mean = mean(mki67)
mki67_sd = sd(mki67)

anova_res = aov(as.double(vis_mat$Ratio_numeric) ~ as.factor(vis_mat$MKI67) )
summary(anova_res)

anova_res = aov(as.double(deconvolution_results$ductal) ~ as.factor(vis_mat$MKI67) )
summary(anova_res)

anova_res = aov(as.double(deconvolution_results$hisc) ~ as.factor(vis_mat$MKI67) )
summary(anova_res)
TukeyHSD(anova_res)

cor_tests = apply(
    expr_raw,
    MARGIN = 1,
    FUN = function(vec){
        test = cor.test(vis_mat$Ratio_numeric,as.double(vec))
        return(test$p.value)
    }
)
#cor_tests = as.double(cor_tests)
cor_tests[is.na(cor_tests)] = 1
hist(cor_tests)
cor_tests = sort(cor_tests, decreasing = F)
cor_tests

#write.table(cor_tests, "~/Deko/Results/Correlation_Baron_GSE98894.tsv",sep="\t",quote=F,row.names= T)
