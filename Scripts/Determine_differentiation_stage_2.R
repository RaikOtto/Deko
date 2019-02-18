library(devtools)
load_all("~/artdeco")
library(stringr)
library(MuSiC)
library(xbioc)

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
#grading = meta_info[colnames(expr_raw),"Grading"]
#expr_raw = expr_raw[,!is.na(grading)]
meta_data = meta_info[colnames(expr_raw),]
meta_data = meta_data[which(meta_data$Location == "Primary"),]
#meta_data = meta_data[which(meta_data$Subtype == "Primary"),]
#meta_data = meta_data[which(meta_data$Location == "pancreas"),]
expr_raw = expr_raw[,rownames(meta_data)]
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]

#write.table(expr,visualization_data_path,sep="\t",quote=F)
#write.table(expr_raw,path_transcriptome_file,sep="\t",quote=F)

#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.tsv"
path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.primary.pancreas.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.tsv"
#path_transcriptome_file = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.tsv"

#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.vis.tsv"
visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/GSE98894.primary.pancreas.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/GSE73338.vis.tsv"
#visualization_data_path = "~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.vis.tsv"

visualization_data = read.table(visualization_data_path, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
transcriptome_file = read.table(path_transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
colnames(visualization_data) = str_replace(colnames(visualization_data), pattern = "^X", "")
colnames(transcriptome_file) = str_replace(colnames(transcriptome_file), pattern = "^X", "")

deconvolution_results = Determine_differentiation_stage(
    transcriptome_file = transcriptome_file,
    deconvolution_algorithm = "bseqsc",
    models = c(
        "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
        "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_progenitor_stanescu_hisc_haber"
        #"Alpha_Beta_Gamma_Delta_Baron",
        #"Alpha_Beta_Gamma_Delta_Baron_progenitor_stanescu_hisc_haber"
        #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe",
        #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe_progenitor_stanescu_hisc_haber"
        #"Alpha_Beta_Gamma_Delta_Segerstolpe",
        #"Alpha_Beta_Gamma_Delta_Segerstolpe_progenitor_stanescu_hisc_haber"
        #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor",
        #"Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor_progenitor_stanescu_hisc_haber"
        #"Alpha_Beta_Gamma_Delta_Lawlor",
        #"Alpha_Beta_Gamma_Delta_Lawlor_progenitor_stanescu_hisc_haber"
        #"Progenitor_Stanescu_HISC_Haber"
    ),
    nr_permutations = 1000,
    output_file = ""
)

if (!(""%in% meta_info[rownames(deconvolution_results),"Grading"]))
    deconvolution_results[,"Grading"] = meta_info[rownames(deconvolution_results),"Grading"]

#deconvolution_results$Confidence_score_dif = log(deconvolution_results$Confidence_score_dif+1)

ki_index = grep(rownames(transcriptome_file),pattern = "MKI67")
if( length(ki_index) != 0 ){
    deconvolution_results[,"MKI67"] = rep(0,nrow(deconvolution_results))
    deconvolution_results[,"MKI67"] = log(as.double(transcriptome_file[ki_index[1],]))
}

vis_mat = create_heatmap_differentiation_stages(
    visualization_data,
    #transcriptome_file,
    deconvolution_results,
    #high_threshold = 10,
    confidence_threshold = .9,
    show_colnames = F,
    aggregate_differentiated_stages = F
)

# 1

mki67 = deconvolution_results[,"MKI67"]
cor.test(mki67,vis_mat$Ratio_numeric)
cor.test(mki67,deconvolution_results$ductal)
cor.test(mki67,deconvolution_results$hisc)

# 2

mki67_cat = as.factor(as.character(vis_mat$MKI67))
chisq.test(mki67_cat,as.factor(as.character(vis_mat$Ratio)))
chisq.test(mki67_cat,as.factor(as.character(vis_mat$ductal)))
chisq.test(mki67_cat,as.factor(as.character(vis_mat$hisc)))

# 3

grading = as.factor(as.character(vis_mat$Grading))
ratio = as.factor(as.character((vis_mat$Ratio)))
hisc = as.factor(as.character(vis_mat$hisc))
ductal = as.factor(as.character(vis_mat$ductal))

chisq.test(grading,ratio)
chisq.test(grading,ductal)
chisq.test(grading,hisc)

# 4

anova_res = aov(as.double(vis_mat$Ratio_numeric) ~ grading )
summary(anova_res)

anova_res = aov(as.double(deconvolution_results$ductal) ~ grading )
summary(anova_res)

anova_res = aov(as.double(deconvolution_results$hisc) ~ grading )
summary(anova_res)
TukeyHSD(anova_res)

#
car::leveneTest(as.double(vis_mat$Ratio_numeric)~grading)
aov_residuals <- residuals(object = anova_res )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )
stats::kruskal.test(as.double(vis_mat$Ratio_numeric)~grading)
