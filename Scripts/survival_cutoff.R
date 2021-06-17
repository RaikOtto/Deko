algorithm = "bseqsc" # NMF # music # bseqsc
type = "hisc"
selector_var ="ductal"
#i = 4

### benchmark runs
# missing_samples = c("105103","130002","PNET08","130003","145702","1344","127403","PNET55")
if (F){

library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library(stringr)
library("bseqsc")
library("MuSiC")

}

models_ductal = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor"))
)
models_hisc = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Segerstolpe")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor"))
)
nr_models = length(models_ductal)

transcriptome_files = list.files("~/Deko_Projekt/Data/Bench_data/",full.names = T,pattern = "[0-9].tsv")
transcriptome_files = as.character(sapply(transcriptome_files,FUN=rep,3))
visualization_files = str_replace_all(transcriptome_files,pattern ="\\.tsv",".vis.tsv")

meta_info = read.table("~/MAPTor_NET///Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko_Projekt/Scripts/Benchmark.R")
fractions <<- matrix( as.character(), ncol = 6)

dataset_query = tail(as.character(unlist(str_split(transcriptome_files[i],pattern = "/"))),1)
dataset_query = str_replace_all(dataset_query,".tsv","")

if (type == "ductal") {
    models = models_ductal#[[1]]
} else if  (type == "hisc") {
    models = models_hisc#[[1]]
}

models = models[((i-1) %% 3) + 1]
dataset_training = as.character(unlist(models))[2]

path_benchmark_files = paste0(
    "~/Deko_Projekt/Results/Cell_fraction_predictions/",
    paste0(
        c(dataset_query,
          dataset_training,
          #paste0(models[2], collapse = ".", sep =""),
          algorithm,"tsv"
        ),
        collapse = "."
    )
)

path_benchmark_files_dec_res = paste0(
    "~/Deko_Projekt/Results/Cell_fraction_predictions/",
    paste0(
        c(dataset_query,
          dataset_training,
          #paste0(models[2], collapse = ".", sep =""),
          algorithm,".dec_res.tsv"
        ),
        collapse = "."
    )
)

transcriptome_file = transcriptome_files[i]
visualization_file = visualization_files[i]
algorithm = str_to_lower(algorithm)
transcriptome_data = read.table(transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
colnames(transcriptome_data) = str_replace_all(colnames(transcriptome_data),pattern="^X","")
meta_data = meta_info[colnames(transcriptome_data),]
if(  !nrow(meta_data) == ncol(transcriptome_data))
    stop("Inconclusive meta data and transcriptome file dimensions")
row_names = as.character(rownames(transcriptome_data))
col_names = as.character(colnames(transcriptome_data))
transcriptome_data = matrix(as.double(as.character(unlist(transcriptome_data))),ncol = length(col_names))
rownames(transcriptome_data) = row_names
colnames(transcriptome_data) = col_names

visualization_data = read.table(visualization_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
colnames(visualization_data) = str_replace_all(colnames(visualization_data),pattern="^X","")

name_query_data = tail(as.character(unlist(str_split(transcriptome_file,pattern ="/"))),1)
name_query_data = str_replace_all(name_query_data,pattern =".tsv","")

name_training_data = tail(as.character(unlist(str_split(dataset_training,pattern ="_"))),1)

name_algorithm_datatype = tail(as.character(unlist(str_split(path_benchmark_files,pattern ="/"))),1)
name_algorithm_datatype = str_replace_all(name_algorithm_datatype,pattern ="(.tsv)|(Benchmark_results.)","")
name_algorithm = head(as.character(unlist(str_split(name_algorithm_datatype,pattern ="\\."))),1)
name_datatype = tail(as.character(unlist(str_split(name_algorithm_datatype,pattern ="\\."))),1)

identifier = tail( as.character(unlist(str_split(dataset_training, pattern = "_"))), 1)
ki_index = which(rownames(transcriptome_data) == "MKI67")

###

model_path = paste0(
    c(
        "~/Deko_Projekt/Data/Bench_data/Models/",
        dataset_query,
        type,
        identifier,
        algorithm,
        "RDS"
    ),
    collapse = "."
)
model_path = str_replace_all(model_path, pattern = "/\\.","/")

###
if (
    file.exists(str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS"))
){
    #benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
    deconvolution_results = readRDS(
        str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS"))
    
} else {
    
    deconvolution_results = Deconvolve_transcriptome(
        transcriptome_data = transcriptome_data,
        deconvolution_algorithm = algorithm,
        models = dataset_training,
        nr_permutations = 1000,
        output_file = ""
    )
    
    if ( type == "ductal")
        deconvolution_results = deconvolution_results[
            grep(deconvolution_results$model, pattern = "Alpha_Beta_Gamma_Delta_Acinar_Ductal", value = F) ,]
    if ( type == "hisc")
        deconvolution_results = deconvolution_results[
            grep(deconvolution_results$model, pattern = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc", value = F) ,]
    
    if (algorithm == "bseqsc")
        deconvolution_results$P_value = as.double(deconvolution_results$P_value)
    
    col_names = colnames(deconvolution_results)
    width = ncol(deconvolution_results)
    height = nrow(deconvolution_results) 
    res_table = matrix(unlist(deconvolution_results),ncol = width, nrow = height)
    colnames(res_table) = col_names
    res_table= as.data.frame(res_table)
    
    saveRDS(
        deconvolution_results,
        str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS")
    )
    
    write.table(
        deconvolution_results,
        path_benchmark_files_dec_res,
        sep="\t",
        quote =F, row.names = F
    )
}

#return(deconvolution_results)

meta_data = meta_info[ rownames(deconvolution_results),]

if (sum( meta_data[rownames(deconvolution_results),"Grading"] != "") > 0){
    
    meta_data = meta_data[meta_data$Grading!="",]
    deconvolution_results = deconvolution_results[rownames(meta_data),]
    deconvolution_results$Grading = meta_data$Grading
    visualization_data  = visualization_data[,rownames(meta_data)]
    transcriptome_data  = transcriptome_data[,rownames(meta_data)]
    
}
#write.table(deconvolution_results,"~/Deco/Results/Cell_fraction_predictions/RepSet_Cibersort_Baron.tsv",sep ="\t", row.names =T , quote=F)

if( length(ki_index) != 0 ){
    
    deconvolution_results[,"MKI67"] = rep(0,nrow(deconvolution_results))
    deconvolution_results[,"MKI67"] = log(as.double(transcriptome_data[ki_index[1],])+1)
    
} else {
    
    deconvolution_results[,"MKI67"] = as.double(meta_data$KI67)
}

deconvolution_results$Strength_de_differentiation = 0
deconvolution_results$Confidence_score_dif = 0
deconvolution_results$Subtype = ""

vis_mat = create_visualization_matrix(
    visualization_data = visualization_data,
    deconvolution_results = deconvolution_results,
    confidence_threshold = 0.05,
    high_threshold = 66,
    low_threshold = 33
)

### survival curve

target_genes = expr_raw[sad_genes[sad_genes %in% rownames(expr_raw)] ,]
target_genes = as.double(apply(target_genes, FUN = mean, MARGIN = 2))

#vis_mat = vis_mat[rownames(deconvolution_results),]
vis_mat = vis_mat[rownames(meta_data),]
vis_mat$OS_Tissue = as.double(str_replace_all(meta_data$OS_Tissue, pattern = ",", "\\."))
vis_mat$OS_Tissue[is.na(vis_mat$OS_Tissue)] = 1
vis_mat$Grading = meta_data$Grading
vis_mat$Zensur = meta_data$Zensur

ratio_m = data.frame(
    #"alpha" = as.double(deconvolution_results[,"alpha"]),
    #"ductal" = as.double(deconvolution_results[,"ductal"]),
    #"hisc" = as.double(deconvolution_results[,"hisc"]),
    "grading" = meta_data$Grading,
    "OS_Tissue" = vis_mat$OS_Tissue,
    "Zensur" = vis_mat$Zensur,
    "MKi67"  = meta_data[,"MKI67"],
    "target_genes" = target_genes 
)
#if ("hisc" %in% colnames(deconvolution_results))
#    ratio_m$hisc = as.double(deconvolution_results[,"hisc"])

ratio_m = ratio_m[!is.na(ratio_m$Zensur),]

#selector_var = "ductal"
#selector_var = "MKi67"
#selector_var = "alpha"
selector_var = "target_genes"
#ratio_m$grading[ratio_m$grading %in% c("G1","G2")] = "G1_G2"

#agg = aggregate(ratio_m[,selector_var], FUN = mean, by = list(ratio_m$grading))
#thresh_low = (agg[1,2] + agg[2,2]) / 2
#thresh_mid = (tail(agg[,2],1) + tail(agg[,2],2)[1]) /2
thresh_mean = mean(target_genes)#agg$x[1]

value = ratio_m[,selector_var]
#thresh_mean = mean(value)
classification = rep("high",length(value))
classification[value <= thresh_mean] = "low_mid"
#classification[value <= thresh_mid] = "mid"
#classification[value <= thresh_low] = "low"
#classification = factor(classification, levels = c("low","mid","high"))
ratio_m[,selector_var] = classification

fit = survival::survfit( survival::Surv( as.double(ratio_m$OS_Tissue), ratio_m$Zensur ) ~ ratio_m[,selector_var], data = ratio_m)
surv_p_value = survminer::surv_pvalue(fit, data = ratio_m)$pval

print(i)
selector_var

surv_p_value

#pdf(graphics_path_survival_hisc,onefile = FALSE)#,width="1024px",height="768px")
#svg(filename = "~/Deko_Projekt/Results/Images/Figure_5_survival_grading_three.svg", width = 10, height = 10)
#svg(filename = "~/Deko_Projekt/Results/Images/Figure_5_survival_ductal_three.svg", width = 10, height = 10)
    print(survminer::ggsurvplot(fit, data = ratio_m, risk.table = F, pval = T, censor.size = 10))
#dev.off()
