### benchmark runs

library(devtools)
load_all("~/artdeco")
library(stringr)
library("MuSiC")
library("xbioc")


models_ductal = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor"))
)
models_hisc = c(
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Segerstolpe")),
    list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor"))
)
models_mixed = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Segerstolpe")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor"))
)
nr_models = length(models_ductal)

path_transcriptome_files = c(
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.Scarpa.tsv",nr_models), # 69
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.tsv",nr_models), #40
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Scarpa.tsv",nr_models), # 29
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/Sadanandam.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/Alvarez.tsv",nr_models)
)

path_visualization_files = c(
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.Scarpa.vis.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Wiedenmann.vis.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Scarpa.vis.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73338/Sadanandam.vis.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE73339/GSE73339.vis.tsv",nr_models),
    rep("~/Deko/Data/Cancer_Pancreas_Bulk_Array/GSE98894/Alvarez.vis.tsv",nr_models)
)

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko/Scripts/Benchmark.R")

algorithm = "nmf"
type = "hisc"
path_benchmark_files = paste0(c("~/Deko/Results/Benchmark_results",algorithm,type,"tsv"),collapse = ".")
high_threshold = 66
low_threshold = 33
confidence_threshold = 1.1

for( i in 1:length(path_transcriptome_files)){
    
    dataset_query = tail(as.character(unlist(str_split(path_transcriptome_files[i],pattern = "/"))),1)
    dataset_query = str_replace_all(dataset_query,".tsv","")
    
    if (type == "ductal") {
        models = models_ductal
    } else if  (type == "hisc") {
        models = models_hisc
    } else if  (type == "mixed") {
        models = models_mixed
    }
    
    dataset_training = as.character(unlist(models[((i-1) %% 3) + 1]))
    
    path_transcriptome_file = path_transcriptome_files[i]
    path_visualization_file = path_visualization_files[i]
    
    print(i)
    
    if (file.exists(path_benchmark_files)){
        benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
        
        if (nrow(benchmark_results_t) >= i)
            next(paste0("Skipping ",i))
    }
    
    run_benchmark(
        dataset_query = dataset_query,
        dataset_training = dataset_training,
        algorithm = algorithm,
        path_transcriptome_file = path_transcriptome_file,
        path_visualization_file = path_visualization_file,
        path_benchmark_files = path_benchmark_files,
        high_threshold = high_threshold,
        low_threshold = low_threshold,
        confidence_threshold = confidence_threshold
    )
}
