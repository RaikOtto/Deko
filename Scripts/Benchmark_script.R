### benchmark runs
missing_samples = c("105103","130002","PNET08","130003","145702","1344","127403","PNET55")

library(devtools)
load_all("~/artdeco")
source("~/Deco/CIBERSORT_package/CIBERSORT.R")
library(stringr)
library("bseqsc")
library("MuSiC")
#library("xbioc")

models_ductal = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"))#,
    #list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe"))#,
    #list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor"))
)
models_hisc = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Hisc_Baron"))#,
    #list(c("Alpha_Beta_Gamma_Delta_Acinar_Hisc_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Segerstolpe"))#,
    #list(c("Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor"))
)
models_mixed = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron"))#,
    #list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Segerstolpe"))#,
    #list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Lawlor"))
)
nr_models = length(models_ductal)

transcriptome_files = list.files("~/Deco/Data/Bench_data/",full.names = T,pattern = "[0-9].tsv")
transcriptome_files = as.character(sapply(transcriptome_files,FUN=rep,3))
transcriptome_files = transcriptome_files[! (transcriptome_files %in% c("/home/ottoraik/Deco/Data/Bench_data//Wiedenmann.S39.tsv","/home/ottoraik/Deco/Data/Bench_data//Wiedenmann.S23.tsv"))]
visualization_files = str_replace_all(transcriptome_files,pattern ="\\.tsv",".vis.tsv")

meta_info = read.table("~/Deco//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deco//Scripts/Benchmark.R")

algorithm = "bseqsc" # NMF # music # bseqsc
type = "hisc"

high_threshold = 66
low_threshold = 33
confidence_threshold = 1.1

fractions <<- matrix( as.character(), ncol = 6)

for( i in 1:length(transcriptome_files)){

    dataset_query = tail(as.character(unlist(str_split(transcriptome_files[i],pattern = "/"))),1)
    dataset_query = str_replace_all(dataset_query,".tsv","")

    if (type == "ductal") {
        models = models_ductal[[1]]
    } else if  (type == "hisc") {
        models = models_hisc[[1]]
    }
    
    #dataset_training = as.character(unlist(models))[((i-1) %% 3) + 1]
    path_benchmark_files = paste0(
        "~/Deco/Results/Cell_fraction_predictions/",
        paste0(
            c(dataset_query,
              paste0(models[2], collapse = ".", sep =""),
              algorithm,"tsv"
            ),
            collapse = "."
        )
    )
    
    path_benchmark_files_dec_res = paste0(
        "~/Deco/Results/Cell_fraction_predictions/",
        paste0(
            c(dataset_query,
              paste0(models[2], collapse = ".", sep =""),
              algorithm,".dec_res.tsv"
            ),
            collapse = "."
        )
    )
    
    
    #if (! str_detect(models[2], pattern = "Baron") )
    #    next()
    
    transcriptome_file = transcriptome_files[i]
    visualization_file = visualization_files[i]
    
    print(i)
    print(dataset_query)
    print(models[2])
    
    if (
        file.exists(path_benchmark_files) |
        file.exists(str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS")))
    {
        #benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
        
        #if (nrow(benchmark_results_t) >= i)
            next(paste0("Skipping ",i))
    }
    
    res = run_benchmark(
        dataset_query = dataset_query,
        dataset_training = models,
        algorithm = algorithm,
        transcriptome_file = transcriptome_file,
        visualization_file = visualization_file,
        path_benchmark_files = path_benchmark_files,
        high_threshold = high_threshold,
        low_threshold = low_threshold,
        confidence_threshold = confidence_threshold,
        path_benchmark_files_dec_res = path_benchmark_files_dec_res
    )
    
    #write.table(res, path_benchmark_files, sep ="\t", quote = F, row.names = F)
    fractions = rbind(fractions,res)
}
