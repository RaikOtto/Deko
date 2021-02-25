### benchmark runs
# missing_samples = c("105103","130002","PNET08","130003","145702","1344","127403","PNET55")

library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")
library("MuSiC")
library("SCDC")

###
#transcriptome_data = read.table("~/Deko_Projekt/Data/Bench_data/Scarpa.S29.tsv",sep = "\t",header = T,row.names = 1)
#colnames(transcriptome_data) = str_replace(colnames(transcriptome_data) , pattern ="^X","")
#transcriptome_data[1:5,1:5]

#models_ductal = list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Hisc_Baron"))

models_ductal = c(
    list(c("Alpha_Beta_Gamma_Delta_Baron","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")),
    list(c("Alpha_Beta_Gamma_Delta_Segerstolpe","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor")),
    list(c("Alpha_Beta_Gamma_Delta_Lawlor","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_Segerstolpe"))
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

#meta_info = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko_Projekt/Scripts/Benchmark.R")

algorithm = "bseqsc" # NMF # music # bseqsc # SCDC
type = "ductal"
#type = "HISC"

high_threshold = 66
low_threshold = 33
confidence_threshold = 1.1
transcriptome_files
i = 4
fractions <<- matrix( as.character(), ncol = 6)

#for( i in 1:length(transcriptome_files)){
#for( i in 16:18){

    dataset_query = tail(as.character(unlist(str_split(transcriptome_files[i],pattern = "/"))),1)
    dataset_query = str_replace_all(dataset_query,".tsv","")

    if (type == "ductal") {
        models = models_ductal#[[1]]
    } else if  (type == "hisc") {
        models = models_hisc#[[1]]
    }
    
    models = models[ (i - 1) %% 3 + 1]
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
    
    print(i)
    print(dataset_query)
    print(dataset_training)
    
    res = run_benchmark(
        dataset_query = dataset_query,
        dataset_training = dataset_training,
        type = type,
        algorithm = algorithm,
        transcriptome_file = transcriptome_file,
        visualization_file = visualization_file,
        path_benchmark_files = path_benchmark_files,
        high_threshold = high_threshold,
        low_threshold = low_threshold,
        confidence_threshold = confidence_threshold,
        path_benchmark_files_dec_res = path_benchmark_files_dec_res
    )
    
    write.table(res, path_benchmark_files, sep ="\t", quote = F, row.names = F)
    #fractions = rbind(fractions,res)
#}
