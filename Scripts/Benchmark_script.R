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

models_exocrine_like = c(
    "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
    "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe",
    "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor")
nr_models = length(models_exocrine_like)

#transcriptome_files = list.files("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_Bulk/",full.names = T,pattern = ".tsv")
#visualization_file =transcriptome_file = "~/Deko_Projekt/Data/Publication_datasets/Fröhling.S20.tsv"
#visualization_file = transcriptome_file = "~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN_NENs/Fröhling.S34.tsv"
visualization_file = transcriptome_file = "/home/ottoraik/Deko_Projekt/Data/Cancer_Pancreas_Bulk_Array/Sato.S35.tsv"
#visualization_file = "~/Deko_Projekt/Data/Publication_datasets/Fröhling.S20.DESeq2.tsv"
transcriptome_files = as.character(sapply(transcriptome_files,FUN=rep,3))
#visualization_files = str_replace_all(transcriptome_files,pattern ="\\.tsv",".vis.tsv")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko_Projekt/Scripts/Benchmark.R")

algorithm = "bseqsc" # NMF # music # bseqsc # SCDC
type = "Exocrine-like"

high_threshold = 66
low_threshold = 33
confidence_threshold = 1.1
transcriptome_files
i = 2
fractions <<- matrix( as.character(), ncol = 6)

#for( i in 1:length(transcriptome_files)){
#for( i in 16:18){

    dataset_query = tail(as.character(unlist(str_split(transcriptome_file,pattern = "/"))),1)
    dataset_query = str_replace_all(dataset_query,".tsv","")

    models = dataset_training = models_exocrine_like[i]

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

    print(i)
    print(dataset_query)
    print(dataset_training)
    "
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
    
    write.table(res, path_benchmark_files, sep ='\t', quote = F, row.names = F)
    "#fractions = rbind(fractions,res)
#}
