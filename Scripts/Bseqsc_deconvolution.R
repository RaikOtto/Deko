library("devtools")
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")
library("stringr")
library("reshape2")
library("bseqsc")
library("dplyr")

remove(props)
datasets = c("Alvarez.S105.tsv","Charite.S23.tsv","Diedisheim.S62.tsv","Master.S20.tsv","Missiaglia.S75.tsv","Sadanandam.S29.tsv","Sato.S18.NEC_NET.tsv","Scarpa.S29.tsv")

dataset_name = datasets[8]
i_filename = "~/Deko_Projekt/Data/Publication_datasets/NEN/Sato.S35.tsv"
i_filename = paste(i_filename, dataset_name, sep ="")
expr_raw = read.table(i_filename,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

#show_models_bseqsc()
#model_name = "Alpha_Beta_Gamma_Delta_Baron"
model_name = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"

props = Deconvolve_transcriptome(
    transcriptome_data = expr_raw,
    deconvolution_algorithm = "bseqsc",
    models = model_name,
    Cibersort_absolute_mode = FALSE,
    nr_permutations = 1000,
    output_file = ""
)

o_filename = "~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron_exocrine/NEN/"
o_filename = paste(o_filename, dataset_name, sep ="/")
write.table(props,o_filename,sep = "\t")
