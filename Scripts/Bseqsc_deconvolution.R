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

dataset_name = datasets[2]
i_filename = "~/Deko_Projekt/Data/Publication_datasets/"
i_filename = paste(i_filename, dataset_name, sep ="")
expr_raw = read.table(i_filename,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

meta_data = meta_info[colnames(expr_raw),]
table(meta_data$Histology_Primary)

#show_models_bseqsc()
model_name = "Alpha_Beta_Gamma_Delta_Baron"
#model_name = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"

print(dataset_name)

props = Deconvolve_transcriptome(
    transcriptome_data = expr_raw,
    deconvolution_algorithm = "bseqsc",
    models = model_name,
    Cibersort_absolute_mode = TRUE,
    nr_permutations = 250,
    output_file = ""
)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

props = cbind(rownames(props),rep(model_name,nrow(props)),props)
colnames(props)[1:2] = c("Sample","Model")

if ("Ductal" %in% colnames(props)){
    props_export = props[,c("Sample","Model","Alpha","Beta","Gamma","Delta","Acinar","Ductal","P_value","Correllation","RMSE")]
} else{
    props_export = props[,c("Sample","Model","Alpha","Beta","Gamma","Delta","P_value","Correlation","RMSE")]
}

#o_filename = "~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Baron_endocrine/"
o_filename = "~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Absolute/Baron_endocrine/"
o_filename = paste(o_filename, dataset_name, sep ="/")
write.table(props_export,o_filename,sep = "\t")
