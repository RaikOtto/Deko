library("devtools")
library("NMF")
library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")
library("SCDC")
library("stringr")
library("reshape2")
library("bseqsc")
library("dplyr")
library("Biobase")

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S56.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_Bulk/Fadista.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
#expr_raw = read.table("~/Dropbox/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep="\t", stringsAsFactors =  F, header = T,as.is = TRUE)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

no_match = expr_raw$Sample %in% meta_info$Sample == F
expr_raw$Sample[no_match] = paste("X",expr_raw$Sample[no_match],sep ="")
no_match = expr_raw$Sample %in% meta_info$Sample == F
table(no_match)

#expr_raw = expr_raw[,meta_data$Histology_Primary %in% c("Pancreatic")]
meta_data = meta_info[colnames(expr_raw),]

table(meta_data$Study)
table(meta_data$NET_NEC_PCA)

fdata = rownames(expr_raw)
pdata = cbind(bulk_sample = colnames(expr_raw))
eset_expr_raw = getESET(expr_raw, fdata = fdata, pdata = pdata)
eset_expr_raw$Grading = meta_data[eset_expr_raw$bulk_sample,"Grading"]

#
source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
i = 13
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
expr[1:5,1:5]
dim(expr)

# ScRNA EXO 

#expr_scrna =  read.table("~/Deko_Projekt//Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.tsv", sep ="\t", header = T)
expr_scrna =  as.data.frame(read.table("~/Dropbox/Tosti.S3573.RepSet_genes.tsv", sep ="\t", header = T))

ell_type_vec = meta_info[colnames(expr_scrna),"NEC_NET_Color"]
cell_type_vec
#expr_scrna = expr_scrna[,!(cell_type_vec %in% c("Ductal","Acinar"))]
table(cell_type_vec)

fdata = rownames(expr_scrna)
pdata = cbind(cellname = colnames(expr_scrna), cell_type_vec = cell_type_vec)
eset_scrna = getESET(expr_scrna, fdata = fdata, pdata = pdata)
#eset_scrna$Subtype = meta_info[eset_scrna$cellname,"NEC_NET_Color"]


eset_scrna$Donor_id = meta_data$patient_ID
table(eset_scrna$Donor_id)


scrna.qc = SCDC_qc(
    eset_scrna,
    ct.varname = "cell_type_vec",
    sample = "Donor_id",
    scsetname = "scRNA",
    ct.sub = unique(eset_scrna$cell_type_vec),
    qcthreshold = 0.7
)
#DemoPlot(eset_scrna, cluster = "Subtype", sample = "Sample", select.ct = unique(eset_scrna$Subtype))
scrna.qc$heatfig

# prediction

scdc_props = SCDC_prop(
    bulk.eset = eset_expr_raw,
    sc.eset = eset_scrna,
    ct.varname = "cell_type_vec",
    sample = "Donor_id",
    ct.sub = unique(eset_scrna$cell_type_vec),
    iter.max = 1000,
    nu = 1e-04,
    epsilon = 0.01,
    truep = NULL,
    weight.basis = T,
    ct.cell.size = NULL,
    Transform_bisque = F
)

props = matrix(scdc_props$prop.est.mvw,nrow = nrow(scdc_props$prop.est.mvw))
colnames(props) = colnames(scdc_props$prop.est.mvw)
rownames(props)  = rownames(scdc_props$prop.est.mvw)
rownames(props) = str_replace(rownames(props), pattern = "^X","")

###

expr_raw = read.table("~/Dropbox/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.tsv",sep="\t", stringsAsFactors =  F, header = T,as.is = TRUE)
expr_raw = read.table("~/Dropbox/Cell_fraction_predictions_visualization/All.endocrine.exocrine.Baron.absolute.with_NENs.tsv",sep="\t", stringsAsFactors =  F, header = T,as.is = TRUE)
colnames(expr_raw)

meta_data = meta_info[expr_raw$Sample,]
expr_raw = expr_raw[meta_data$Grading != "Unknown",]
meta_data = meta_info[expr_raw$Sample,]
expr_raw = expr_raw[!( meta_data$NET_NEC_PCA %in% c("Unknown","MiNEN")),]
meta_data = meta_info[expr_raw$Sample,]
Grading_binary = meta_data$Grading
Grading_binary[Grading_binary %in% c("G1","G2")] = "G1_G2"

table(meta_data$NET_NEC_PCA)
expr_raw$Grading_binary = Grading_binary
expr_raw$NET_NEC = meta_data$NET_NEC_PCA
expr_raw = expr_raw[expr_raw$Model %in% c("Endocrine_exocrine_like","Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"),]
rownames(expr_raw) = expr_raw$Sample
expr_raw = expr_raw[ ,!(colnames(expr_raw) %in% c("Sample","Model"))]
dim(expr_raw)
table(expr_raw$NET_NEC)
meta_data = meta_info[rownames(expr_raw),]

expr_raw$Grading_binary[expr_raw$Grading_binary == "G1_G2"] = 1
expr_raw$Grading_binary[expr_raw$Grading_binary == "G3"] = 0

expr_raw_repset = expr_raw[meta_data$Study %in% c("Charite","Master","Scarpa"),]
dim(expr_raw_repset)

expr_raw_g3_only = expr_raw[meta_data$Grading %in% c("G3"),]
dim(expr_raw_g3_only)

write.table(expr_raw,"~/Dropbox/testproject/Datasets/Deconvolution.Absolute.Exocrine.PanNEN.NEN.S281.tsv",sep = "\t")
table(meta_data[,c("Study","Grading")])

###
props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S103.tsv",sep = "\t",as.is = F, stringsAsFactors = F)

show_models_bseqsc()
model_name = "Alpha_Beta_Gamma_Delta_Baron"


props = Deconvolve_transcriptome(
    transcriptome_data = expr_raw,
    deconvolution_algorithm = "bseqsc",
    models = model_name,
    Cibersort_absolute_mode = FALSE,
    nr_permutations = 1000,
    output_file = ""
)

#write.table(props,"~/Deko_Projekt/Results/Cell_fraction_predictions/Fadista_S89_Cibersort_Baron_metaplastic_only.tsv",sep = "\t")

