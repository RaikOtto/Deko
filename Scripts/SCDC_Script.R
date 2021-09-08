library("SCDC")
library("stringr")
library("reshape2")
library("dplyr")
library("Biobase")

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info = read.table("~/MAPTor_NET///Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$Subtype
res_scdc = as.data.frame(meta_info)

#expr_raw = read.table("~/Deko_Projekt/Data/JGA/Sato.S35.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
expr_raw = read.table("~/Deko_Projekt/Data/Diedisheim.S66.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
no_match

meta_data = meta_info[colnames(expr_raw),]
table(meta_data$Study)

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
expr_scrna =  as.data.frame(read.table("~/Deko_Projekt//Data/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.tsv", sep ="\t", header = T))

cell_type_vec = meta_info[colnames(expr_scrna),"NEC_NET_Color"]
#expr_scrna = expr_scrna[,!(cell_type_vec %in% c("Ductal","Acinar"))]
table(cell_type_vec)

fdata = rownames(expr_scrna)
pdata = cbind(cellname = colnames(expr_scrna), subjects = cell_type_vec)
eset_scrna = getESET(expr_scrna, fdata = fdata, pdata = pdata)
eset_scrna$Subtype = meta_info[eset_scrna$cellname,"NEC_NET_Color"]

sample_id = rep("",length(eset_scrna$Subtype))
sample_id[grep(colnames(expr_scrna),pattern = "human1",value = F)] = "human1"
sample_id[grep(colnames(expr_scrna),pattern = "human2",value = F)] = "human2"
sample_id[grep(colnames(expr_scrna),pattern = "human3",value = F)] = "human3"
sample_id[grep(colnames(expr_scrna),pattern = "human4",value = F)] = "human4"
table(sample_id)
eset_scrna$Sample = sample_id

scrna.qc = SCDC_qc(
    eset_scrna,
    ct.varname = "Subtype",
    sample = "Sample",
    scsetname = "scRNA",
    ct.sub = unique(eset_scrna$Subtype),
    qcthreshold = 0.7
)
#DemoPlot(eset_scrna, cluster = "Subtype", sample = "Sample", select.ct = unique(eset_scrna$Subtype))
scrna.qc$heatfig

# prediction

scdc_props = SCDC_prop(
    bulk.eset = eset_expr_raw,
    sc.eset = eset_scrna,
    ct.varname = "Subtype",
    sample = "Sample",
    ct.sub = unique(eset_scrna$Subtype),
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
#colnames(props) = paste("SCDC",colnames(scdc_props$prop.est.mvw),sep = "_")
rownames(props)  = rownames(scdc_props$prop.est.mvw) 

#write.table(deconvolution_results,"~/Deko_Projekt/Results/Cell_fraction_predictions/Sato.S35.SCDC.tsv",sep = "\t")
props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Diedisheim_S66.tsv",sep = "\t",as.is = F, stringsAsFactors = F)
###
library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

deconvolution_results = Deconvolve_transcriptome(
    transcriptome_data = expr_raw[,],
    deconvolution_algorithm = "bseqsc",
    #models = "Alpha_Beta_Gamma_Delta_Baron",
    models = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
    Cibersort_absolute_mode = TRUE,
    nr_permutations = 1000,
    output_file = ""
)

#write.table(deconvolution_results,"~/Deko_Projekt/Results/Cell_fraction_predictions/Sato.S29.Ninon.tsv",sep = "\t")

#props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Archive/RepSet_Cibersort_Baron.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Alvarez.S104.Cibersort.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[rownames(props),]

###

props = props[(meta_data$Histology == "pancreas") | (meta_data$NEC_NET_Color != "Primary") ,]
meta_data = meta_info[rownames(props),]
props = props[(meta_data$Study == "Alvarez") ,]
meta_data = meta_info[rownames(props),]

selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")
exocrines = as.double(rowSums(props[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(props[,c("Alpha","Beta","Gamma","Delta")]))

meta_data$Ratio = log((exocrines+.1) / (endocrines+.1))
#meta_data$Ratio = ((exocrines+.1) / (endocrines+.1))
#meta_data[,selection] = props[,selection]

#matcher = match(rownames(meta_data),rownames(meta_info))
#meta_info[ matcher,"Ratio"] = meta_data$Ratio

#meta_info[rownames(meta_data), selection] = meta_data[,selection]
#write.table(meta_info,"~/Deko_Projekt/Misc/Meta_information.tsv",sep ="\t",quote =F,row.names = F)
###

vis_mat = props[,selection]
vis_mat$Exocrine = vis_mat$Acinar + vis_mat$Ductal
vis_mat = vis_mat[,!(colnames(vis_mat) %in% c("Acinar","Ductal"))]
#vis_mat = props[ !(rownames(props) %in% outlier),selector]
#vis_mat[,"Ratio"] = as.double(meta_data[rownames(vis_mat),"Ratio"])

correlation_matrix = cor(t(vis_mat))

pcr = prcomp(t(correlation_matrix))

p = pheatmap::pheatmap(
    t(vis_mat),
    annotation_col = meta_data[,c("Grading","Cluster","Functionality","NEC_NET_Ori")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = F,
    treeheight_row = 0,
    cluster_rows = F,
    legend = F,
    fontsize_row = 14,
    clustering_method = "average"
)
p +  theme(legend.position="top",axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.text=element_text(size=18),legend.title=element_text(size=18))

### PCA

p = ggbiplot::ggbiplot(
    prcomp(t(vis_mat)),
    groups = as.character(meta_data$Grading),
    var.axes = F,
    ellipse = TRUE
)
ki_67_vec = as.double(expr_raw["MKI67",])*.5

p = p + geom_point( aes( shape = meta_data$Grading, color = meta_data$Grading ), size = ki_67_vec ) # Fig 4
p = p + scale_color_manual( values = c("darkgreen","yellow","red") ) #Fig 4 Master
p = p + guides(fill=FALSE) + scale_fill_discrete(guide=FALSE)+ theme(legend.position="none")
p

