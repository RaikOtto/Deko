library("SCDC")
library("stringr")
library("reshape2")
library("dplyr")
library("Biobase")

meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$Subtype
res_scdc = as.data.frame(meta_info)

#meta_info$Patient = rep("",nrow(meta_info))
#h1_detect = str_detect(meta_info$Name,pattern ="human1")
#h2_detect = str_detect(meta_info$Name,pattern ="human2")
#h3_detect = str_detect(meta_info$Name,pattern ="human3")
#h4_detect = str_detect(meta_info$Name,pattern ="human4")
#meta_info[h1_detect,"Patient"] = "human_1"
#meta_info[h2_detect,"Patient"] = "human_2"
#meta_info[h3_detect,"Patient"] = "human_3"
#meta_info[h4_detect,"Patient"] = "human_4"

#write.table(meta_info,"~/Deko_Projekt/Misc/Meta_information.tsv", sep ="\t")

expr_raw = read.table("~/MAPTor_NET/BAMs_new/RepSet_S84.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
colnames(expr_raw)[no_match] = paste("X",colnames(expr_raw)[no_match],sep ="")
no_match = colnames(expr_raw) %in% meta_info$Sample == F
no_match

meta_data = meta_info[colnames(expr_raw),]
meta_data$Sample
#expr_raw = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_Bulk/GSE142720_rma_norm_log2_matrix.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

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

cell_type_vec = meta_info[colnames(expr_scrna),"Subtype"]
#expr_scrna = expr_scrna[,!(cell_type_vec %in% c("Ductal","Acinar"))]
table(cell_type_vec)

fdata = rownames(expr_scrna)
pdata = cbind(cellname = colnames(expr_scrna), subjects = cell_type_vec)
eset_scrna = getESET(expr_scrna, fdata = fdata, pdata = pdata)
eset_scrna$Subtype = meta_info[eset_scrna$cellname,"Subtype"]

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

#write.table(deconvolution_results,"~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet.S84.Cibersort.tsv",sep = "\t")
props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Alvarez.S104.SCDC.tsv",sep = "\t",as.is = F, stringsAsFactors = F)
###
library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

deconvolution_results = Deconvolve_transcriptome(
    transcriptome_data = expr_raw,
    deconvolution_algorithm = "bseqsc",
    models = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
    nr_permutations = 1000,
    output_file = ""
)

#write.table(deconvolution_results,"~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet.S80.exocrine.CIBERSORT.tsv",sep = "\t")
props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet.S84.Cibersort.tsv",sep = "\t", as.is = T, stringsAsFactors = F)
rownames(props) = str_replace(rownames(props), pattern = "^X", "")
no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
meta_data = meta_info[rownames(props),]

###

colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")
exocrines = as.double(rowSums(props[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(props[,c("Alpha","Beta","Gamma","Delta")]))

meta_data$Ratio = log((exocrines+.1) / (endocrines+.1))
#meta_data$Ratio = ((exocrines+.1) / (endocrines+.1))
meta_data[,selection] = props[,selection]
#max_selection = max(meta_data[,selection])
meta_data[,selection] = meta_data[,selection] / max(meta_data[,selection])

#meta_data$Ratio = meta_data$Ratio/max(meta_data$Ratio)

###


for ( i in 1:nrow(meta_data)){
    meta_data[i,selection] = meta_data[i,selection] / max(meta_data[i,selection] )
    #res_scdc[rownames(props),colname] = as.double(res_scdc[rownames(props),colname])
}
#### visualization

expr = expr[,str_detect(colnames(expr),pattern = "_",negate = T)]
vis_mat = vis_mat[str_detect(rownames(vis_mat),pattern = "_",negate = T),]

correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

#svg(filename = "~/Deko_Projekt/Results/Images/SM_Figure_4_Correlation_Heatmap_RepSet.svg", width = 10, height = 10)
pheatmap::pheatmap(
    correlation_matrix,
    #annotation_col = vis_mat[,c("Alpha","SCDC_Alpha","Beta","SCDC_Beta","Gamma","SCDC_Gamma","Delta","SCDC_Delta","Acinar","SCDC_Acinar","Ductal","SCDC_Ductal","Grading")],
    #annotation_col = vis_mat[,c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","Grading")],
    annotation_col = vis_mat[,c("alpha","beta","gamma","delta","acinar","ductal","NEC_NET","Grading")],
    annotation_colors = aka3,
    show_rownames = F,
    show_colnames = F,
    #treeheight_col = 0,
    treeheight_row = 0,
    legend = T,
    fontsize_col = 7,
    clustering_method = "median"
)

#

cell_m = as.data.frame(props[colnames(expr_raw),])
cell_m$MKI67 = as.double(round(expr_raw["MKI67",rownames(cell_m)] / max(expr_raw["MKI67",rownames(cell_m)]) * 100,1))
cell_m$Sample = rownames(cell_m)
cell_m$Grading = meta_data$Grading
cell_m$NEC_NET = meta_data$Subtype
cell_m$Histology = meta_data$Histology
cell_m$Grading[(cell_m$NEC_NET == "NEC") & (cell_m$Grading == "G3")] = "G3_NEC"
cell_m$Grading[(cell_m$NEC_NET == "NET") & (cell_m$Grading == "G3")] = "G3_NET"

cell_m = cell_m %>% filter(str_detect(cell_m$Sample, pattern = "_",negate = T))

cell_m_exo = cell_m[,c("Grading","NEC_NET","alpha","beta","gamma","delta","acinar","ductal")] %>% melt() 
colnames(cell_m_exo) = c("Grading","NEC_NET","Celltype","Proportion")
cell_m_exo = cell_m_exo %>% filter(!( Celltype %in%  c("MKI67","P_value")))
cell_m_exo$Celltype = sapply( as.character(cell_m_exo$Celltype), FUN = function(vec){return(tail(as.character(unlist(str_split(vec,pattern = "_"))),1))})

cell_m_exo_g1 = cell_m_exo[cell_m_exo$Grading == "G1",]
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] + 1
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_exo_g1 = aggregate(cell_m_exo_g1$Proportion, by = list(cell_m_exo_g1$Celltype), FUN = sum)
vis_mat_exo_g1$x = round(vis_mat_exo_g1$x / sum(vis_mat_exo_g1$x) * 100, 1 )
vis_mat_exo_g1$Grading = rep("G1",nrow(vis_mat_exo_g1))
cell_m_exo_g2 = cell_m_exo[cell_m_exo$Grading == "G2",]
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] + .5
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_exo_g2 = aggregate(cell_m_exo_g2$Proportion, by = list(cell_m_exo_g2$Celltype), FUN = sum)
vis_mat_exo_g2$x = round(vis_mat_exo_g2$x / sum(vis_mat_exo_g2$x)  * 100, 1 )
vis_mat_exo_g2$Grading = rep("G2",nrow(vis_mat_exo_g2))

cell_m_exo_g3_NET = cell_m_exo[ (cell_m_exo$Grading == "G3_NET"),]
vis_mat_exo_g3_NET = aggregate(cell_m_exo_g3_NET$Proportion, by = list(cell_m_exo_g3_NET$Celltype), FUN = sum)
vis_mat_exo_g3_NET$x = round(vis_mat_exo_g3_NET$x / sum(vis_mat_exo_g3_NET$x)  * 100, 1 )
vis_mat_exo_g3_NET$Grading = rep("G3_NET",nrow(vis_mat_exo_g3_NET))
cell_m_exo_g3_NEC = cell_m_exo[ (cell_m_exo$Grading == "G3_NEC"),]
vis_mat_exo_g3_NEC = aggregate(cell_m_exo_g3_NEC$Proportion, by = list(cell_m_exo_g3_NEC$Celltype), FUN = sum)
vis_mat_exo_g3_NEC$x = round(vis_mat_exo_g3_NEC$x / sum(vis_mat_exo_g3_NEC$x)  * 100, 1 )
vis_mat_exo_g3_NEC$Grading = rep("G3_NEC",nrow(vis_mat_exo_g3_NEC))

vis_mat_exo = rbind(vis_mat_exo_g1,vis_mat_exo_g2,vis_mat_exo_g3_NET,vis_mat_exo_g3_NEC)
colnames(vis_mat_exo) = c("Celltype","Proportion","Grading")
vis_mat_exo$Grading = factor(vis_mat_exo$Grading, levels = c("G1","G2","G3_NET","G3_NEC"))

library("ggplot2")
p_exo = ggplot(
    data = vis_mat_exo,
    aes(
        x = Grading,
        y = Proportion
    )
) + geom_bar(
    aes(
        y = Proportion,
        x = Grading,
        fill = Celltype
    ),
    data = vis_mat_exo,
    stat="identity",
    colour="black"
) + scale_fill_manual(values = c("cyan", "blue","yellow","purple","darkred","orange")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_exo
