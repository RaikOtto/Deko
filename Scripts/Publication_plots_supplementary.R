library("umap")
library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("dplyr")
library("grid")

meta_info_maptor = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info_maptor) = meta_info_maptor$Sample
colnames(meta_info_maptor) = str_replace(colnames(meta_info_maptor),pattern = "\\.","_")

meta_info_maptor$OS_Tissue = as.double(str_replace(meta_info_maptor$OS_Tissue,pattern = ",","."))

meta_info = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

matcher = match(meta_info_maptor$Sample,meta_info$Sample, nomatch = 0)
meta_info[matcher,"OS_Tissue"] = meta_info_maptor[matcher != 0,"OS_Tissue"]

### Supplementary Figure 4 and 5

#expr_raw = read.table("~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN/Riemer_Scarpa.S52.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
expr_raw = read.table("~/Deko_Projekt/Data/Publication_datasets/Combinations_PanNEN/Riemer_Scarpa_Master_Diedisheim.DESeq2.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "\\.", "")
expr_raw[1:5,1:5]
dim(expr_raw)

meta_data = meta_info[colnames(expr_raw),]
dim(meta_data)

cell_type_predictions = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions_visualization/Absolute/Baron_exocrine/All.exocrine.Baron.absolute.S356.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = TRUE)
rownames(cell_type_predictions) = cell_type_predictions$Sample
cell_type_predictions = cell_type_predictions[meta_data$Sample,]

source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko_Projekt/Misc/Stem_signatures.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)

i = 7

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
length(sad_genes)
meta_data = meta_info[colnames(expr_raw),]
cell_type_predictions$Exocrine_like = rowSums(cell_type_predictions[,c("Acinar","Ductal")])
meta_data[,c("Alpha","Beta","Gamma","Delta","Exocrine_like")] = cell_type_predictions[,c("Alpha","Beta","Gamma","Delta","Exocrine_like")]
meta_data$Mki_67 = log(meta_data$Mki_67)

expr = expr_raw[rownames(expr_raw) %in% sad_genes[],]
correlation_matrix = cor((expr))
meta_data$MKi67 = meta_data$Mki_67

#svg(filename = "~/Dropbox/Figures/Supplementary/SM_Figure_4_NEC_NET_PCA_new.svg", width = 10, height = 10)
p  =pheatmap::pheatmap(
  correlation_matrix,
  #expr,
  annotation_col = meta_data[,c("MKi67","NEC_NET","Grading","Primary_Metastasis","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  clustering_method = "average"
)
dev.off()

# SM Figure 4 <- here be changes for supervised versus unsupervised

#data_t = read.table("~/Deko_Projekt/Results/Figure_4.tsv",sep="\t",header = T,stringsAsFactors = F)
data_t = read.table("~/Deko_Projekt/Results/SM_Figure_3.tsv",sep="\t",header = T,stringsAsFactors = F)
data_t = reshape2::melt(data_t)
data_t$variable = str_replace(data_t$variable,"\\.","")
data_t$Type = factor(data_t$Type)
data_t = data_t %>% dplyr::rename("Parameter" = "variable")
data_t = data_t %>% dplyr::filter(Parameter %in% c("Accuracy"))
data_t$Dataset = factor(data_t$Dataset, levels = c("Missiaglia","Scarpa","RepSet","Riemer","Sadanandam","Average"))
data_t$Type = factor(data_t$Type, levels = c("Unsupervised","Supervised"))

p = ggplot( 
  data = data_t,
  aes( 
    x = Dataset,
    y = value,
    fill = Type
  )
)
p = p + geom_bar(stat="identity", position=position_dodge())
#p = p + theme(axis.text.x = element_text(angle = 45, vjust = .5))
p = p + xlab("Dataset") + ylab("Accuracy") + theme(legend.position = "top")
p = p + scale_fill_manual(values = c("blue","red"))
p = p + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))


#svg(filename = "~/Deco/Results/Images/SM_Figure_4.svg", width = 10, height = 10)
p 
dev.off()

#### PCA Ductal Hisc gene signature

props = read.table("~/Deko_Projekt/Results/All.S200.CIBERSORT.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T)
rownames(props) = props$Sample
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

colnames(expr_raw)[colnames(expr_raw) %in% rownames(props) == F] = paste("X",colnames(expr_raw)[colnames(expr_raw) %in% rownames(props) == F],sep ="")
colnames(expr_raw)[colnames(expr_raw) %in% rownames(props) == F] = str_replace(colnames(expr_raw)[colnames(expr_raw) %in% rownames(props) == F],pattern ="^XX","")
colnames(expr_raw)[colnames(expr_raw) %in% rownames(props) == F]

matcher = match( colnames(expr_raw), rownames(props), nomatch = 0)
props = props[matcher,]
dim(props)
meta_data = meta_info[rownames(props),]

###

props = read.table("~/Deko_Projekt/Results/All.S200.CIBERSORT.tsv",sep ="\t", header = T, as.is=TRUE)
rownames(props ) = props$Sample

no_matcher = which(!(colnames(expr_raw) %in% rownames(props)))
colnames(expr_raw)[no_matcher] = str_replace(colnames(expr_raw)[no_matcher], pattern ="^X","")
no_matcher = which(!(colnames(expr_raw) %in% rownames(props)))
no_matcher
props = props[colnames(expr_raw),]

props %>% filter()

selection = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal")
exocrines = as.double(rowSums(props[,c("Ductal","Acinar")]))
endocrines = as.double(rowSums(props[,c("Alpha","Beta","Gamma","Delta")]))

meta_data$Ratio = log((exocrines+.1) / (endocrines+.1))

expr = cbind(props[,c(selection,"P_value","Correlation")],meta_data$Ratio)
expr = matrix(as.double(as.character(unlist(expr))), ncol = 9,nrow = nrow(expr))
colnames(expr) = c(selection,"P_value","Correlation","Ratio")
rownames(expr) = props$Sample

correlation_matrix = cor(t(expr))
pcr = prcomp(t(correlation_matrix))

meta_data$P_value = props$P_value
meta_data$P_value[meta_data$P_value >= 0.05] = "not_sig"
meta_data$P_value[meta_data$P_value != "not_sig"] = "sig"
rownames(meta_data)[!(rownames(meta_data) %in% colnames(correlation_matrix))] = str_replace(rownames(meta_data)[!(rownames(meta_data) %in% colnames(correlation_matrix))], pattern ="^X","")
rownames(meta_data)[!(rownames(meta_data) %in% colnames(correlation_matrix))] 

#svg(filename = "~/Downloads/Heatmap.svg", width = 10.5, height = 10)
p  = pheatmap::pheatmap(
  correlation_matrix,
  annotation_col = meta_data[,c("Grading","NEC_NET","Study")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  treeheight_row = 0,
  legend = T,
  fontsize_col = 7,
  clustering_method = "complete"
)
dev.off()

###

meta_t = read.table("~/Deko_Projekt/GSEA/metaplastic_and cancer signatures from Schlesinger_Zhibo_Hendley.transpose.gmt.tsv",sep ="\t", header = F)
dim(meta_t)

genes_9 = c( genes_9, rep("",ncol(meta_t)-length(genes_9)))
genes_10 = c( genes_10, rep("",ncol(meta_t)-length(genes_10)))
genes_11 = c( genes_11, rep("",ncol(meta_t)-length(genes_11)))
genes_14 = c( genes_14, rep("",ncol(meta_t)-length(genes_14)))
genes_15 = c( genes_15, rep("",ncol(meta_t)-length(genes_15)))
genes_17 = c( genes_17, rep("",ncol(meta_t)-length(genes_17)))

stem_t = rbind(stem_t,genes_10)
stem_t = rbind(stem_t,genes_11)
stem_t = rbind(stem_t,genes_14)
stem_t = rbind(stem_t,genes_15)
stem_t = rbind(stem_t,genes_17)


meta_t[,1:10]
meta_t = rbind(meta_t, genes_1, genes_9, genes_10, genes_11, genes_14,genes_15,genes_17)

write.table(meta_t,"~/Deko_Projekt/GSEA//Stem_signatures.gmt.tsv",sep ="\t",quote = F, row.names = FALSE,col.names = FALSE)

#### umap

### UMAP

props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/RepSet_S103.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)

#props = read.table("~/Deko_Projekt/Results/Cell_fraction_predictions/Alvarez.S104.Cibersort.tsv",sep = "\t", as.is = T, stringsAsFactors = F, header = T,row.names = 1)
colnames(props)[colnames(props) == "alpha"] = "Alpha";colnames(props)[colnames(props) == "beta"] = "Beta";colnames(props)[colnames(props) == "gamma"] = "Gamma";colnames(props)[colnames(props) == "delta"] = "Delta";colnames(props)[colnames(props) == "acinar"] = "Acinar";colnames(props)[colnames(props) == "ductal"] = "Ductal"
colnames(props) = str_replace(colnames(props),pattern = "\\.","-")
colnames(props) = str_replace(colnames(props),pattern = "\\.ductal","ductal")
colnames(props) = str_replace(colnames(props),pattern = "-reg\\.","-reg+")
colnames(props) = str_replace(colnames(props),pattern = "muc5b-","muc5b+-")

no_match = rownames(props) %in% meta_info$Sample == F
rownames(props)[no_match] = paste("X",rownames(props)[no_match],sep ="")
no_match = rownames(props) %in% meta_info$Sample == F
sum(no_match)

dim(props)
meta_data = meta_info[rownames(props),]

selection = colnames(props)[!(colnames(props) %in% c("model","Sig_score","P_value","Correlation","RMSE","Subtype","Strength_subtype","hisc"))]
props = as.data.frame(props)
vis_mat = props[,colnames(props) %in% selection]
vis_mat = vis_mat[meta_data$NET_NEC_PCA != "Unknown",]
meta_data = meta_info[rownames(vis_mat),]

correlation_matrix = cor(t(vis_mat));pcr = prcomp(t(correlation_matrix))
#vis_mat = vis_mat[order(vis_mat$endocrine_fully_differentited,decreasing = T),]

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data$NET_NEC_PCA,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$NET_NEC_PCA) ))
#umap_p = umap_p+geom_text(size= 4,aes(label=rownames(meta_data),color = as.character(meta_data$Cluster)),hjust=0, vjust=0)
umap_p = umap_p +  xlab("") + ylab("")  
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$NET_NEC_PCA), level=.5, type ="t", size=1.5)
#umap_p = umap_p + scale_color_manual( values = c("darkgreen","orange"),name="Drug_treatment") ##33ACFF ##FF4C33
umap_p
genes_of_interest_hgnc_t$V1[i]

###

meta_data = meta_info[meta_info$Study!="Fadista",]
meta_data = meta_data[meta_data$Grading != "Control",]
meta_data = meta_data[meta_data$Primary_Metastasis != "Outlier",]
meta_data_sato = meta_data[meta_data$Study == "Sato",]
meta_data_sato = meta_data_sato[meta_data_sato$NEC_NET != "MiNEN",]
meta_data_sato_pan = meta_data_sato[meta_data_sato$Site_of_primary == "Pancreatic",]
meta_data_sato_gep = meta_data_sato[meta_data_sato$Site_of_primary != "Pancreatic",]
#meta_data = meta_data[meta_data$Study != "Sato",]
meta_data_pan = meta_data[meta_data$Site_of_primary == "Pancreatic",]
meta_data_gep = meta_data[meta_data$Site_of_primary != "Pancreatic",]

table(meta_data_pan[meta_data_pan$NEC_NET == "NET","Study"],meta_data_pan[meta_data_pan$NEC_NET == "NET","Grading"])
table(meta_data_pan[meta_data_pan$NEC_NET == "NEC","Study"],meta_data_pan[meta_data_pan$NEC_NET == "NEC","Grading"])
table(meta_data_pan[meta_data_pan$NEC_NET == "Ambiguous","Study"],meta_data_pan[meta_data_pan$NEC_NET == "Ambiguous","Grading"])

table(meta_data_pan$Grading)
table(meta_data_pan[meta_data_pan$NEC_NET == "NET","Grading"])
table(meta_data_pan[meta_data_pan$NEC_NET == "NEC","Grading"])
table(meta_data_pan[meta_data_pan$NEC_NET == "Ambiguous","Grading"])

table(meta_data_gep[meta_data_gep$NEC_NET == "NET","Study"],meta_data_gep[meta_data_gep$NEC_NET == "NET","Grading"])
table(meta_data_gep[meta_data_gep$NEC_NET == "NEC","Study"],meta_data_gep[meta_data_gep$NEC_NET == "NEC","Grading"])
table(meta_data_gep[meta_data_gep$NEC_NET == "Ambiguous","Study"],meta_data_gep[meta_data_gep$NEC_NET == "Ambiguous","Grading"])

table(meta_data_gep$Grading)
table(meta_data_gep[meta_data_gep$NEC_NET == "NET","Grading"])
table(meta_data_gep[meta_data_gep$NEC_NET == "NEC","Grading"])

table(meta_data[meta_data$NEC_NET == "NET","Study"],meta_data[meta_data$NEC_NET == "NET","Grading"])
table(meta_data[meta_data$NEC_NET == "NEC","Grading"],meta_data[meta_data$NEC_NET == "NEC","Study"])
table(meta_data[meta_data$NEC_NET == "Ambiguous","Grading"],meta_data[meta_data$NEC_NET == "Ambiguous","Study"])

table(meta_data_sato_pan$Grading,meta_data_sato_pan$NEC_NET)
table(meta_data_sato_gep$Grading,meta_data_sato_gep$NEC_NET)
