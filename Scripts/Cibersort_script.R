library("bseqsc")
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
#library("pheatmap")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T)
#bam_data = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep ="\t", header = T)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")

condition = (meta_info$Included == "Yes") & (meta_info$NEC_NET %in% c("NEC","NET")) #& (meta_info$Grading %in% c("G3"))
include_list  = meta_info$Name[ condition ]
exclude_list  = meta_info$Name[ ! condition ]

bam_data = bam_data[, ( colnames(bam_data) %in% include_list )]
dim(bam_data)
bam_data[1:5,1:5]

#count_data = read.table("~/Deko/Data/baron_data.tsv",sep ="\t", header = T, stringsAsFactors = F)
count_data = read.table("~/Deko/Data/Count_data.Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)

#subtypes = read.table(  "~/Deko/Data/baron_subtypes.tsv",sep ="\t", header = T, stringsAsFactors = F)[,1]
seg_meta = read.table("~/Deko/Data/Segerstolpe_Meta_info.tsv", sep ="\t", header = T)
s_match = match(colnames(s_t), seg_meta$Extract.Name, nomatch = 0)

subtypes = as.character(seg_meta$Characteristics.cell.type.)[s_match]
subtypes = str_replace_all(subtypes, pattern = " cell", "")
table(subtypes)

#s_match = match(colnames(count_data), meta_info$Name, nomatch = 0)
#meta_data = meta_info[s_match,]
#rownames(meta_data) = meta_data$Name

### normalization

row_var = apply(count_data, FUN = var, MARGIN = 1)
col_var = apply(count_data, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
count_data = count_data[row_var != 0,col_var != 0]
count_data = count_data[rowSums(count_data) >= 1,]
count_data = t(count_data)
dim(count_data)


### load_data

marker_genes = read.table(
  "~/Deko/Data/Baron_pancreas_marker.tsv",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

pancreasMarkers = list(
  "Alpha" = marker_genes$alpha[marker_genes$alpha != ""],
  "Beta" = marker_genes$beta[marker_genes$beta != ""],
  "Gamma" = marker_genes$gamma[marker_genes$gamma != ""],
  "Delta" = marker_genes$delta[marker_genes$delta != ""],
  "Ductal" = marker_genes$ductal[marker_genes$ductal != ""],
  "Acinar" = marker_genes$acinar[marker_genes$acinar != ""]
)

# count data baron

# count_data based training

table(str_to_upper(subtypes))
table(str_to_upper(names(pancreasMarkers)))
count_data = count_data[,str_to_upper(subtypes) %in% str_to_upper(names(pancreasMarkers)) ]
subtypes = subtypes[str_to_upper(subtypes) %in% str_to_upper(names(pancreasMarkers)) ]
dim(count_data)
length(subtypes)

eislet = new("ExpressionSet", exprs = as.matrix(count_data))
subtypes[subtypes == "alpha"] = "Alpha"
subtypes[subtypes == "beta"] = "Beta"
subtypes[subtypes == "gamma"] = "Gamma"
subtypes[subtypes == "delta"] = "Delta"
subtypes[subtypes == "acinar"] = "Acinar"
subtypes[subtypes == "ductal"] = "Ductal"

fData(eislet) = data.frame( subtypes )
pData(eislet) = data.frame( subtypes )

B = bseqsc_basis(
  eislet,
  pancreasMarkers,
  clusters = 'subtypes',
  samples = colnames(exprs(eislet)),
  ct.scale = FALSE
)

### run

eset = new("ExpressionSet", exprs=as.matrix(bam_data));
fit = bseqsc_proportions(eset, B, verbose = TRUE)

###

res_coeff = t(fit$coefficients)
res_cor   = fit$stats

###

meta_match = match( colnames(exprs(eset)), meta_info$Name, nomatch = 0 ) 
res_coeff_match = match( rownames(res_coeff), meta_info$Name, nomatch = 0 )
res_cor_match = match( rownames(res_cor), meta_info$Name, nomatch = 0 )

maxi = apply( res_coeff, FUN = which.max, MARGIN = 1 )
cell_type = colnames(res_coeff)[maxi]
#cell_type[res_cor[,"Correlation"] <= .3] = "Not_sig"
meta_data$Cell_type = cell_type[match(rownames(meta_data), rownames(res_coeff))]

meta_info$Alpha = rep("",nrow(meta_info))
meta_info$Beta = rep("",nrow(meta_info))
meta_info$Gamma = rep("",nrow(meta_info))
meta_info$Delta = rep("",nrow(meta_info))
meta_info$Acinar = rep("",nrow(meta_info))
meta_info$Ductal = rep("",nrow(meta_info))
meta_info$Res_cor = rep("",nrow(meta_info))
meta_info$Alpha[res_coeff_match] = as.data.frame(res_coeff)$Alpha
meta_info$Beta[res_coeff_match] = as.data.frame(res_coeff)$Beta
meta_info$Gamma[res_coeff_match] = as.data.frame(res_coeff)$Gamma
meta_info$Delta[res_coeff_match] = as.data.frame(res_coeff)$Delta
meta_info$Acinar[res_coeff_match] = as.data.frame(res_coeff)$Acinar
meta_info$Ductal[res_coeff_match] = as.data.frame(res_coeff)$Ductal
meta_info$Res_cor[res_cor_match] = as.data.frame(res_cor)$Correlation

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv",sep="\t", row.names = F, quote = F)

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0
res_coeff[ res_cor[,"Correlation"] <= .3, ] = .0
res_cor[ res_cor[,"Correlation"] <= .3, "Correlation" ] = 0.0

pheatmap::pheatmap(
  t(res_coeff),
  annotation_col = meta_data[c("Subtype_Sadanandam","MKI67","NEC_NET","Grading","Location","Histology","Study")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  show_colnames = T,
  show_rownames = T,
  gaps_row = 20
)

