library("bseqsc")
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
#library("pheatmap")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

bam_data = read.table("~/Deko/Data/GSE73338.HGNC.tsv",sep ="\t", header = T)
bam_data = read.table("~/Deko/Data/TPMs.Not_normalized.Controls_Groetzinger_Scarpa.88S.tsv",sep ="\t", header = T)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")
bam_data[1:5,1:5]

#condition = (meta_info$Included == "Yes") & (meta_info$NEC_NET %in% c("NEC","NET")) #& (meta_info$Grading %in% c("G3"))
#include_list  = meta_info$Name[ condition ]
#exclude_list  = meta_info$Name[ ! condition ]

#bam_data = bam_data[, ( colnames(bam_data) %in% include_list )]
#dim(bam_data)
#bam_data[1:5,1:5]

count_data = read.table("~/Deko/Data/Count_data.Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)
count_data[1:5,1:5]

gene_length_t = read.table("~/Deko/Misc/gene_length.tsv",sep ="\t", header = T, stringsAsFactors = F)
l_match = match(rownames(count_data), gene_length_t$hgnc_symbol, nomatch = 0)
count_data = count_data[ l_match != 0,]
l_match = match( gene_length_t$hgnc_symbol, rownames(count_data), nomatch = 0)
gene_length = gene_length_t$transcript_length[l_match]
x = count_data / gene_length
count_data = t(x) * 1e6 / colSums(x)

seg_meta = read.table("~/Deko/Data/Segerstolpe_Meta_info.tsv", sep ="\t", header = T)

count_data = t(count_data)
s_match = match( colnames(count_data), seg_meta$Extract.Name, nomatch = 0)

subtypes = as.character(seg_meta$Characteristics.cell.type.)[s_match]
subtypes = str_replace_all(subtypes, pattern = " cell", "")
table(subtypes)

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

#pancreasMarkers = list("Alpha" = marker_genes$alpha[marker_genes$alpha != ""],"Beta" = marker_genes$beta[marker_genes$beta != ""],"Gamma" = marker_genes$gamma[marker_genes$gamma != ""],"Delta" = marker_genes$delta[marker_genes$delta != ""],"Ductal" = marker_genes$ductal[marker_genes$ductal != ""],"Acinar" = marker_genes$acinar[marker_genes$acinar != ""])
pancreasMarkers = list(
  "Alpha" = marker_genes$alpha[1:3],
  "Beta" = marker_genes$beta[1:3],
  "Gamma" = marker_genes$gamma[1:3],
  "Delta" = marker_genes$delta[1:3]
  #"Ductal" = marker_genes$ductal[1],
  #"Acinar" = marker_genes$acinar[1]
)

# count data baron

# count_data based training

#table(str_to_upper(subtypes))
#table(str_to_upper(names(pancreasMarkers)))
count_data = count_data[ which(str_to_upper(subtypes) %in% str_to_upper(names(pancreasMarkers))), ]
subtypes = subtypes[ str_to_upper(subtypes) %in% str_to_upper(names(pancreasMarkers)) ]
dim(count_data)
length(subtypes)

eislet = new("ExpressionSet", exprs = 
               t(as.matrix(count_data))
)
subtypes[subtypes == "alpha"] = "Alpha"
subtypes[subtypes == "beta"] = "Beta"
subtypes[subtypes == "gamma"] = "Gamma"
subtypes[subtypes == "delta"] = "Delta"
subtypes[subtypes == "epsilon"] = "Epsilon"
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

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0
res_coeff[ res_cor[,"Correlation"] <= .3, ] = .0
res_cor[ res_cor[,"Correlation"] <= .3, "Correlation" ] = 0.0


###

meta_match = match( colnames(exprs(eset)), meta_info$Name, nomatch = 0 ) 

meta_info$Deco_cor[meta_match]  = as.character(res_cor[,2])
res_coeff_match = match( rownames(res_coeff), meta_info$Name, nomatch = 0 )
res_cor_match = match( rownames(res_cor), meta_info$Name, nomatch = 0 )

maxi = apply( res_coeff, FUN = which.max, MARGIN = 1 )
cell_type = colnames(res_coeff)[maxi]
rownames(meta_info) = meta_info$Name
cell_type[res_cor[,"Correlation"] <= .3] = "Not_sig"
meta_info$Deco_type[meta_match] = cell_type
meta_data = meta_info[meta_match,]

#meta_data = data.frame(
#  "Sample" = colnames(exprs(eset)),
#  "Deco_type" = cell_type
#)
#rownames(meta_data) = meta_data$Sample

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv",sep="\t", row.names = F, quote = F)

pheatmap::pheatmap(
  t(res_coeff),
  annotation_col = meta_data[c("Deco_type")],
  annotation_colors = aka3,
  annotation_legend = T,
  treeheight_col = 0,
  show_colnames = T,
  show_rownames = T,
  gaps_row = 20
)

#write.table(meta_info, "~/Deko/Misc/Deko_table.tsv",sep ="\t", row.names = F, quote = F)