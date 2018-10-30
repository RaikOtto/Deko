library("bseqsc")
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
#library("pheatmap")

#write.table(meta_info, "~/Deko/Misc/Meta_information.tsv", row.names = F, quote =F , sep = "\t")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###

#bam_data = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Groetzinger_57.tsv",sep ="\t", header = T)
bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T, row.names = 1)
rownames(bam_data) = str_to_upper( rownames( bam_data) )

#bam_data = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Groetzinger_57.tsv",sep ="\t", header = T)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")
hgnc_list = rownames(bam_data)
hgnc_list_uni = unique(hgnc_list)

source("~/Deko/Scripts/Variance_selection.R")

bam_data[1:5,1:5]
dim(bam_data)

#condition = (meta_info$Included == "Yes") & (meta_info$NEC_NET %in% c("NEC","NET")) #& (meta_info$Grading %in% c("G3"))
#include_list  = meta_info$Name[ condition ]
#exclude_list  = meta_info$Name[ ! condition ]

#bam_data = bam_data[, ( colnames(bam_data) %in% include_list )]
#dim(bam_data)
#bam_data[1:5,1:5]

#count_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Baron_human.tsv",sep ="\t", header = T, stringsAsFactors = F)
count_data = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)
#count_data = read.table("~/Deko/Data/Merge_mat_HSC_Stanescu_Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)
colnames(count_data) = str_replace(colnames(count_data), pattern = "\\.", "_")
count_data[1:5,1:5]
subtypes = meta_info[colnames(count_data),"Subtype"]
names(subtypes) = meta_info[colnames(count_data),"Name"]

marker_genes = read.table(
    "~/Deko/Misc/Baron_pancreas_marker.tsv",
    sep = "\t",
    header = T,
    stringsAsFactors = F
)
delimiter = 1:100

pancreasMarkers = list(
    "Alpha" = marker_genes$Alpha[delimiter],
    #"Alpha_five" = marker_genes$Alpha[delimiter],
    "Beta" = marker_genes$Beta[delimiter & (marker_genes$Beta != "")],
    "Gamma" = marker_genes$Gamma[delimiter & (marker_genes$Gamma != "")],
    "Delta" = marker_genes$Delta[delimiter & (marker_genes$Delta != "")]#,
    #"Ductal" = marker_genes$Ductal[delimiter],
    #"Acinar" = marker_genes$Acinar[delimiter]#,
   # "E13.5" = marker_genes$E13.5[delimiter & (marker_genes$E13.5 != "")],
#    "HSC" = marker_genes$HSC[delimiter & (marker_genes$HSC != "")]
    #E17.5 = names(E17.5),
    #P0 = names(P0),
    #P3 = names(P3),
    #P9 = names(P9),
    #P15 = names(P15),
    #P18 = names(P18),
    #P60 = names(P60),
    #Pancreas = names(Pancreas)
)

cands = names(subtypes)[ str_to_upper( subtypes ) %in% str_to_upper( names(pancreasMarkers))]
count_data = count_data[, cands]
dim(count_data)
count_data[1:5,1:5]
sub_list = as.character( subtypes[ cands] )
table( sub_list )

# TPM count transformation

gene_length_t = read.table("~/Deko/Misc/gene_length.tsv",sep ="\t", header = T, stringsAsFactors = F)
l_match = match(rownames(count_data), gene_length_t$hgnc_symbol, nomatch = 0)
count_data = count_data[ l_match != 0,]
l_match = match( gene_length_t$hgnc_symbol, rownames(count_data), nomatch = 0)
gene_length = gene_length_t$transcript_length[l_match]
x = count_data / gene_length
count_data = t(x) * 1e6 / colSums(x)
count_data = t(count_data)

#seg_meta = read.table("~/Deko/Misc/Segerstolpe_Meta_info.tsv", sep ="\t", header = T)

#count_data = t(count_data)

### normalization

row_var = apply(count_data, FUN = var, MARGIN = 1)
col_var = apply(count_data, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
count_data = count_data[row_var != 0,col_var != 0]
count_data = count_data[rowSums(count_data) >= 1,]
#count_data = t(count_data)
dim(count_data)

### load_data

#for( type in names(pancreasMarkers)){
#    cell_type = (eval(paste(type)))
#    genes = as.character(unlist(pancreasMarkers[cell_type]))
#    genes = genes[genes %in% rownames(count_data)]
#    genes = rownames(count_data)
#    pancreasMarkers[cell_type] = list(genes)
#}

table(as.character(unlist(pancreasMarkers)) %in% rownames(count_data))
table( !( as.character(unlist(pancreasMarkers)) %in% rownames(count_data) ))

#write.table(pancreasMarkers,"~/Deko/Data/Mouse_progenitor_pancreas_scRNA/Pancreas_marker_Zang.tsv",sep = "\t", quote = F)

