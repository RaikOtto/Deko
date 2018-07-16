library("bseqsc")
library("stringr")
library("grid")     ## Need to attach (and not just load) grid package
#library("pheatmap")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

###
# "~/Deko/Data/baron_data.tsv"
bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T)
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "^X","")
colnames(bam_data) = str_replace(colnames(bam_data),pattern = "\\.","_")

condition = (meta_info$Included == "Yes") & (meta_info$NEC_NET %in% c("NEC","NET")) #& (meta_info$Grading %in% c("G3"))
include_list  = meta_info$Name[ condition ]
exclude_list  = meta_info$Name[ ! condition ]

bam_data = bam_data[, ( colnames(bam_data) %in% include_list )]
dim(bam_data)
bam_data[1:5,1:5]

count_data = read.table("~/Deko/Data/baron_data.tsv",sep ="\t", header = T, stringsAsFactors = F)
subtypes = read.table(  "~/Deko/Data/baron_subtypes.tsv",sep ="\t", header = T, stringsAsFactors = F)[,1]
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
dim(count_data)


dds_raw = DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_data,0),
  colData = meta_data,
  design = ~ NEC_NET,
  tidy = F
)
dds_raw = DESeq2::estimateSizeFactors(dds_raw)

dds_n = DESeq2::varianceStabilizingTransformation(dds_raw)
library("DESeq2")
count_data = assay( dds_n )

vsn::meanSdPlot( count_data)

###