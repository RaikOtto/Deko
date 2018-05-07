
### contamination SM Figures

library("estimate")
librarY(stringr)

ncol = ncol(expr_raw)
nrow = nrow(expr_raw)
rownames = rownames(expr_raw)
colnames = colnames(expr_raw)
expr_raw = as.double(as.character(as.matrix(expr_raw)))
expr_raw = round(expr_raw,1)
expr_raw[expr_raw < 1 ] = 1
expr_raw = matrix(expr_raw, ncol = ncol, nrow = nrow)
expr_raw = cbind( rownames,rownames, expr_raw )
colnames = c("NAME","Description", colnames)
colnames(expr_raw) = colnames

expr_raw[1:5,1:5]
#write.table(expr_raw,"~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.estimate.gct",sep="\t",quote = F, row.names = F)

setwd("/home/ottoraik/R/x86_64-pc-linux-gnu-library/3.4/estimate/extdata/")
estimateScore( 
  "~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.estimate.gct",
  "~/MAPTor_NET/BAMs/Kallisto_three_groups/estimate.tsv"
)


