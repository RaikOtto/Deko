library("org.Hs.eg.db")
### hgnc

expr_raw = bam_data

gene_ids = row.names(expr_raw)

hgnc_list = as.character( mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first"
) )



expr_raw = expr_raw[ ! is.na(hgnc_list),]
hgnc_list = hgnc_list[!is.na(hgnc_list)]
hgnc_list_uni = as.character( unique(hgnc_list) )

max_list = as.integer( sapply( hgnc_list_uni, FUN = function(gene){
  
  var_match = which( hgnc_list == gene )
  if (length(var_match) > 1){
    row_var_max = which.max( apply( as.matrix(expr_raw[var_match,]), FUN = var, MARGIN = 1)  );
  } else{
    row_var_max = 1
  }
  return( var_match[row_var_max] )
}))

expr_raw = expr_raw[max_list,]
rownames(expr_raw) = hgnc_list[max_list]
dim(expr_raw)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "" )
rownames(expr_raw) = str_to_upper(rownames(expr_raw))
expr_raw[1:5,1:5]

bam_data = expr_raw
#write.table(expr_raw,"~/Deko/Data/Mouse_progenitor_pancreas_scRNA/Zhang.tsv",sep="\t",quote =F , row.names = T)
