expr_raw = bam_data

hgnc_list = row.names(expr_raw)
expr_raw = expr_raw[ ! is.na(hgnc_list),]
expr_raw = expr_raw[ hgnc_list != "",]
hgnc_list = hgnc_list[!is.na(hgnc_list)]
hgnc_list = hgnc_list[hgnc_list != ""]

hgnc_list = str_to_upper(hgnc_list)
#expr_raw = expr_raw[ hgnc_list != "",]
#hgnc_list = hgnc_list[ hgnc_list != ""]

#hgnc_list = str_replace_all( hgnc_list,pattern= "(\\.)|(-)|(_)","")

hgnc_list_uni = as.character( unique(hgnc_list) )

max_list <<- as.integer( sapply( hgnc_list_uni, FUN = function(gene){
  
  var_match = which( hgnc_list == gene )
  if (length(var_match) > 1){
        row_var_max = which.max( 
            apply( 
                matrix(
                    as.integer(
                        as.character(
                            unlist(
                                expr_raw[var_match,]
                            )
                        )
                    ),
                    ncol = length(var_match)
                ),
                FUN = var,
                MARGIN = 2
            )
        );
  } else{
    row_var_max = 1
  }
  return( var_match[row_var_max] )
}))

max_list[is.na(max_list)] = 0
expr_raw = expr_raw[max_list,]
rownames(expr_raw) = hgnc_list[as.integer(max_list)]

expr_raw[(nrow(expr_raw)-10):nrow(expr_raw),1:5]
dim(expr_raw)
length(hgnc_list_uni)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "" )

bam_data = expr_raw
