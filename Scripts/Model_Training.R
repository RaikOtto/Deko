subtypes = meta_info[colnames(count_data),"Subtype"]
names(subtypes) = meta_info[colnames(count_data),"Name"]
table(subtypes)

marker_genes = read.table(
    "~/Deko/Misc/Baron_pancreas_marker.tsv",
    sep = "\t",
    header = T,
    stringsAsFactors = F
)
delimiter = 1:100
marker_genes = marker_genes[delimiter,]

pancreasMarkers = list(
    "Alpha" = marker_genes$Alpha[ (marker_genes$Alpha %in% rownames(count_data))],
    "Beta" = marker_genes$Beta[ (marker_genes$Beta != "") & (marker_genes$Beta %in% rownames(count_data))],
    "Gamma" = marker_genes$Gamma[ (marker_genes$Gamma != "")& (marker_genes$Gamma %in% rownames(count_data))],
    "Delta" = marker_genes$Delta[ (marker_genes$Delta != "")& (marker_genes$Delta %in% rownames(count_data))],
    #"Ductal" = marker_genes$Ductal[ (marker_genes$Ductal != "")& (marker_genes$Ductal %in% rownames(count_data))],
    #"Acinar" = marker_genes$Acinar[ (marker_genes$Acinar != "")& (marker_genes$Acinar %in% rownames(count_data))]#,
    "Progenitor" = marker_genes$Progenitor[ (marker_genes$Progenitor != "")& (marker_genes$Progenitor %in% rownames(count_data))],
    "HISC" = marker_genes$HESC[(marker_genes$HISC != "")& (marker_genes$HISC %in% rownames(count_data))]
    #E17.5 = names(E17.5),
)

cands = names(subtypes)[ str_to_upper( subtypes ) %in% str_to_upper( names(pancreasMarkers))]
length(cands)
count_data = count_data[, cands]
subtypes = subtypes[cands]
dim(count_data)
table(subtypes)

### normalization

row_var = apply(count_data, FUN = var, MARGIN = 1)
col_var = apply(count_data, FUN = var, MARGIN = 2)
table(row_var == 0)
table(col_var == 0)
count_data = count_data[row_var != 0,col_var != 0]
count_data = count_data[rowSums(count_data) >= 1,]
dim(count_data)

table(as.character(unlist(pancreasMarkers)) %in% rownames(count_data))
for(marker in names(pancreasMarkers)){
    genes = as.character(unlist(pancreasMarkers[marker]))
    genes = genes[genes %in% rownames(count_data)]
    pancreasMarkers[marker] = list(genes)
}

names(pancreasMarkers) = str_to_lower(names(pancreasMarkers))
eislet = new("ExpressionSet", exprs = as.matrix(count_data))

sub_list = str_to_lower(subtypes)
names(sub_list) = names(subtypes)
fData(eislet) = data.frame( sub_list  )
pData(eislet) = data.frame( sub_list )
names(pancreasMarkers)[!(names(pancreasMarkers) %in% sub_list)]

B = bseqsc_basis(
    eislet,
    pancreasMarkers,
    clusters = 'sub_list',
    samples = colnames(exprs(eislet)),
    ct.scale = FALSE
)
plotBasis(B, pancreasMarkers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

#saveRDS(B,"~/ArtDeco/inst/Models/Four_differentiation_stages_Segerstolpe_Progenitor_HISC.RDS")
