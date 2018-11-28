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