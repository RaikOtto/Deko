library("stringr")
library("limma")

identify_marker_genes = function(
    expression_training_mat,
    label_vector,
    label,
    nr_marker_genes
){

    print(paste("Calculating marker genes for subtype: ",label, sep =""))
    
    groups = label_vector
    groups[groups == label] = "CASE"
    groups[groups != "CASE"] = "CTRL"
    design <- model.matrix(~0 + groups)
    colnames(design) = c("Case","Ctrl")
    
    vfit <- lmFit(expression_training_mat,design)
    contr.matrix = makeContrasts( contrast = 
        Case - Ctrl,  levels = design )
    vfit <- contrasts.fit( vfit, contrasts = contr.matrix)
    efit <- eBayes(vfit)
    
    result_t = topTable(
        efit,
        coef = "contrast",
        number  = nrow(expression_training_mat),
        adjust  ="none",
        p.value = 1,
        lfc = 0
    )
    result_t$hgnc_symbol = rownames(result_t)
    colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")
    
    result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
    result_t = result_t[order(result_t$P_value, decreasing = F),]
    result_t$Log_FC = round(result_t$Log_FC, 1)
    result_t$Average_Expr = round(result_t$Average_Expr, 1)
    result_t = result_t[order(result_t$Log_FC,decreasing = T),]
    
    #result_t$ABS_FC = abs(result_t$Log_FC)
    result_t = result_t[order(result_t$Log_FC, decreasing = T),]
    marker_genes = as.character(result_t$HGNC)
    marker_genes = marker_genes[marker_genes %in% rownames(expression_training_mat)]
    marker_genes = marker_genes[1:nr_marker_genes]
    
    #write.table("~/Deko/Results/Dif_Exp/Dif_exp_HISC_vs_Not_HISC_Differentiated_Segerstolpe_Prog_Hisc_scRNA.tsv", x = result_t, sep = "\t", quote = F, row.names = F)
    return(marker_genes)
}