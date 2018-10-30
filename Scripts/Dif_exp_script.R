library("stringr")
library("limma")

meta_data = meta_info[colnames(bam_data),]
meta_data[1:5,1:5]

groups = meta_data$Subtype
groups[groups == "delta"] = "CASE"
groups[groups != "CASE"] = "CTRL"
design <- model.matrix(~0 + groups)
colnames(design) = c("delta", "Not_delta")

vfit <- lmFit(bam_data,design)
contr.matrix = makeContrasts( contrast = delta - Not_delta ,  levels = design )
vfit <- contrasts.fit( vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))

result_t = topTable( efit, coef = "contrast", number  = nrow(expr_raw), adjust  ="none", p.value = 1, lfc = 0)
result_t$hgnc_symbol = rownames(result_t)
colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")

result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
result_t = result_t[order(result_t$P_value, decreasing = F),]
result_t$Log_FC = round(result_t$Log_FC, 1)
result_t$Average_Expr = round(result_t$Average_Expr, 1)
result_t = result_t[order(result_t$Log_FC,decreasing = T),]

result_t$ABS_FC = abs(result_t$Log_FC)
result_t = result_t[order(result_t$ABS_FC, decreasing = T),]

#write.table("~/Deko//Results/Dif_Exp/Dif_exp_delta_vs_Not_delta_Merge_mat_HSC_Stanescu_Segerstolpe.tsv", x = result_t, sep = "\t", quote = F, row.names = F)

