library("stringr")
library("limma")

groups = meta_data$Deco_type
groups[groups != "Not_sig"] = "Sig"
groups[groups != "Sig"] = "Not_sig"
design <- model.matrix(~0 + groups)
colnames(design) = c("Not_Sig", "Sig")

vfit <- lmFit(expr_raw,design)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#plotSA(efit)

contr.matrix = makeContrasts( contrast = Not_Sig - Sig ,  levels = design )
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

write.table("~/Deko//Results/Dif_Exp/Not_Sig_Minus_Sig.tsv", x = result_t, sep = "\t", quote = F, row.names = F)

