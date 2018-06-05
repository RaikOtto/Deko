library("stringr")
library("caret")

r_mat = read.table("~/Koop_Klinghammer/Results/Linear_regression/r_mat.tsv",sep="\t", header = T)
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv", sep ="\t", stringsAsFactors = F, header = T)
meta_info = subset( meta_info, Included == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
meta_info$OS = as.double(meta_info$OS)

f_mat = meta_info[ ,c("Subtype","Geschlecht","OS","Loc_Primaertumor","Treatment","Chemozyklen","Best_response","Vortherapie","Raucher","Alkohol","Anzahl_py")]
f_mat$OS = str_replace(f_mat$OS, pattern = ",",".")
f_mat$OS = as.double(f_mat$OS)

###

s_match = match(rownames(t(r_mat)), meta_info$Name, nomatch = 0)
r_mat = r_mat[s_match != 0,]
meta_data = meta_info[ s_match, ]

r_mat = cbind( as.character( meta_data$Best_response ) , t(pure_data))
colnames(r_mat)[1] = "Best_response"
#r_mat = cbind( as.character( rownames(r_mat) ) , r_mat)
#colnames(r_mat)[1] = "Sample"
#lambda <- 10^seq(10, -2, length = 100)

splitIndex = caret::createDataPartition( 1:nrow(r_mat), p = .75, list = FALSE, times = 1)
trainDF = r_mat[ splitIndex,]
testDF  = r_mat[-splitIndex,]

predictorsNames = colnames(r_mat)[2:(ncol(r_mat))]

tune_grid = tuneGrid = expand.grid(
  .alpha=1,
  .lambda=seq(0, 100, by = 0.1)
)

objControl <- caret::trainControl(method='cv', number=3, returnResamp='none', classProbs = TRUE)

trainControl = caret::trainControl(
  method = "cv",
  number = 10#,
  #summaryFunction = prSummary,
  #classProbs = T
)

prop.table(table(r_mat[,1]))

modelFit = caret::train(
  trainDF[,predictorsNames],
  trainDF[,"Best_response"],
  method='gbm', 
  trControl=objControl,  
  metric = "ROC",
  preProc = c("center", "scale")
)

predictions <- predict(object = modelFit, testDF[,predictorsNames], type='raw')
table(as.character(predictions), testDF[,1])
#auc <- pROC::roc( testDF[,outcomeName], predictions)
#print(auc$auc)

d=varImp(modelFit,scale=T)
imp_vars = d$importance
vis_vars = imp_vars[order(imp_vars[,1], decreasing = T),]
names(vis_vars) = rownames(imp_vars)[order(imp_vars[,1], decreasing = T)]

plot( vis_vars[1:10] )

aggregate( pure_data[rownames(pure_data) == "SPRR3",], by = list(meta_data$Subtype), FUN = mean)
###

#create test and training sets

vimp <- varImp(modelFit, scale=F)
results <- data.frame(row.names(vimp$importance),vimp$importance$Overall)
results$VariableName <- rownames(vimp)
colnames(results) <- c('VariableName','Weight')
results <- results[order(results$Weight),]
results <- results[(results$Weight != 0),]

par(mar=c(5,15,4,2)) # increase y-axis margin. 
xx <- barplot(results$Weight, width = 0.85, 
              main = paste("Variable Importance -",outcomeName), horiz = T, 
              xlab = "< (-) importance >  < neutral >  < importance (+) >", axes = FALSE, 
              col = ifelse((results$Weight > 0), 'blue', 'red')) 
axis(2, at=xx, labels=results$VariableName, tick=FALSE, las=2, line=-0.3, cex.axis=0.6)  

###

write.table(r_mat, "~/Koop_Klinghammer/Results/Linear_regression/r_mat.tsv",quote = F, row.names = F, sep ="\t")
