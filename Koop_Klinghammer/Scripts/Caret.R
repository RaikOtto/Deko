library("stringr")
library("caret")

r_mat = read.table("~/Koop_Klinghammer/Results/Linear_regression/r_mat.tsv",sep="\t", header = T)
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv", sep ="\t", stringsAsFactors = F, header = T)
meta_info = subset( meta_info, Included == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
colnames(meta_info) = str_replace_all(colnames(meta_info), pattern = "\\.AF8\\.","_")
meta_info$OS = str_replace_all( meta_info$OS, pattern = ",", "." )
meta_info$OS = as.double(meta_info$OS)

f_mat = meta_info[ ,c("Subtype","Geschlecht","OS","Loc_Primaertumor","Treatment","Chemozyklen","Best_response","Vortherapie","Raucher","Alkohol","Anzahl_py")]
f_mat$OS = str_replace(f_mat$OS, pattern = ",",".")
f_mat$OS = as.double(f_mat$OS)

s_match = match(rownames(r_mat), meta_info$Name, nomatch = 0)
r_mat = r_mat[s_match != 0,]
meta_data = meta_info[ s_match, ]

tuneGrid = expand.grid(
  .alpha=1,
  .lambda=seq(0, 100, by = 0.1)
)

trainControl = trainControl(
  method = "cv",
  number = 10#,
  #summaryFunction = prSummary,
  #classProbs = T
)

trainControl = trainControl(method='cv', number=3, returnResamp='none')

r_mat = cbind( as.double( meta_data$OS ) , r_mat)
colnames(r_mat)[1] = "OS"
r_mat = cbind( as.character( rownames(r_mat) ) , r_mat)
colnames(r_mat)[1] = "Sample"
lambda <- 10^seq(10, -2, length = 100)

splitIndex = createDataPartition( 1:nrow(r_mat), p = .75, list = FALSE, times = 1)
trainDF = r_mat[ splitIndex,]
testDF  = r_mat[-splitIndex,]

predictorsNames = colnames(r_mat)[3:ncol(r_mat)]

expand.grid(alpha = seq(0.1, 1, length = len),
            lambda = seq(.1, 3, length = 3 * len))

modelFit = train(
  trainDF[,predictorsNames],
  trainDF[,"OS"],
  method = "glmnet",
  lambda = 0.0,
  #tuneLength = 10,
  #trControl = trainControl,
  metric = "RMSE",
  preProc = c("center", "scale")
  #family="binomial"
)

predictions <- predict(object = modelFit, testDF[,predictorsNames], type='raw')
plot(predictions, testDF$OS)
#auc <- pROC::roc( testDF[,outcomeName], predictions)
#print(auc$auc)

plot(varImp(modelFit,scale=T))
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
