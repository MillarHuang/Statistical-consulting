library(dplyr)
library(Seurat)
library(patchwork)
require(caret)
#load data
data <- readRDS("scRNAseq_SCT_final_annotated (1).rds")
#Preprocess data
#############################################################################
#1/Select the 3000 representative features that exhibit high cell-to-cell variation 
pbmc <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
# Visualization of the 3000 selected genes
top10 <- head(VariableFeatures(pbmc), 10)#mark the top 10 highly expressed genes
top10
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#2/Perform PCA on the 3000 selected features(genes) 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) #uses the 3000 selected features(genes) to run PCA (has 50 PCs)

#3/Construct the plot of variance proportion explained by principal components
pc=1:50
std_pc = pbmc@reductions$pca@stdev#The standard deviations explained of each PC
var_ratio = std_pc^2/sum(std_pc^2)#the variance proportion explained by each PC
cum_var_ratio = cumsum(var_ratio)#cumulative variance proportion explained:0.6918382 for 10 PCs
#The plot of Variance proportion explained by each principal component
plot(y = var_ratio,x= pc,ylab='Variance proportion explained',xlab='PC number',main='Variance proportion explained',type='l')
abline(v=10,col='red')
#The plot of Cumulative Variance proportion explained
plot(y = cum_var_ratio,x= pc,ylab='Variance proportion explained',xlab='Number of PC',main='Cumulative Variance proportion explained',type='l')
abline(v=10,add=T,col='red')

#4/Perform Umap clustering on top 10 principal components
pbmc <- RunUMAP(pbmc, dims = 1:10,reduction = 'pca')#use PC as input, dims means how many PCs as input features
DimPlot(data, reduction = "umap",label=TRUE)#Dimension reduction plot on a 2D scatter plot(reduced to 2 dimension in umap, each point is a cell)

#save the preprocess data
saveRDS(pbmc, file = "data_preprocess.rds")

#Modeling part
##############################################################
#1/GLMM model with 10 PCs
library(lme4)
pc10 = pbmc@reductions$pca@cell.embeddings[,1:10] #The 10 PCs values
ct = pbmc@meta.data$celltype.labels #cell labels
brain_reg = pbmc@meta.data$brainregion #brain region
brain_reg = factor(brain_reg)
y = pbmc@meta.data$mousemodel #response
y[y!='wt'] = 1 #mice with deficiency of dystrophin isoforms
y[y=='wt'] =0 #mice without deficiency of dystrophin isoforms
y = as.numeric(y)

#1.1 Fit GLMM with 10 PCs : Use random intercept
model = glmer(y ~ pc10+ct+ (1|brain_reg), family = "binomial",nAGQ = 10, control = glmerControl(optimizer="bobyqa"))
summary(model)
ranef(model) #random effect measure
coef_model = data.frame(coef(model)[[1]]) #cofficients of model

#1.2 Calculate Prediction accuracy on training set: 0.6319091
probabilities <- predict(model)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0) #set classification threshold as 0.5
glmm_acc = mean(predicted.classes == y)
glmm_acc

#1.3 Calculate 10-fold CV accuracy:0.6453727
df = data.frame(cbind(y,pc10,ct,brain_reg))
test_idx <- createFolds(1:nrow(df), k = 10, list = TRUE)
acc = rep(NA,10)
for(i in 1:10){
  test = df[test_idx[[i]],]
  train_ = df[-test_idx[[i]],]
  model_ = glmer(y ~ PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+ PC_7 + PC_8 + PC_9 + PC_10 + ct+ (1|brain_reg), family = "binomial",nAGQ = 10, control = glmerControl(optimizer="bobyqa"),data = train_)
  probabilities <- predict(model_,test)
  predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
  # Model accuracy
  acc[i] = mean(predicted.classes == test$y)
}
mean(acc)

#1.4 Pick the top 10 associated genes for each of top 10 PC
weights_pc = pbmc@reductions$pca@feature.loadings#Extract the weights of each gene on each PC
pc_top10_weights = matrix(NA,nrow = 10, ncol=10) #Weights value for top 10 genes on the 10 PC
pc_top10_genes = data.frame(Top1=c(),Top2=c(),Top3=c(),Top4=c(),Top5=c(),Top6=c(),Top7=c(),Top8=c(),Top9=c(),Top10=c()) #top 5 genes' names
for(i in 1:10){
  pc_top10 = weights_pc[,i][order(abs(weights_pc[,i]),decreasing = TRUE)][1:10]
  pc_top10_weights[i,] = pc_top10
  top10_genes = names(pc_top10)
  pc_top10_genes = rbind(pc_top10_genes,data.frame(Top1=c(top10_genes[1]),Top2=c(top10_genes[2]),Top3=c(top10_genes[3]),Top4=c(top10_genes[4]),Top5=c(top10_genes[5]),
                                                   Top6=c(top10_genes[6]),Top7=c(top10_genes[7]),Top8=c(top10_genes[8]),Top9=c(top10_genes[9]),Top10=c(top10_genes[10])))
}
row.names(pc_top10_genes) = paste('PC',1:10,sep='')
pc_top10_genes

#1.5 Pick the top 10 associated genes for the 6 significant PCs in model
#The counts for the top 10 associated genes in the 6 significant PCs
pc_sig = pc_top10_genes[c(1,2,3,6,7,10),]
sort(table(unlist(pc_sig)),decreasing = TRUE)
#Visualization of the 6 significant principal component with their top 10 associated genes
VizDimLoadings(pbmc, dims = c(1,2,3,6,7,10), reduction = "pca",nfeatures = 10)

###################################################################
#2/ GLMM combined with Lasso regression (use lambda = 0.001073725): give best prediction accuracy
library(glmnet)
#sct data for the 3000 selected genes
sct3000 = pbmc@assays$SCT[VariableFeatures(pbmc)]
sct3000 = t(sct3000)
#encode categorical data
encode = function(x){
  x_ = unique(x)
  x_matrix = matrix(0,nrow = length(x),ncol=length(x_))
  for(i in 1:length(x_)){
    x_matrix[,i][x==x_[i]] = 1
  }
  colnames(x_matrix) = x_
  return(x_matrix)
}
#encode the cell type matrix(0 or 1 for each cell type)
ct_matrix = encode(ct)
#encode the brain region matrix
brain_matrix = encode(brain_reg)
#combine all data as train_data
x_train = cbind(sct3000,ct_matrix[,-1],brain_matrix[,-1]) #remove out one dummy variable in ct and brain_region respectively to avoid perfect collinearity

#2.1 Do 10-fold CV on training set to find the lambda that gives best prediction accuracy
cv.lasso <- cv.glmnet(x_train, y, alpha = 1, family = "binomial",nfolds=10) #alpha = 1 means just using L1 penalty(which is Lasso)
#2.2  Fit the final model on the training data
model_lasso <- glmnet(x_train, y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.min) #cv.lasso$lambda.min == 0.001073725

#2.3 The selected features with their coefficients: 1771 features selected
gene_lasso_df = data.frame(Variables = row.names(coef(model_lasso)), Coefficients = coef(model_lasso))
colnames(gene_lasso_df)=c('Variables','Coefficients')
gene_lasso_df[order(gene_lasso_df$Coefficients,decreasing = TRUE),][1:10,] #Top 10 associated features among 1771 features

#2.4 Prediction accuracy on training set:0.9019623
probabilities <- predict(model_lasso,newx = x_train)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
lasso_cv_acc = mean(predicted.classes == y)
lasso_cv_acc
##################################################################
#3/ Check the relationship between the lambda and accuracy & lambda and number of non-zero features: find the elbow point at lambda = 0.015
lam = seq(0.001,0.035,0.001)
num_lasso = rep(NA,length(lam))
acc_lasso = rep(NA,length(lam))
for(i in 1:length(lam)){
  model_lasso <- glmnet(x_train, y, alpha = 1, family = "binomial",
                        lambda = lam[i]) 
  gene_lasso_df_ = data.frame(Variables = row.names(coef(model_lasso)), Coefficients = coef(model_lasso))
  colnames(gene_lasso_df_)=c('Variables','Coefficients')
  gene_lasso_df_ = gene_lasso_df_[gene_lasso_df_$Coefficients>0,] #remove out the features that has 0 coefficient
  top_features_lasso = gene_lasso_df_[order(gene_lasso_df_$Coefficients,decreasing = TRUE),] #Most associated genes
  num_lasso[i] = nrow(top_features_lasso)#number of features that have non-zero coefficient
  #Prediction
  probabilities <- predict(model_lasso,newx = x_train)
  predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
  # Model accuracy
  acc_lasso[i] = mean(predicted.classes == y)
}
par(mfrow=c(1,2))
plot(acc_lasso~lam,type='l', xlab = 'Lambda',ylab = 'Prediction accuracy')
plot(num_lasso~lam,type='l', xlab = 'Lambda',ylab = 'Number of non-zero features')
abline(v=0.015,col='red') #elbow point at lambda = 0.015

#######################################################################
#4/ GLMM combined with Lasso regression (use lambda = 0.035) : Lasso regression selection validation
#4.1 Do Lasso regression with lambda=0.035: 10 non-zero features selected
model_lasso <- glmnet(x_train, y, alpha = 1, family = "binomial",
                      lambda = 0.035) 
#4.2 The 10 selected genes with their coefficients
gene_lasso_df_ = data.frame(Variables = row.names(coef(model_lasso)), Coefficients = coef(model_lasso))
colnames(gene_lasso_df_)=c('Variables','Coefficients')
gene_lasso_df_ = gene_lasso_df_[gene_lasso_df_$Coefficients>0 & gene_lasso_df_$Variables != '(Intercept)',] #remove out the features that has 0 coefficient
top10_features_lasso = gene_lasso_df_[order(gene_lasso_df_$Coefficients,decreasing = TRUE),] #Top 10 associated genes
nrow(top10_features_lasso) #number of non-zero features in the model
top10_features_lasso

#4.3 Fit the 10 non-zero features to GLMM: all 10 features significant
top10_lasso = x_train[, top10_features_lasso$Variables]#remove out the intercept term(8)
top10_lasso = as.matrix(top10_lasso)
model_lasso10 = glmer(y ~ top10_lasso+ct+ (1|brain_reg), family = "binomial",nAGQ = 10, control = glmerControl(optimizer="bobyqa"))
summary(model_lasso10)# Model summary

#4.4 10-fold CV accuracy :0.6721114
df = data.frame(cbind(y,top10_lasso,ct,brain_reg))
test_idx <- createFolds(1:nrow(df), k = 10, list = TRUE)
acc = rep(NA,10)
for(i in 1:10){
  test = df[test_idx[[i]],]
  train_ = df[-test_idx[[i]],]
  model_ = glmer(y ~ Gm10076+Gm42418+Rorb+Slc7a5+Ccnd2+Lars2+ Itgb1 + Uba52 + Hist2h2aa1 + Ttc3 + ct+ (1|brain_reg), family = "binomial",nAGQ = 10, control = glmerControl(optimizer="bobyqa"),data = train_)
  probabilities <- predict(model_,test)
  predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
  # Model accuracy
  acc[i] = mean(predicted.classes == test$y)
}
mean(acc)#0.6721114

#################################################################
#5/ GLMM combined with Lasso regression (use lambda = 0.015) :elbow point
#5.1 Do Lasso regression with lambda=0.015: 55 non-zero features selected
model_lasso <- glmnet(x_train, y, alpha = 1, family = "binomial",
                      lambda = 0.015) 
#5.2 The 55 selected features with their coefficients
gene_lasso_df_ = data.frame(Variables = row.names(coef(model_lasso)), Coefficients = coef(model_lasso))
colnames(gene_lasso_df_)=c('Variables','Coefficients')
gene_lasso_df_ = gene_lasso_df_[gene_lasso_df_$Coefficients>0 & gene_lasso_df_$Variables != '(Intercept)',] #remove out the features that has 0 coefficient and intercerpt term
top10_features_lasso = gene_lasso_df_[order(gene_lasso_df_$Coefficients,decreasing = TRUE),] #Top 10 associated genes
top10_features_lasso #55 non-zero features
top10_lasso = x_train[, top10_features_lasso$Variables]#pick the non-zero features
top10_lasso = as.matrix(top10_lasso)
#5.3 Fit 55 selected features in GLMM: 50 features are significant
model_lasso10 = glmer(y ~ top10_lasso+ct+ (1|brain_reg), family = "binomial",nAGQ = 10, control = glmerControl(optimizer="bobyqa"))
summary(model_lasso10)
#5.4 Prediction accuracy:0.7367938
probabilities <- predict(model_lasso10)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
glmm_acc = mean(predicted.classes == y)
glmm_acc
ranef(model_lasso10) #random effect
#save GLMM model results summary
sink('output_top10.txt')
print(summary(model_lasso10))
sink()



