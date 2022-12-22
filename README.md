# Statistical-consulting-code
This R code file contains 2 parts: data preprocessing part and the modeling part. Below is the structure of this file. You can find the corresponding code by this order in this R code file.

## Data preprocessing: 

1.Select the 3000 representative features that exhibit high cell-to-cell variation

2.Perform PCA on the 3000 selected features(genes) 

3.Construct the plot of variance proportion explained by principal components

4.Perform Umap clustering on top 10 principal components

## Modeling part:

#### 1.GLMM model with 10 PCs

1.1 Fit GLMM with 10 PCs : Use random intercept

1.2 Calculate Prediction accuracy on training set

1.3 Calculate 10-fold CV accuracy

1.4 Find the top 10 associated genes for each of the top 10 PCs

1.5 Find the top 10 associated genes for each of the 6 significant PCs in model

#### 2.GLMM combined with Lasso regression (use lambda = 0.001073725): give best the prediction accuracy

2.1 Do 10-fold CV on training set to find the lambda that gives best prediction accuracy (0.001073725)

2.2 Do Lasso regression with lambda = 0.001073725 (1771 features selected)

2.3 The selected features by Lasso regression with their coefficients 

2.4 Prediction accuracy on training set

#### 3.Check the relationship between the lambda and prediction accuracy & lambda and the number of non-zero features (find the elbow point at lambda = 0.015)

#### 4.GLMM combined with Lasso regression (use lambda = 0.035) : Lasso regression selection validation

4.1 Do Lasso regression with lambda=0.035 (10 features selected)

4.2 The 10 features selectedby Lasso regression with their coefficients

4.3 Fit the 10 selected features to GLMM (all 10 features are significant)

4.4 calculate 10-fold CV accuracy

#### 5.GLMM combined with Lasso regression (use lambda = 0.015) :elbow point

5.1 Do Lasso regression with lambda=0.015 (55 non-zero features selected)

5.2 The 55 features selecedt by Lasso regression with their coefficients

5.3 Fit 55 selected features in GLMM (50 features are significant)

5.4 Prediction accuracy








