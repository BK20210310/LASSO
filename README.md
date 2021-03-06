# LASSO (Least Absolute Shrinkage and Selection Operator)

LASSO is one method for feature selection. LASSO adds L1 penalty on GLM (Generalized Linear Model) to select important features. 

## Introduction

Here, I used gene expression data as an example to implement LASSO.

GSE130563 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130563) is a microarray data which analyzed 23786 genes of skeletal muscle samples from PC (Pancreatic Cancer) patients with and without cachexia.

## Data visualization
First, that's take a look of 23786 genes (features).

Volcano plot is a data visualization tool for microarray data. Fold change and p-value are used to compare gene expression between cachexia and non-cachexia.
The horizontal axis is log fold change and vertical axis is -log p-value. P-value was calculated by t-test.

![volcano plot](https://user-images.githubusercontent.com/80352910/140676210-5211e0c7-bfb6-40ff-94f0-63abd70c72ed.png)

## Logistic regression

For example, I simply picked 8 independent variables (genes) to run logistic regression and used AIC ( Akaike Information Criterion) to evaluate the fitness of model.

The AIC of logistic regression was 48.765.

![logistic regression](https://user-images.githubusercontent.com/80352910/141036895-da7e82a9-1c0c-433b-b6f3-dee7162aee9a.PNG)

## LASSO

Then added L1 penalty to logistic regression. It can select important independent variables from 8 genes. When λ increased, the coefficients of independent variables would shrinkage to zero.

![LASSO plot](https://user-images.githubusercontent.com/80352910/140690682-c1e770c7-fccc-4428-864f-f443d1e5fe24.png)

In this case, I used 10-fold cross validation to choose optimal λ.

![MAE plot](https://user-images.githubusercontent.com/80352910/140690674-d9e2d699-6f41-40b2-8b83-d0491700fca1.png)

## Post-LASSO

After that, I used 4 independent variables which were selected by LASSO to run logistic regression again.

The AIC of this logistic regression was 41.055. It was smaller than previous AIC.

![post-LASSO](https://user-images.githubusercontent.com/80352910/141037921-ff84e56c-6a0b-47c9-8e97-83e568f53cd2.PNG)

Finally, using 4 independent variables was better than 8 independent variables to construct logistic regression. 
