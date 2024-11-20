# autoMrP 1.1.1

+ implements knn as an additional classifier
+ best.subset.deep and pca.deep now implements deep interactions for best subset and pca.

# autoMrP 1.1.0

+ implements stacking

# autoMrP 1.0.6

+ implements Deep MrP by Gopelrud as presented in https://doi.org/10.1017/S0003055423000035
+ Set argument deep.mrp = TRUE to include Deep MrP in the ensemble 

# autoMrP 1.0.5

+ drops missing values on y, L1.x, L2.x, L2.unit, L2.reg. Missing values on the DV would previously lead to errors in SVM
+ works with continuous DV.

# autoMrP 0.93

+ block sampling in bootstrapping instead of state-stratified sampling

# autoMrP 0.91

+ bootstrapping returns GB prediction
+ predictions do not fail if census data contains more factor levels than training data for SVM and Lasso
+ svm post-stratification uses the user-specified formula instead of all information
+ lasso post-stratification uses correct user-specified context level variables if L2.x and lasso.L2.x differ
+ parallel processing loops are replicable now
