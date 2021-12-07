# autoMrP 0.93

+ block sampling in bootstrapping instead of state-stratified sampling

# autoMrP 0.91

+ bootstrapping returns GB prediction
+ predictions do not fail if census data contains more factor levels than training data for SVM and Lasso
+ svm post-stratification uses the user-specified formula instead of all information
+ lasso post-stratification uses correct user-specified context level variables if L2.x and lasso.L2.x differ
+ parallel processing loops are replicable now
