fit_gam_model <- function(train_data, test_data, kval, outcome, predictor) {
  
  formula_string <- paste0("log(", outcome, ") ~ s(", predictor, ", bs='cc', k = ", kval, ")")
  formula <- as.formula(formula_string)
  
  gam_model <- gam(formula, data=train_data) # fit
  predictions <- predict(gam_model, newdata = test_data) # predict on test data
  rmse <- sqrt(mean((log(test_data[,outcome] %>% c() %>% unlist) - predictions)^2)) # performance metric
  return(rmse)
} 

mygroupKFold = function(x){
  myfold = list()
  for(kf in 1:length(unique(x))){
    myfold[[kf]] = which(x!=unique(x)[kf])
  }
  return(myfold)
}

fit_gam_model_kout <- function(data, outcome, predictor) {
  # tmp = data %>% dplyr::select(all_of(outcome)) %>% c() %>% unlist
  # folds <- createFolds(tmp, k = 326, list = TRUE) # from the caret package, this k is in the leave-k-out context, not the basis dimension of the spline
  # Also, if k=length of data, this is equivalent to leave-one-out validation; computational time will suffer
  folds <- mygroupKFold(data$month)
  
  k_values <- 3:11 # k values to test (basis dimension)
  
  results <-  rep(NA, max(k_values)) # Initialize results vector
  
  # Perform cross-validation for each k value
  for (kvali in k_values) {
    fold_errors <- sapply(folds, function(fold) {
      train_data <- data[fold, ]
      test_data <- data[-fold, ]
      rmse <- fit_gam_model(train_data, test_data, kvali, outcome, predictor)
      return(rmse)
    })
    mean_rmse <- mean(fold_errors)
    results[[kvali]] <-  mean(fold_errors)
  }
  
  return(results)
}

# Random iterations out of a spline function
gam_predictions = function(model, newpred, iterates){
  tmp = predict.gam(model, newdata=newpred, type="lpmatrix") # se.fit ignored when lpmatrix is the type of prediction
  somebetas = mgcv::rmvn(n=iterates, coef(model), vcov(model))
  # mvtnorm::rmvnorm(n=iterates, coef(model), vcov(model)) #se.fit taken care of here.
  someiterates = (tmp %*% t(somebetas))
  allpred = someiterates
  return(allpred)
  
  # Quick and dirty plots to check the outputs
  # matplot(allpred, type="l", col=rgb(0,0,0,alpha=.1), lty=1)
  # with(data, points(newpred[,1], model[["y"]], pch=20, col="red"))
  
  # Any back-transformation (exponentiation, etc) to correct for transforming y must happen outside of this function
}