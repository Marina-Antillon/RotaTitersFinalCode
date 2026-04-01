extract_Int_CI = function(model, alpha=0.05){
  intercept <- coef(model)[1] # Extract the intercept (the first coefficient)
  se <- sqrt(vcov(model)[1,1]) # Get the standard error for the intercept
  z_value <- qnorm(1 - alpha/2)
  
  ests = c(intercept, intercept - z_value * se, intercept + z_value * se)
  return(ests)
}