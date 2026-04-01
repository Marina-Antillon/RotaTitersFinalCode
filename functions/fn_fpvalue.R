# Getting the p-value from an f
# broom::glance(fit1A) could do the same, if broom is loaded.
fpvalue <- function(model) {
  model_sum = summary(model)
  model_sum_F = model_sum$fstatistic["value"]
  model_sum_numdf = model_sum$fstatistic["numdf"]
  model_sum_dendf = model_sum$fstatistic["dendf"]
  out = list("Fp" = (pf(model_sum_F, model_sum_numdf, model_sum_dendf, lower.tail = F) %>% unname))
  return(out)
}
