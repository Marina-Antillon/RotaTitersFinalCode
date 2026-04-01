# Normal Q-Q plots
# in theory, qqplotr has a function to do the same, but when I tried to install it, it had dependencies that somehow would not install themselves...
# https://rdrr.io/cran/qqplotr/man/stat_qq_band.html
gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  
  # Code for this function here: https://gist.github.com/rentrop/d39a8406ad8af2a1066c
  # in theory, qqplotr has a function to do the same, but when I tried to install it, it had dependencies that somehow would not install themselves...
  
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value = coef[1] + coef[2] * df$z
  df$upper = fit.value + zz * SE
  df$lower = fit.value - zz * SE
  jb_lbl = my_smallvalue(jarque.test(x)$p.value, 3)
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) + 
    theme_bw() + 
    xlab("Theoretical quantiles") + ylab("Standardized residuals") # + 
  # annotation_custom(
  #   grob = textGrob(paste0("Jarque-Bera Normality Test\nP-value = ", jb_lbl), 
  #                   x = unit(0.25, "npc"), y = unit(0.15, "npc"), 
  #                   just = "center",
  #                   gp = gpar(col = "black", fontsize = 10)),
  #   xmin = -Inf, xmax = Inf, ymin = Inf, ymax = Inf)
  
  if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  return(p)
  # coef
}