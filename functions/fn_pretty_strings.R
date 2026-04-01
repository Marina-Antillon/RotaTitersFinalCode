## functions to make confidence interval summaries presentable
# Express the number in pretty print as <0.001 if it's too small
# useful for p-values and such
# Only feed single numbers to this function. ifelse is picky about that.
# Within a dplyr mutate function, this function should work fine.
my_smallvalue = function(x, dec){
  ifelse(x<10^-dec, paste("<", 10^-dec, sep=""), format(round(x, dec), nsmall=dec))
}

ci_string_dec = function(myvector, dec){
  result = paste(my_smallvalue(sort(myvector)[2], dec), " (", 
                 my_smallvalue(min(myvector), dec), ", ", 
                 my_smallvalue(max(myvector), dec), ")", sep="")
  return(result)
}

time_range_concat = function(x){
  # Check that the input is a numeric vector of length 3
  if (length(x) != 3) {stop("Input must be a vector of exactly 3 elements.")}
  
  # Generate the formatted string
  result <- paste0(x[1], " (", x[2], ", ", x[3], ")")
  
  # Issue a warning about vector order responsibility
  warning("Ensure the input vector is in the correct order: (lower CI, estimate, upper CI).")
  
  return(result)
}

# To remove leading white space in formatted numbers. Occasionally this is an issue.
my_format = function(x, dec){trimws(format(round(x, dec), nsmall=dec))}

# not quite sig-figs, but close
my_format_mod = function(x,dec,bigx){
  ifelse(bigx<=x, round(x), my_format(x, dec))
}

range_pretty = function(x){
  return(paste0("[", round(x[1],0), "-", round(x[2],0), "]"))
}
