GMT <- function (titers) {
  exp(mean(log(titers)))
}

GMT_CI <- function (titers) {
  lt <- log(titers)
  sdlt <- sd(lt) # should this be SE?
  mlt <- mean(lt)
  lci <- exp(mlt - 1.96*(sdlt))
  uci <- exp(mlt + 1.96*(sdlt))
  return(c(lci, uci))
}

GMT_CI_low = function (titers) {GMT_CI(titers)[[1]]}
GMT_CI_high = function (titers) {GMT_CI(titers)[2]}

GMT_CI_se <- function (titers) {
  lt <- log(titers)
  sdlt <- sd(lt)/sqrt(length(titers))
  mlt <- mean(lt)
  lci <- exp(mlt - 1.96*(sdlt))
  uci <- exp(mlt + 1.96*(sdlt))
  return(c(lci, uci))
}

GMT_CI_se_low = function (titers) {GMT_CI_se(titers)[[1]]}
GMT_CI_se_high = function (titers) {GMT_CI_se(titers)[2]}

# standard error function
GMT_SE <- function (titers) {
  lt <- log(titers)
  lse <- sd(lt)/sqrt(length(titers))
  exp(lse)
}