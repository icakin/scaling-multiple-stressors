# growth curve models taken from nlsMicrobio

baranyi_without_lag <- function(log10_nmax, log10_n0, mumax, t){
  return(log10_nmax - log10(1 + (10^(log10_nmax - log10_n0) - 1) * exp(-mumax * t)))
}
baranyi <- function(log10_nmax, log10_n0, mumax, t, lag){
  return(log10_nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(log10_nmax - log10_n0))))
}
baranyi_without_nmax <- function(log10_n0, mumax, t, lag){
  return(log10_n0 + mumax * t/log(10) + log10(exp(-mumax * t) * (1 - exp(-mumax * lag)) + exp(-mumax * lag)))
}
gompertz <- function(log10_nmax, log10_n0, mumax, t, lag){
  log10_n0 + (log10_nmax - log10_n0) * exp(-exp(mumax * exp(1) * (lag - t)/((log10_nmax - log10_n0) * log(10)) + 1))
}
buchanan <- function(log10_nmax, log10_n0, mumax, t, lag){
  log10_n0 + (t >= lag) * (t <= (lag + (log10_nmax - log10_n0) * log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (log10_nmax - log10_n0) * log(10)/mumax)) * (log10_nmax - log10_n0)
}
buchanan_without_lag <- function(log10_nmax, log10_n0, mumax, t){
  log10_n0 + (t <= ((log10_nmax - log10_n0) * log(10)/mumax)) * mumax * t/log(10) + (t > ((log10_nmax - log10_n0) * log(10)/mumax)) *  (log10_nmax - log10_n0)
}

# calculate Topt
get_topt <- function(Eh, Th, Ea){
  return((Eh * Th)/(Eh + (8.62e-05 * Th * log((Eh/Ea) - 1))))
} 

# model for sharpe-schoolfield model
sharpeschoolhigh_1981 <- function(temp_k, rtref, e, eh, th, tref){
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- rtref*exp(e/k * (1/tref - 1/temp_k))
  inactivation.term <- 1/(1 + exp(eh/k * (1/th - 1/temp_k)))
  return(boltzmann.term * inactivation.term)
}