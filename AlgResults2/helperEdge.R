sdLog <- function(Parray) {
  Parray[Parray == 0] <- 1e-11
  Plog <- log(Parray)
  Pm  <- apply(Plog, c(1, 2), mean, na.rm = TRUE)
  Psd <- apply(Plog, c(1, 2), sd,   na.rm = TRUE)
  return(list(Pmean = exp(Pm), 
              Plow = exp(Pm - Psd), Phi = exp(Pm + Psd)))
}

sdReg <- function(Parray) {
  Pmean <- apply(Parray, c(1, 2), mean, na.rm = TRUE)
  Psd   <- apply(Parray, c(1, 2), sd,   na.rm = TRUE)
  return(list(Pmean = Pmean, Plow = Pmean - Psd, Phi = Pmean + Psd)) 
}


