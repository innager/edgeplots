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

getCond <- function(Parray) {
  cond.mat <- Parray[1, c(2:5, 7:10), ] 
  cond <- t(!is.na(cond.mat))   # rows are samples
  cond01 <- apply(cond, 1, function(x) paste(as.integer(x), collapse = ""))
  cond01 <- paste(substr(cond01, 1, 4), substr(cond01, 5, 8))
  tbl.cond <- sort(table(cond01), decreasing = TRUE)
  uniq.cond <- names(tbl.cond)
  combNA <- NULL
  for (cond in uniq.cond) {
    combNA <- cbind(combNA, cond01 == cond)
  }    
  colnames(combNA) <- uniq.cond
  return(list(comb = combNA, tbl = tbl.cond))
}

getUniqCond <- function(Parray) {
  cond.mat <- Parray[1, c(2:5, 7:10), ] 
  cond <- t(!is.na(cond.mat))   # rows are samples
  cond01 <- apply(cond, 1, function(x) paste(as.integer(x), collapse = ""))
  cond01 <- paste(substr(cond01, 1, 4), substr(cond01, 5, 8))
  tbl.cond <- sort(table(cond01), decreasing = TRUE)
  return(names(tbl.cond))
}

load("ParrayGam3lr.RData")
uniq.cond <- getUniqCond(Parray)
a <- getCond(Parray)
