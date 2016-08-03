#'Calculate AUPR and AUROC
#'
#'\code{calc.aupr} takes a list of predictions and gold standard 
#'and returns AUPR, AUROC, precision at 10% recall and 
#'sensitivity at 1% false positive rate.
#'@param pred a numeric vector of predictions 
#'@param gs a logical vector of gold standards specifying 
#'whether the predictions are in the gold standard
#'@export
calc.aupr <- function(pred, gs) {
  ord.idx <- order(pred, decreasing = T)

  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
  rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
  fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate

  prec <- c(prec[1], prec)
  rec <- c(0, rec)
  fpr <- c(0, fpr)
  rec10 <- mean(prec[which(round(rec,1)==0.1)])
  fpr1 <- mean(rec[which(round(fpr,2)==0.01)])

  aupr <- areaUnderCurve(rec, prec)
  auroc <- areaUnderCurve(fpr, rec)

  return(list(prec=prec, rec=rec, fpr = fpr, AUPR = aupr, AUROC = auroc, rec10 = rec10, fpr1 = fpr1))
}

#'Internal function
areaUnderCurve <- function(x, y) {
 dx <- diff(x)
 my <- y[1:(length(y) - 1)] + diff(y) / 2
 return(sum(dx * my))
}
