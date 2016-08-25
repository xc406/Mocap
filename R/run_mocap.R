#'Internal function
#'
#'Parse params df to list
#'@param params b x 6 data frame
getparams <- function(params){
	if (is.data.frame(params)){
		if (dim(params)[1]!=1) message("Warning: Input params should be a one-column data frame")
	}else{
		params <- do.call(rbind,params)
	}
	b <- dim(params)[1]#number of bootstrap
	p <- colMeans(params)
	return(list(mu1.est=p[1],
		theta1.est=p[2],
		mu2.est=p[3],
		theta2.est=p[4],
		prob0.est=p[5],
		prob1.est=p[6]))
}

#'generate motif/genomic region-assoicated accessibility scores
#'
#' \code{run.MocapG} takes as input a BED-4 formatted data frame and \code{em}
#'returned parameter set and output a cutoff value, a vector of binary classifications of 
#'genomic (motif) regions and a vector of log-likelihood accessibility scores.
#'  
#'@param bedCount a n x 4 data frame, where columns specify chr, start, end, count
#'@param params a b x 6 data frame containing \code{em} output parameter set
#'@return a list containing an integer cutOff value, 
#'a numeric vector of binary accessibility classifications and
#'a numeric vector of accessibility scores
#'@export
run.MocapG <- function(bedCount,params){
	params.list <- getparams(params)
	count.seq <- seq(0,range(bedCount[,4])[2])
	pred.close <- dnbinom(count.seq,mu=params.list$mu1.est,size=params.list$theta1.est)
	pred.open <- dnbinom(count.seq,mu=params.list$mu2.est,size=params.list$theta2.est)
	ll <- log(((1-params.list$prob1.est)*pred.open)/(params.list$prob1.est*pred.close))
	d <- cbind(count.seq,ll)

	for (i in 2:length(d[,2])){
	        if (d[i,2] > log(2) && d[i,2] > d[i-1,2]){
	                cutoff <- d[i,1]
	             	names(cutoff) <- NULL
			break
	        }
	}
	acces.score <- unlist(sapply(bedCount[,4], function(x) ll[which(count.seq==round(x))]))
	##correct Inf
	highest <- d[which(ll==max(acces.score[which(acces.score!=Inf)])),1]
	acces.score[which(acces.score==Inf)] <- max(acces.score[which(acces.score!=Inf)])-highest+bedCount[which(acces.score==Inf),4]
	lowest <- d[which(ll==min(ll)),1]
	##if zero-inflated assign zeros with the lowest score
	if (params.list$prob0.est > 1e-3){
		acces.score[which(bedCount[,4]<lowest)] <- min(acces.score)-lowest+bedCount[which(bedCount[,4]<lowest),4]
	}
	acces <- as.numeric(bedCount[,4] >= cutoff)
	return(list(cutoff=cutoff,acces=acces,acces.score=acces.score))
}




