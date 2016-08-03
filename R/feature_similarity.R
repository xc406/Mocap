#'Chi-square permutation test
#'
#'\code{chisq.perm} function tests for categorical feature similarity and returns an empirical p-value of significance.
#'@param df1 a one-column data frame storing feature values
#'@param df2 a one-column data frame storing feature values
#'@param size an integer specifying sample size; default size to 10000
#'@param nperm an integer specifying number of permutations; default to 10000
#'@export
#'@examples
#'p.list <- chisq.perm(df1,df2)

chisq.perm <- function(df1,df2,size=10000,nperm=10000){
        chisqv <- vector()#testing as one mixture population
        chisqv.t <- vector()#testing as two separate populations
        pool <- c(df1,df2)
        for (i in 1:nperm){
                samp1 <- sample(1:length(df1),size,replace=FALSE)
                samp2 <- sample(1:length(df2),size,replace=FALSE)
                chisq.t <- chisq.test(cbind(table(df1[samp1]),table(df2[samp2])))
                chisqv.t <- c(chisqv.t,chisq.t$p.value)
        }
        for (i in 1:nperm){
                samp <- sample(1:length(pool),size*2,replace=FALSE)
                samp1 <- sample(samp,size,replace=FALSE)
                samp2 <- samp[-samp1]
                chisq <- chisq.test(cbind(table(pool[samp1]),table(pool[samp2])))
                chisqv <- c(chisqv,chisq$p.value)
        }
        p.emp <- (length(which(chisqv<=mean(chisqv.t)))+1)/(nperm+1)
        return(list(p = p.emp, chisqv = chisqv, chisqv.t = chisqv.t))
}

#'KS permutation test
#'
#'\code{ks.perm} function tests for similarity between continous features and returns an empirical p-value of significance.
#'@param df1 a one-column data frame storing feature values
#'@param df2 a one-column data frame storing feature values
#'@param size an integer specifying sample size; default size to 500
#'@param nperm an integer specifying number of permutations; default nperm to 10000
#'@export
#'@examples
#'p.list <- ks.perm(df1,df2)

ks.perm <- function(df1,df2,size=500,nperm=10000){
        ksv <- vector()#testing as one mixture population
        ksv.t <- vector()#testing as two separate populations
        pool <- c(df1,df2)
        for (i in 1:nperm){
                samp1 <- sample(1:length(df1),size,replace=FALSE)
                samp2 <- sample(1:length(df2),size,replace=FALSE)
                ks.t <- suppressWarnings(ks.test(df1[samp1],df2[samp2]))
                ksv.t <- c(ksv.t, ks.t$p.value)
        }
        for (i in 1:nperm){
                samp <- sample(1:length(pool),size*2,replace=FALSE)
                samp1 <- sample(samp,size,replace=FALSE)
                samp2 <- samp[-samp1]
                ks <- suppressWarnings(ks.test(pool[samp1],pool[samp2]))
                ksv <- c(ksv,ks$p.value)
        }
        p.emp <- (length(which(ksv<=mean(ksv.t)))+1)/(nperm+1)

        return(list(p = p.emp, ksv = ksv, ksv.t = ksv.t))
}

#'Wrapper function to assess feature similarity. 
#'
#'\code{feature.similarity} is a wrapper function to test for similarities between feature 
#'columns in two data frames d1 and d2.
#'@param d1 data frame 1
#'@param d2 data frame 2
#'@param verbose logical, optional parameter
#'@return a data frame of p-values
#'@export
#'@examples
#'feature.similarity(d1,d2)
feature.similarity <- function(d1,d2,verbose=TRUE) {
	#create an empty vector to hold p values
	pv <- vector()
	#remove chr,start,end,gs from data frame 
	d1 <- d1[,-c(1,2,3,ncol(d1))]
	d2 <- d2[,-c(1,2,3,ncol(d2))]
	nv <- intersect(colnames(d1),colnames(d2))
	for (f in nv){
      		df1 <- d1[,f]
      		df2 <- d2[,f]
		if (verbose) {
			cat("calculating similarity for feature", f,"\n")
		}
      		if (length(table(df1))<=6){
                	p <- chisq.perm(df1,df2)[[1]]
      		}else{
                	p <- ks.perm(df1,df2)[[1]]
        	}
        	names(p) <- f
        	pv <- c(pv,p)
    	}
	return(as.data.frame(pv))
}

