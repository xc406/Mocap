#'cross-sample mapping function
#'
#'\code{sample.mapping} makes cross-sample mappings between models and samples. 
#'@param xdir a path to directory storing the x files 
#'@param ydir a path to directory storing the y files
#'@param gamma a numeric value
#'@param ind.test an integer
#'@param cores an integer
#'@param exp a character string, either "ATAC-Seq" or "DNase-Seq"
#'@param verbose logical
#'@export
#'@import caret parallel MASS
#'@examples
#'sample.mapping(xdir,ydir,gamma,ind.test,cores)
sample.mapping <- function(xdir,ydir,gamma,ind.test,cores,exp="ATAC-Seq",verbose=TRUE){

	tfct <- get_tfct(ydir)
	ntf <- length(tfct)
	ind.all <- seq(1,ntf)		

	##get original data for testing
 	ori.data <- select_mat(xdir,ydir,1,0)[[2]]
	##get original training data for cv steps
	ori.aupr.mat.scale.all <- select_mat(xdir,ydir,1,ind.test)[[2]]
	##get gamma-weighted training data
        data.train <- select_mat(xdir,ydir,gamma,ind.test)

        pvalv <- vector()
        goodv <- vector()
        badv <- vector()
        dsumv <- vector()
        bestPval <- 1
        bestRatio <- 1
	bestAvg <- 0.5
	bestDsum <- 0

	##hold-one-out cross-validation
	if (ind.test==0){
		cvflds <- createFolds(ind.all, k = ntf, list = TRUE, returnTrain = FALSE)
	}else{
        	cvflds <- createFolds(ind.all[-ind.test], k = cores, list = TRUE, returnTrain = FALSE)#98-length(ind.test), list = TRUE, returnTrain = FALSE)
	}

	try.c <- rev(seq(2,15,0.2))
        for (c in try.c){
                ##train data to select hyperparameter c through 10-fold cv
                kres <- mclapply(1:length(cvflds),FUN=cv,c=c,ind.test=ind.test,cvflds=cvflds,
			xdir=xdir,ydir=ydir,data.train=data.train,
			ori.aupr.mat.scale.all=ori.aupr.mat.scale.all,gamma=gamma,
			mc.cores=cores,mc.set.seed=FALSE,mc.preschedule=FALSE)
		##combine results from all cv-folds
                res <- do.call(rbind,kres)
		##compare cross-predictions with mocap defaults
                good <- length(which(res[,2]>res[,4]))
                bad <- length(which(res[,2]<res[,4]))

                avg <- colMeans(res[which(res[,2]!=res[,4]),])[2]
                dsum <- sum(res[,2]-res[,4])
                pval <- tryCatch({
                        t.test(res[which(res[,2]!=res[,4]),2],res[which(res[,2]!=res[,4]),4],paired=T,alternative="greater")$p.value
                },error=function(e){
                        1#print(paste("Error",length(which(res[,2]!=res[,4]))))
                })
                if (verbose) {cat(gamma,c,good,bad,avg,dsum,pval,colMeans(res),"\n")}
                ratio <- bad/good
                pvalv <- c(pvalv,pval)
                goodv <- c(goodv,good)
                badv <- c(badv,bad)
                dsumv <- c(dsumv,dsum)
                if (dsum > bestDsum){
                #if (pval < bestPval){
			bestRatio <- ratio
		        bestPval <- pval
		        bestC <- c
		        bestAvg <- avg
		        bestDsum <- dsum
		}
	}

	#use bestC to derive trained beta and w
	tf.df.t.train <- data.train[[1]]
	aupr.mat.scale.all.train <- data.train[[2]]

	if (exp=="ATAC-Seq"){
		tf.df.t.train.subset <- subset(tf.df.t.train, select = -c(motif.score,footprint.score))
	}else{
		tf.df.t.train.subset <- subset(tf.df.t.train, select = -c(motif.score))
	}
	rlm <- rlm(y~.,data=tf.df.t.train.subset,psi=psi.bisquare,c=bestC,maxit=1000)
	beta <- rlm$coefficients[-1]
	intercept <- rlm$coefficients[1]
	w <- rlm$w

	#get all data for predictions
	data <- select_mat(xdir,ydir,gamma,0)
	##get the training+testing X mat
	tf.df.t <- data[[3]]

	##no use of Ys
	if (exp=="ATAC-Seq"){
		tf.df.t.subset <- subset(tf.df.t,select=-c(motif.score,footprint.score,y))
	}else{
		tf.df.t.subset <- subset(tf.df.t,select=-c(motif.score,y))
	}
	d.tf <- dist(scale(tf.df.t.subset,scale=TRUE,center=TRUE))
	all.pairs <- rownames(tf.df.t)
	names(w) <- attr(rlm$x,"dimnames")[[1]]
	known.pairs <- names(w)
	unknown.pairs <- rownames(tf.df.t)[which(!rownames(tf.df.t)%in%known.pairs)]

	d.tf.mat <- as.matrix(d.tf)
	d.tf.known <- d.tf.mat
	for (i in c(1:dim(d.tf.mat)[1])){
		d.tf.known[i,i] <- NA
	}

	d.tf.known[unknown.pairs,] <- NA

	##assign a local weight
	nn <- apply(d.tf.known,2,function(x) rownames(d.tf.mat)[which(x==min(x,na.rm=T))])
	wv <- vector()
	for (i in c(1:length(unknown.pairs))){
		if (length(nn[[i]])>1){
			wv <- c(wv,estimate_mode(w[nn[[i]]]))
		}else{
			wv <- c(wv,w[nn[[i]]])
		}
	}
	names(wv) <- unknown.pairs
	nw <- c(w,wv)
	nw <- nw[sort(names(nw))]

	y.fit <- as.data.frame(intercept + nw * as.matrix(tf.df.t.subset) %*% beta)
	rownames(y.fit) <- rownames(tf.df.t)

	if (ind.test!=0){

		aupr.fitted.mat <- matrix(NA,ntf,ntf)
		rownames(aupr.fitted.mat) <- tfct
		colnames(aupr.fitted.mat) <- tfct
		for (i in c(1:length(tfct))){
			for (j in c(1:length(tfct))){
				if (tfct[i] > tfct[j]){
					aupr.fitted.mat[i,j] <- y.fit[paste(tfct[i],tfct[j],sep="_"),1]
				}else{
					aupr.fitted.mat[i,j] <- y.fit[paste(tfct[j],tfct[i],sep="_"),1]
		        	}
			}
		}

		rownames(aupr.fitted.mat) <- sapply(rownames(aupr.fitted.mat), function(x) paste(x,"fit",sep="_"))
		ind.rest <- ind.all[which(!ind.all %in% ind.test)]
		aupr.fitted.mat <- as.matrix(aupr.fitted.mat[ind.rest,])
		#default.scale <- aupr.mat.scale.all[nrow(aupr.mat.scale.all),]
		default.scale <- aupr.mat.scale.all.train[nrow(aupr.mat.scale.all.train),]

		tfct.d <- tfct[ind.rest]
		tfct.dfit <- c(tfct[ind.rest],"default_fit")
		tfl <- apply(aupr.fitted.mat,2,function(x) paste(tfct.d[which(x==max(x,na.rm=T))],"fit",sep="_"))
		tfl.values <- apply(aupr.fitted.mat,2,function(x) x[which(x==max(x,na.rm=T))])
		for (i in 1:length(tfl)){
			if (length(tfl[[i]])>1){
				default.values <- sapply(tfl[[i]], function(x) default.scale[[strsplit(x,"_fit")[[1]]]])
		        	tfl[[i]] <- tfl[[i]][which(default.values==max(default.values))]
		        	tfl.values[[i]] <- tfl.values[[i]][which(default.values==max(default.values))]
			}
			default.value <- default.scale[[strsplit(tfl[[i]],"_fit")[[1]]]]#max(c(mean(default.scale[ind.rest]),default.scale[[strsplit(tfl[[i]],"_fit")[[1]]]]))
			tfl.value <- tfl.values[[i]]/gamma#-0.17)
			if(default.value>tfl.value){tfl[[i]] <- "default_fit"}
		}
	
		##non-gamma-weighted training+testing data to compare results 
		ori.aupr.mat.scale <- ori.data[-nrow(ori.data),]
		ori.default.scale <- ori.data[nrow(ori.data),]
		aupr.real.mat.scale <- ori.aupr.mat.scale
		for (i in c(1:ntf)){
			aupr.real.mat.scale[i,i] <- NA
		}
		aupr.real.mat.d <- rbind(aupr.real.mat.scale[ind.rest,],ori.default.scale)
		rownames(aupr.real.mat.d) <- c(rownames(aupr.real.mat.scale)[ind.rest],"default")
		tfbl <- apply(aupr.real.mat.d,2,function(x) paste(tfct.d[which(x==max(x,na.rm=T))],"fit",sep="_"))

		default_fit <- default <- ori.default.scale
		aupr.all.mat.scale <- rbind(aupr.fitted.mat,ori.aupr.mat.scale,default,default_fit)
		res <- matrix(NA,ntf,6)##fitted,cross-validated,self,default,random
		for (i in c(1:length(tfl))){
			if (length(tfl[[i]])==1){
				res[i,1] <- aupr.all.mat.scale[tfl[[i]],names(tfl)[i]]
				res[i,2] <- aupr.all.mat.scale[strsplit(tfl[[i]],"_fit")[[1]][1],names(tfl)[i]]
			}else{
		        	res[i,1] <- mean(aupr.all.mat.scale[tfl[[i]],names(tfl)[i]])
				temp <- vector()
				for (j in c(1:length(tfl[[i]]))){
					temp <- c(temp,aupr.all.mat.scale[strsplit(tfl[[i]][j],"_fit")[[1]][1],names(tfl)[i]])
				}
				res[i,2] <- mean(temp)
			}
			res[i,3] <- aupr.all.mat.scale[names(tfl)[i],names(tfbl)[i]]
			res[i,4] <- aupr.all.mat.scale[nrow(aupr.all.mat.scale),names(tfl)[i]]
			rand.tf <- sample(1:length(tfct.d),1)
			res[i,5] <- mean(aupr.all.mat.scale[tfct.d[rand.tf],names(tfl)[i]])
			if (length(tfbl[[i]])==1){
				res[i,6] <- aupr.all.mat.scale[strsplit(tfbl[[i]],"_fit")[[1]][1],names(tfbl)[i]]
			}else{
				temp <- vector()
				for (j in c(1:length(tfbl[[i]]))){
					temp <- c(temp,aupr.all.mat.scale[strsplit(tfbl[[i]][j],"_fit")[[1]][1],names(tfbl)[i]])
				}
				res[i,6] <- mean(temp)
			}
		}

		rownames(res) <- tfct
		#pval <- t.test(res[which(res[,2]!=res[,4]),2],res[which(res[,2]!=res[,4]),4],paired=T,alternative="greater")$p.value
		res.test <- matrix(res[ind.test,],length(ind.test),6)
		rownames(res.test) <- rownames(res)[ind.test]
		tfl.test <- tfl[ind.test]

		assign(paste("res.test",ind.test,sep=""), res.test)
		assign(paste("res.fld",ind.test,sep=""), res)
		assign(paste("tfl.test",ind.test,sep=""), tfl.test)
		assign(paste("tfl.fld",ind.test,sep=""), tfl)
		res.all <- list(get(paste("res.fld",ind.test,sep="")),
			get(paste("tfl.fld",ind.test,sep="")),
			get(paste("res.test",ind.test,sep="")),
			get(paste("tfl.test",ind.test,sep="")),
			cvFolds=length(cvflds),bestC=bestC,bestDsum=bestDsum,bestPval=bestPval,c=try.c,dsum=dsumv,pval=pvalv,better=goodv,worse=badv)
		assign(paste("res",ind.test,sep=""),res.all)
		res.out <- res.all
		#save.image(paste("/data/cgsb/bonneau/xchen/motifs/xpv4/done/new/cvbisquare_98_97_",param3,"_",ind.test,"_meandm1c1",adjust,".RData",sep=""))
	}else{
		new.rows <- c(4754:4851)
		y.fit.value <- y.fit[new.rows,1]
        	y.fit.test <- rownames(y.fit)[new.rows]
        	res <- y.fit.test[which(y.fit[new.rows,1]==max(y.fit[new.rows,1]))]
        	res.value <- y.fit.value[which(y.fit[new.rows,1]==max(y.fit[new.rows,1]))]

		default.scale <- aupr.mat.scale.all[nrow(aupr.mat.scale.all),]
		default.value <- default.scale[[strsplit(res,"_")[[1]][1]]]
        	#default.value <- max(c(mean(default.scale),default.scale[[strsplit(res,"_")[[1]][1]]]))
        	#default.value <- default.scale[which(y.fit[c(4754:4851),1]==max(y.fit[c(4754:4851),1]))]
        	if (default.value > (res.value/(gamma))){
                	res.out <- paste("default",strsplit(res,"_")[[1]][2],sep="_")
        	}
		#save.image(paste("/data/cgsb/bonneau/xchen/motifs/xpv4/done/new/cvbisquare_98_98_",param3,"_",res,"_atac_",adjust,"_c1max.RData",sep=""))
	}
	return(res.out)
}

#'Internal function
#'
#'cross-validation function
#'@import MASS
cv <- function(k,c,ind.test,cvflds,xdir,ydir,gamma,data.train,ori.aupr.mat.scale.all,exp="ATAC-Seq"){
        tf.df.t.train <- data.train[[3]]
        aupr.mat.scale.all.train <- data.train[[2]]
        ind.rest <- c(1:98)[-ind.test]
        ind.cv <- cvflds[[k]]
        data.cv <- select_mat(xdir,ydir,gamma,c(ind.rest[ind.cv],ind.test))#get original ind in c(1:98)
        tf.df.t.cv <- data.cv[[1]]
        aupr.mat.scale.all.cv <- data.cv[[2]]

	##remove motif scores and footprint scores
	if (exp=="ATAC-Seq"){
		tf.df.t.cv.subset <- subset(tf.df.t.cv, select = -c(motif.score,footprint.score))
	}else{
		tf.df.t.cv.subset <- subset(tf.df.t.cv, select = -c(motif.score)) 
	}
        rlm <- rlm(y~.,data=tf.df.t.cv.subset,psi=psi.bisquare,c=c,maxit=1000)
        beta <- rlm$coefficients[-1]
        intercept <- rlm$coefficients[1]
        w <- rlm$w

	##get all training data for predictions
	if (exp=="ATAC-Seq"){
		tf.df.t.train.subset <- subset(tf.df.t.train,select = -c(motif.score,footprint.score,y))
	}else{
		tf.df.t.train.subset <- subset(tf.df.t.train,select = -c(motif.score,y))
	}
        d.tf.mat <- as.matrix(dist(scale(tf.df.t.train.subset,scale=TRUE,center=TRUE)))
        all.pairs <- rownames(tf.df.t.train)
        names(w) <- attr(rlm$x,"dimnames")[[1]]
        known.pairs <- names(w)
        unknown.pairs <- rownames(tf.df.t.train)[which(!rownames(tf.df.t.train)%in%known.pairs)]

        d.tf.known <- d.tf.mat
        for (i in c(1:dim(d.tf.mat)[1])){
                d.tf.known[i,i] <- NA
        }

        d.tf.known[unknown.pairs,] <- NA

        ##assign a local weight
        nn <- apply(d.tf.known,2,function(x) rownames(d.tf.mat)[which(x==min(x,na.rm=T))])
        wv <- vector()
        for (i in c(1:length(unknown.pairs))){
                if (length(nn[[i]])>1){
                        wv <- c(wv,estimate_mode(w[nn[[i]]]))
                }else{
                        wv <- c(wv,w[nn[[i]]])
                }
        }
        names(wv) <- unknown.pairs
        nw <- c(w,wv)
        nw <- nw[sort(names(nw))]

        y.fit <- as.data.frame(intercept + nw * as.matrix(tf.df.t.train.subset) %*% beta)
        rownames(y.fit) <- rownames(tf.df.t.train)

	tfct <- get_tfct(ydir)
	tfct <- tfct[-ind.test]

        aupr.fitted.mat <- matrix(NA,length(tfct),length(tfct))
        rownames(aupr.fitted.mat) <- tfct
        colnames(aupr.fitted.mat) <- tfct
        for (i in c(1:length(tfct))){
                for (j in c(1:length(tfct))){
                        if (tfct[i] > tfct[j]){
                                aupr.fitted.mat[i,j] <- y.fit[paste(tfct[i],tfct[j],sep="_"),1]
                        }else{
                                aupr.fitted.mat[i,j] <- y.fit[paste(tfct[j],tfct[i],sep="_"),1]
                        }
                }
        }

        rownames(aupr.fitted.mat) <- sapply(rownames(aupr.fitted.mat), function(x) paste(x,"fit",sep="_"))
        aupr.fitted.mat <- as.matrix(aupr.fitted.mat[-ind.cv,])

        #default.scale <- aupr.mat.scale.all.train[nrow(aupr.mat.scale.all.train),]
	default.scale.cv <- aupr.mat.scale.all.cv[nrow(aupr.mat.scale.all.cv),]

        tfct.d <- tfct[-ind.cv]
        tfct.dfit <- c(tfct[-ind.cv],"default_fit")
	##get the best fitted predictive models and values
        tfl <- apply(aupr.fitted.mat,2,function(x) paste(tfct.d[which(x==max(x,na.rm=T))],"fit",sep="_"))
        tfl.values <- apply(aupr.fitted.mat,2,function(x) x[which(x==max(x,na.rm=T))])
        for (i in 1:length(tfl)){
		##if best fitted models are more than one, choose the one with higher default value
                if (length(tfl[[i]])>1){
                        default.values <- sapply(tfl[[i]], function(x) default.scale.cv[[strsplit(x,"_fit")[[1]]]])
                        tfl[[i]] <- tfl[[i]][which(default.values==max(default.values))]
                        tfl.values[[i]] <- tfl.values[[i]][which(default.values==max(default.values))]
                }
                default.value <- default.scale.cv[[strsplit(tfl[[i]],"_fit")[[1]]]]#max(c(mean(default.scale.cv),default.scale.cv[[strsplit(tfl[[i]],"_fit")[[1]]]]))
                tfl.value <- tfl.values[[i]]/(gamma)
                if(default.value>tfl.value){tfl[[i]] <- "default_fit"}
        }

        ##get original training data (no gamma weighting)
        ori.aupr.mat.scale <- ori.aupr.mat.scale.all[-nrow(ori.aupr.mat.scale.all),]
        ori.default.scale <- ori.aupr.mat.scale.all[nrow(ori.aupr.mat.scale.all),]
        aupr.real.mat.scale <- ori.aupr.mat.scale
        for (i in c(1:length(tfct))){
                aupr.real.mat.scale[i,i] <- NA
        }

	##get original cv data (no gamma weighting)
        aupr.real.mat.d <- rbind(aupr.real.mat.scale[-ind.cv,],ori.default.scale)
        rownames(aupr.real.mat.d) <- c(rownames(aupr.real.mat.scale)[-ind.cv],"default")
	##true best
        tfbl <- apply(aupr.real.mat.d,2,function(x) paste(tfct.d[which(x==max(x,na.rm=T))],"fit",sep="_"))

	default_fit <- default <- ori.default.scale
        aupr.all.mat.scale <- rbind(aupr.fitted.mat,ori.aupr.mat.scale,default_fit,default)

        res <- matrix(NA,length(tfct),6)##fitted,cross-validated,self,default,random
        for (i in c(1:length(tfl))){
                if (length(tfl[[i]])==1){
                        res[i,1] <- aupr.all.mat.scale[tfl[[i]],names(tfl)[i]]
                        res[i,2] <- aupr.all.mat.scale[strsplit(tfl[[i]],"_fit")[[1]][1],names(tfl)[i]]
                }else{
                        res[i,1] <- mean(aupr.all.mat.scale[tfl[[i]],names(tfl)[i]])
                        temp <- vector()
                        for (j in c(1:length(tfl[[i]]))){
                                temp <- c(temp,aupr.all.mat.scale[strsplit(tfl[[i]][j],"_fit")[[1]][1],names(tfl)[i]])
                        }
                        res[i,2] <- mean(temp)
                }
                res[i,3] <- aupr.all.mat.scale[names(tfl)[i],names(tfbl)[i]]
                res[i,4] <- aupr.all.mat.scale[nrow(aupr.all.mat.scale),names(tfl)[i]]
                rand.tf <- sample(1:length(tfct.d),1)
                res[i,5] <- mean(aupr.all.mat.scale[tfct.d[rand.tf],names(tfl)[i]])
                if (length(tfbl[[i]])==1){
                        res[i,6] <- aupr.all.mat.scale[strsplit(tfbl[[i]],"_fit")[[1]][1],names(tfbl)[i]]
                }else{
                        temp <- vector()
                        for (j in c(1:length(tfbl[[i]]))){
                                temp <- c(temp,aupr.all.mat.scale[strsplit(tfbl[[i]][j],"_fit")[[1]][1],names(tfbl)[i]])
                        }
                        res[i,6] <- mean(temp)
                }
        }

        rownames(res) <- tfct
        resl <- as.data.frame(res)[ind.cv,]
        rownames(resl) <- tfct[ind.cv]
        return(resl)
}#cv function

#'Internal function
get_tfct <- function(ydir){
	lf.y <- list.files(ydir,pattern="x4.txt")
	diagf <- vector()
	for (f in lf.y){
		tf.names <- strsplit(f, "x4.txt")[[1]][1]
		tf1 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][1]
		tf2 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][2]
		ct1 <- strsplit(strsplit(tf.names, tf1)[[1]][2],tf2)[[1]][1]
		ct2 <- tail(strsplit(tf.names, tf2)[[1]],n=1)
		        
		if (ct1==ct2 & tf1==tf2){
		        diagf <- c(diagf,tf.names)
		        #cat(paste(f,"\n",sep=""))
		}
	}

	tfct <- sapply(diagf,function(x) regmatches(x, regexpr("K562|Hepg2|A549", x), invert = TRUE)[[1]][2])
	names(tfct) <- NULL
	return(tfct)
}

#'Internal fucntion
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
