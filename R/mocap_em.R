#'simulate count data
#'
#'\code{sim.count} simulate cut count data from a mixture of zero-inflated negative 
#'binomial distributions
#'@param num an integer specifying the total number of counts
#'@param prob0.true a numeric value < 0.3
#'@param prob1.true a numeric value > 0.8
#'@param mu1.true a numeric value
#'@param theta1.true a numeric value
#'@param mu2.true a numeric value, where mu1.true << mu2.true
#'@param theta2.true a numeric value
#'@return a numeric vector of simulated cut counts
#'@export
sim.count <- function(num,prob0.true,prob1.true,mu1.true,theta1.true,mu2.true,theta2.true){
        ##simulate inaccessible component
        count1 <- rnbinom(num*(1-prob0.true)*prob1.true,mu=mu1.true,size=theta1.true)
        ##simulate accessible component
        count2 <- rnbinom(num*(1-prob0.true)*(1-prob1.true),mu=mu2.true,size=theta2.true)
        ##simulate zero component
        zero <- rep(0,num*prob0.true)

        count.all <- c(zero,count1,count2)
        prob0.zero <- length(zero)/length(which(count.all==0))
        return(count.all)
}

#'Internal function
#'
#'log-likelihood function; setting 0s to .Machine$double.xmin to avoid numerical underflow of log(0)
loglikfun <- function(count,mu1,theta1,mu2,theta2,prob0,prob1,zi0,zi1){
	zero <- .Machine$double.xmin
        prob0[which(prob0==0)] <- zero
        prob0.rest <- 1-prob0
        prob0.rest[which(prob0.rest==0)] <- zero
        prob1[which(prob1==0)] <- zero
        prob2 <- 1-prob1
        prob2[which(prob2==0)] <- zero
        zi0[which(zi0==0)] <- zero
        zi1[which(zi1==0)] <- zero
        zi2 <- 1-zi0-zi1
        zi2[which(zi2==0)] <- zero
        loglik0 <- log(prob0*(count==0)/zi0)
        loglik1 <- log(prob0.rest*prob1*dnbinom(count,mu=mu1,size=theta1,log=F)/zi1)
        loglik2 <- log(prob0.rest*prob2*dnbinom(count,mu=mu2,size=theta2,log=F)/zi2)
        loglik0[which(loglik0==-Inf)] <- log(zero)
        loglik1[which(loglik1==-Inf)] <- log(zero)
        loglik2[which(loglik2==-Inf)] <- log(zero)
        loglik <- -(sum(loglik0*zi0+loglik1*zi1+loglik2*(1-zi0-zi1)))
        return(list(loglik,-2*(loglik2-loglik1)))##- sign: a positive value to be minimized
}

#' Internal fucntion
#'
#' \code{estep} updates zis of each EM iteration.
estep <- function(count,prob0,prob1,theta1,mu1,theta2,mu2){
        zi0_estep <- (prob0*(count==0))/ (( prob0*(count==0)) + (1-prob0)*(prob1*dnbinom(count,size=theta1,mu=mu1) +
                (1-prob1)*dnbinom(count,size=theta2,mu=mu2) ) )
        zi1_estep <- (1-prob0)*prob1*dnbinom(count,size=theta1,mu=mu1) / ( prob0*(count==0) +
                (1-prob0)*(prob1*dnbinom(count,size=theta1,mu=mu1) +
                (1-prob1)*dnbinom(count,size=theta2,mu=mu2) ) )
        zi2_estep <- (1-prob0)*(1-prob1)*dnbinom(count,size=theta2,mu=mu2) / ( prob0*(count==0) +
                (1-prob0)*(prob1*dnbinom(count,size=theta1,mu=mu1) +
                (1-prob1)*dnbinom(count,size=theta2,mu=mu2) ) )
        return(list(zi0_estep,zi1_estep,zi2_estep))
}

#' EM algorithm to derive zero-inflated negative binomial mixture model parameters
#'
#' \code{em} function models cut count as a mixture of zero-inflated negative binomial distributions. 
#'\code{em} takes as input a set of initial parameter guesses and iteratively updates 
#'parameter estimation to maximize the log-likelihood function until convergence criteria 
#'\code{reltol} is met. 
#'@param count a list of numeric vectors of cut counts from motif/genomic sampling
#'@param b an integer specifying the bootstrap number, where 1 <= b <= length(count)
#'@param mu1 a numeric value given as the initial guess of the mean parameter of inaccessible component
#'@param theta1 a numeric value given as the initial guess of the shape paramter of inaccessible component
#'@param mu2 a numeric value given as the initial guess of the mean parameter of accessible component, requires mu1 << mu2
#'@param theta2 a numeric value given as the initial guess of the shape paramter of accessible component
#'@param prob0 a numeric value given as the initial guess of the probability of zero component
#'@param prob1 a numeric value given as the initial guess of the probability of inaccessible component
#'@param reltol a small numeric value specifying the relative tolerance for convergence
#'@return a list of numeric values containing estimates of mu1, theta1, mu2, theta2, prob0, prob1
#'@import MASS parallel
#'@export
#'@examples
#'count.all <- sim.count(num=10000,prob0.true=0.01, prob1.true=0.89, mu1.true=10, theta1.true=3, mu2.true=98, theta2.true=5)
#'##run as a single process
#'count <- vector("list",1)
#'size <- length(count.all)/2
#'count[[1]] <- count.all[sample(length(count.all),size,replace=FALSE)]
#'res <- em(count,b=1,mu1.init=2,theta1.init=2,mu2.init=200,theta2.init=5,prob0.init=0.05,prob1.init=0.95,1e-10)
#'
#'##run with multicores
#'num.boots<-detectCores()
#'mu1.init <- 0.1
#'mu2.init <- 10000
#'theta1.init <- 5.1
#'theta2.init <- 2.1
#'prob0.init <- 0.1
#'prob1.init <- 0.8
#'count.resp <- vector("list", num.boots)
#'for (bootstrap in 1:num.boots) {count.resp[[bootstrap]] <- count.all[sample(length(count.all),size,replace=TRUE)]}
#'res.mat <- mclapply(1:num.boots,FUN = em, count = count.resp,
#'mu1 = mu1.init,theta1 = theta1.init,mu2 = mu2.init,theta2 = theta2.init,prob0 = prob0.init,prob1=prob1.init,reltol=1e-10,
#'mc.cores=num.boots,mc.set.seed=FALSE,mc.preschedule=FALSE)

em <- function(count,b,mu1,theta1,mu2,theta2,prob0,prob1,reltol,verbose=FALSE){
        if (verbose) {cat("initial estimates",'bootstrap=',b,'mu1=',mu1,'theta1=',theta1,'mu2=',mu2,'theta2=',theta2,'pi0=',prob0,'pi1=',prob1,'\n')}
        weights <- rep(1,length(count[[b]]))
        zi <- estep(count[[b]],prob0,prob1,theta1,mu1,theta2,mu2)
        zi0 <- zi[[1]]
        zi1 <- zi[[2]]
        zi2 <- zi[[3]]
        ll_new <- loglikfun(count[[b]],mu1,theta1,mu2,theta2,prob0,prob1,zi0,zi1)[[1]]
        ll_old <- 2 * ll_new
        iter <- 0
        while(abs((ll_old - ll_new)/ll_old) > reltol) {
        	iter <- iter+1
                ll_old <- ll_new
                model_count1 <- suppressWarnings(glm.nb(count[[b]]~1 , weights = weights * zi1, start = log(mu1) , init.theta = theta1))
                model_count2 <- suppressWarnings(glm.nb(count[[b]]~1 , weights = weights * zi2, start = log(mu2), init.theta = theta2))
                mu1 <- exp(coefficients(model_count1))
                theta1 <- model_count1$theta
                mu2 <- exp(coefficients(model_count2))
                theta2 <- model_count2$theta

                model_zero <- suppressWarnings(glm.fit(as.integer(count[[b]]==0), zi0, weights = weights*(count[[b]]==0),family = binomial()))
                zi0 <- model_zero$fitted#sum(zi0)/sum(rep(1,length(zi0)))
                zi0[count[[b]]>0] <- 0
                ##update prob0, prob1
                prob0 <- sum(zi0)/sum(rep(1,length(count[[b]])))
                prob1 <- sum(zi1)/sum(1-zi0)
                ##update zi0, zi1, zi2
                zi0 <- estep(count[[b]],prob0,prob1,theta1,mu1,theta2,mu2)[[1]]
                zi1 <- estep(count[[b]],prob0,prob1,theta1,mu1,theta2,mu2)[[2]]
                zi2 <- estep(count[[b]],prob0,prob1,theta1,mu1,theta2,mu2)[[3]]
                ll_new <- loglikfun(count[[b]],mu1,theta1,mu2,theta2,prob0,prob1,zi0,zi1)[[1]]
                if (verbose) {cat(".",b,mu1,theta1,mu2,theta2,prob0,prob1,ll_new,iter,'\n')}
        }
        df <- data.frame(mu1,theta1,mu2,theta2,prob0,prob1)
        #write.table(df,file=output.file,append=TRUE,row.names=FALSE,col.names=FALSE,sep='\t')
        if (verbose) {cat("final estimates",'bootstrap=',b,'mu1=',mu1,'theta1=',theta1,'mu2=',mu2,'theta2=',theta2,'pi0=',prob0,'pi1=',prob1,'\n')}
        return(df)
}
