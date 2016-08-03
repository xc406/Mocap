#'wrapper function to calculate area under the PR and ROC curves 
#'
#'\code{get.pr.roc} with the option to plot PR and ROC curves 
#'@param d a data frame with feature columns and gold standard column as the last column
#'@param tfct a character string specifying tf name and cell type
#'@param m sparse logistic regression model inherited from the LiblineaR package
#'@param plot logical
#'@param dir specifying a path to store plots in
#'@export
#'@examples
#'get.pr.roc(d,tfct,slr.model,plot=FALSE)
#'get.pr.roc(d,tfct,slr.model,plot=TRUE,dir=".")
#'@import LiblineaR
get.pr.roc <- function(d,tfct,m,plot=TRUE,dir="."){

        y=factor(d[,ncol(d)])
        s=scale(d[,c(1:ncol(d)-1)],center=TRUE,scale=TRUE)
        p=predict(m,s,proba=FALSE,decisionValues=TRUE)

        pred <- (-1)*p$decisionValues[,1]
        gs <- as.numeric(y) == 2
        aupr <- calc.aupr(pred,gs)
        res=table(p$predictions,y)

	if (plot){
        	pdf(paste(dir,"/pr",tfct,".pdf",sep=""))
        	par(mfrow=c(1,1))
        	par(mar=c(5.1,5.1,2.1,2.1),font.axis=2,font.lab=2,oma=c(0,0,2,0))
        	plot(aupr$rec,aupr$prec,type="l",lwd=3,col="blue",ylab="Precision",xlab="Recall")
       		legend("bottomleft",bty='n',legend=c(paste("AUPR=",round(aupr$AUPR,2)," Rec10=", round(aupr$rec10,2))),cex=1.2)
        	title(paste(tfct,sep=" "))
        	dev.off()

        	pdf(paste(dir,"/roc",tfct,".pdf",sep=""))
        	par(mar=c(5.1,5.1,2.1,2.1),font.axis=2,font.lab=2,oma=c(0,0,2,0))
        	plot(aupr$fpr,aupr$rec,type="l",lwd=3,col="blue",ylab="TPR (sensitivity)",xlab="FPR (1-specificity)")
        	legend("bottomright",bty='n',legend=c(paste("AUROC=",round(aupr$AUROC,2), " Fpr1=", round(aupr$fpr1,2))),cex=1.2)
       		title(paste(tfct,sep=" "))
        	dev.off()
	}
        return(c(aupr$AUPR,aupr$rec10,aupr$AUROC,aupr$fpr1))
}
