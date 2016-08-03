library(Mocap)##this will automatically load the example datasets to calculate X and Ys with

##command line prompts
print.usage <- function(){
        cat("
        USAGE:

        getXYmat.R --tf_name_cell_type RESTK562 \\
        --target_tf_name_cell_type MAXK562 \\
        --xdir /home/user/xdir \\
        --ydir /home/user/ydir

        ")
}

cmd.args <- commandArgs(trailingOnly = T)

args <- sapply(strsplit(cmd.args," "), function(i) i)
if (length(args) != 8){
        ##prompt usage
        stop(print.usage())
}else{
        tfct <- args[2]
        ptfct <- args[4]
        xdir <- args[6]
        ydir <- args[8]
}

##generate X mat
if(tfct>=ptfct){
        d1 <- get(paste0("data.",tfct))#model
        d2 <- get(paste0("data.",ptfct))#target
        xs <- feature.similarity(d1$data,d2$data)
        colnames(xs) <- paste(tfct,ptfct,sep="_")
        write.table(xs,file=paste(xdir,"/",tfct,ptfct,"X.txt",sep=""),row.names = TRUE,col.names = TRUE)
}

##generate Y mat: predicting d2 data with d1 features and d1-trained models
dx <- add.int(d2$data[,d1$fv+3])[,d1$fv.i]
##cross-sample predictions
aupr.x <- get.pr.roc(dx,paste(tfct,ptfct),d1$slr.model,plot=FALSE)

##MocapG predictions
pred.mocap <- d2$data$acces.score*d2$data$acces
gs <- d2$data$gs
aupr.mocap <- calc.aupr(pred.mocap,gs)

##MocapS predictions
d2.i <- add.int(d2$data[,d2$fv+3])[,d2$fv.i]
aupr.res <- get.pr.roc(d2.i,ptfct,d2$slr.model,plot=FALSE)

ys <- as.data.frame(c(aupr.res,aupr.x,aupr.mocap$AUPR,aupr.mocap$rec10,aupr.mocap$AUROC,aupr.mocap$fpr1))
colnames(ys) <- paste(tfct,ptfct,sep="_")
write.table(ys,file=paste(ydir,"/",tfct,ptfct,"Y.txt",sep=""),row.names = FALSE,col.names = TRUE)

