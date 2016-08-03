#'select matrices
#'
#'\code{select_mat} constructs  matrices for weighted least squares regression.
#'@param xdir a directory containing pair-wise feature similarity files
#'@param ydir a directory containing cross-prediction performance files
#'@param gamma a numeric value specifying a hyperparameter
#'@param ind.test a vector of integers specifying hold-out data
#'@import LiblineaR parallel caret MASS
#'@export
select_mat <- function(xdir,ydir,gamma,ind.test){

	##design y matrix
	lf.y <- list.files(ydir,pattern="x4.txt")
	##extract diagonal Ys
	diagf <- vector()##store diagonal file names, e.g. "TFCtY.txt"
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

	diagtf.df <- data.frame()##store diagonal Ys
	for (fname in diagf){
        	tf <- read.table(paste(ydir,"/",fname,"x4.txt",sep=""), header=T)
                tfb <- log(tf[5,])##cross-prediction AUPR
                names(tfb) <- names(tf)
                if (length(diagtf.df)==0){
                        diagtf.df <- tfb
                }else{
                        diagtf.df <- c(diagtf.df,tfb)
                }
	}

	tfct <- paste(sapply(names(diagtf.df), function(x) strsplit(x,"_")[[1]][1]),
		sapply(names(diagtf.df), function(x) strsplit(x,"_")[[1]][2]),sep="")
	names(tfct) <- NULL

	tfct.all <- tfct
	if(ind.test[1]!=0){
		ind.all <- seq(1,length(tfct))
		tfct <- tfct.all[-ind.test]
	}
	tf.df.y <- data.frame()
	for (f in lf.y){
        	tf.names <- strsplit(f, "x4.txt")[[1]][1]
                tf1 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][1]
                tf2 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][2]
                ct1 <- strsplit(strsplit(tf.names, tf1)[[1]][2],tf2)[[1]][1]
                ct2 <- tail(strsplit(tf.names, tf2)[[1]],n=1)

                tfct1 <- paste(tf1,ct1,sep="")
                tfct2 <- paste(tf2,ct2,sep="")
                if (!tf.names %in% diagf & tfct2 %in% tfct & tfct1 %in% tfct){
                        tf <- read.table(paste(ydir,f,sep="/"), header=T)
			tfb <- log(tf[5,])
                        names(tfb) <- paste(paste(strsplit(names(tf),"_")[[1]][1],strsplit(names(tf),"_")[[1]][2],sep=""),
				      paste(strsplit(names(tf),"_")[[1]][3],strsplit(names(tf),"_")[[1]][4],sep=""),sep="_")
                	if (length(tf.df.y)==0){
                        	tf.df.y <- tfb
                	}else{
                        	tf.df.y <- c(tf.df.y,tfb)
                	}
        	}
	}

	##design and normalize x matrix
	lf.x <- list.files(xdir,pattern="pv3.txt")
	diagtf.df.x <- data.frame()
	for (fname in diagf){
        	tf1 <- strsplit(fname, "K562|Hepg2|A549")[[1]][1]
        	tf2 <- strsplit(fname, "K562|Hepg2|A549")[[1]][2]
        	ct1 <- strsplit(strsplit(fname, tf1)[[1]][2],tf2)[[1]][1]
        	ct2 <- tail(strsplit(fname, tf2)[[1]],n=1)

        	tfct1 <- paste(tf1,ct1,sep="")
        	tfct2 <- paste(tf2,ct2,sep="")
                tf <- read.table(paste(xdir,"/",fname,"pv3.txt",sep=""), header=F, row.names=1)
                if (length(diagtf.df.x)==0){
                        diagtf.df.x <- log(tf)
                        tf.names <- fname
                }else{
                        tf.names <- c(tf.names,fname)
                        diagtf.df.x <- cbind(diagtf.df.x,log(tf[,1]))
                }
	}

	colnames(diagtf.df.x) <- tf.names

	tf.df.x <- data.frame()
	for (f in lf.x){
        	tf.names <- strsplit(f, "pv3.txt")[[1]][1]
        	tf1 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][1]
        	tf2 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][2]
        	ct1 <- strsplit(strsplit(tf.names, tf1)[[1]][2],tf2)[[1]][1]
        	ct2 <- tail(strsplit(tf.names, tf2)[[1]],n=1)
        	tfct1 <- paste(tf1,ct1,sep="")
        	tfct2 <- paste(tf2,ct2,sep="")
        	if (!tf.names %in% diagf & tfct1 %in% tfct & tfct2 %in% tfct){
                	tf <- read.table(paste(xdir,f,sep="/"), header=F, row.names=1)
                	colnames(tf) <- paste(tfct1,tfct2,sep="_")
                	diag.means <- rowMeans(diagtf.df.x)
                	diag.sds <- apply(diagtf.df.x,1,sd)
                	v <- ifelse((log(tf) >= diag.means-diag.sds*8),1,0)
                	colnames(v) <- paste(paste(tf1,ct1,sep=""),paste(tf2,ct2,sep=""),sep="_")
                	if (length(tf.df.x)==0 & tfct1 > tfct2){
                        	tf.df.x <- v
                	}else{
                        	if (tfct1 > tfct2){
                                	tf.df.x <- cbind(tf.df.x,v)
                        	}
                	}
        	}
	}

	#normalize y matrices
	ntf <- length(tfct)
	aupr.mat <- matrix(NA,ntf,ntf)
	rownames(aupr.mat) <- tfct
	colnames(aupr.mat) <- tfct
	for (i in c(1:length(tfct))){
        	for (j in c(1:length(tfct))){
                	aupr.mat[i,j] <- tf.df.y[paste(tfct[i],tfct[j],sep="_")]
        	}
	}

	tf.df.diag <- data.frame()
	weight <- vector()
	default <- vector()
	for (f in lf.y){
        	tf.names <- strsplit(f, "x4.txt")[[1]][1]
        	tf1 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][1]
        	tf2 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][2]
        	ct1 <- strsplit(strsplit(tf.names, tf1)[[1]][2],tf2)[[1]][1]
        	ct2 <- tail(strsplit(tf.names, tf2)[[1]],n=1)
        	tfct1 <- paste(tf1,ct1,sep="")
        	tfct2 <- paste(tf2,ct2,sep="")
        	if (tf.names %in% diagf & tfct1 %in% tfct & tfct2 %in% tfct){
                	tf <- read.table(paste(ydir,f,sep="/"), header=T)
                	tfb <- log(tf[5,])
                	tfd <- log(tf[9,])##Mocap performance
			if (tf[11,] < 0.8){##control for general auroc
				tfw <- gamma##binarize Y mat
			}else{
				tfw <- ifelse(as.numeric(tf[5,]/tf[9,]>1.1),1,gamma)
			}
                	names(tfb) <- paste(paste(strsplit(names(tf),"_")[[1]][1],strsplit(names(tf),"_")[[1]][2],sep=""),
				      paste(strsplit(names(tf),"_")[[1]][3],strsplit(names(tf),"_")[[1]][4],sep=""),sep="_")
                	names(tfd) <- paste(paste(strsplit(names(tf),"_")[[1]][1],strsplit(names(tf),"_")[[1]][2],sep=""),
				      paste(strsplit(names(tf),"_")[[1]][3],strsplit(names(tf),"_")[[1]][4],sep=""),sep="_")
                	if (length(tf.df.diag)==0){
                        	tf.df.diag <- tfb
                        	default <- tfd
				weight <- tfw
                	}else{
                        	tf.df.diag <- c(tf.df.diag,tfb)
                        	default <- c(default,tfd)
				weight <- c(weight,tfw)
                        	#cat(paste(f,"\n",sep=""))
                	}
        	}
	}

	for (i in c(1:length(tfct))){
        	aupr.mat[i,i] <- tf.df.diag[paste(tfct[i],tfct[i],sep="_")]
	}

	names(default) <- colnames(aupr.mat)

	aupr.mat.weight <- aupr.mat*weight	
	aupr.mat.scale.all.weight <- scale(rbind(aupr.mat.weight,default),center = TRUE, scale = FALSE)
	##select sample pairs for regression
	aupr.mat.scale.select <- apply(aupr.mat.scale.all.weight,2,function(x) ifelse(x > -2, x, -300))	
	aupr.mat.scale.select[nrow(aupr.mat.scale.select),] <- aupr.mat.scale.all.weight[nrow(aupr.mat.scale.all.weight),]

	aupr.mat.scale.all <- scale(rbind(aupr.mat,default),center = TRUE, scale = FALSE)#original mat		
	default.scale <- aupr.mat.scale.all[nrow(aupr.mat.scale.all),]
	names(default.scale) <- colnames(aupr.mat.scale.all)
	aupr.mat.scale <- aupr.mat.scale.all[-nrow(aupr.mat.scale.all),]

	tf.df.y.scale <- vector()
	tf.df.y.scale.all <- vector()
	for (f in lf.y){
        	tf.names <- strsplit(f, "x4.txt")[[1]][1]
        	tf1 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][1]
        	tf2 <- strsplit(tf.names, "K562|Hepg2|A549")[[1]][2]
        	ct1 <- strsplit(strsplit(tf.names, tf1)[[1]][2],tf2)[[1]][1]
        	ct2 <- tail(strsplit(tf.names, tf2)[[1]],n=1)

        	tfct1 <- paste(tf1,ct1,sep="")
        	tfct2 <- paste(tf2,ct2,sep="")
        	if (!tf.names %in% diagf & tfct2 %in% tfct & tfct1 %in% tfct){
                	tfb <- aupr.mat.scale.select[tfct1,tfct2]
                	names(tfb) <- paste(tfct1,tfct2,sep="_")
                	if (tfct1 < tfct2 & !paste(tfct2,tfct1,sep="_") %in% names(tf.df.y.scale)){
                        	tf.df.y.scale[paste(tfct2,tfct1,sep="_")] <- aupr.mat.scale.select[tfct1,tfct2]
				tf.df.y.scale.all[paste(tfct2,tfct1,sep="_")] <- aupr.mat.scale.all[tfct1,tfct2]
                	}else if (tfct1 > tfct2 & !paste(tfct1,tfct2,sep="_") %in% names(tf.df.y.scale)){
                        	tf.df.y.scale[paste(tfct1,tfct2,sep="_")] <- aupr.mat.scale.select[tfct1,tfct2]
   				tf.df.y.scale.all[paste(tfct1,tfct2,sep="_")] <- aupr.mat.scale.all[tfct1,tfct2]
	             	}else if (tfct1 < tfct2 & paste(tfct2,tfct1,sep="_") %in% names(tf.df.y.scale)){
                        	tf.df.y.scale[paste(tfct2,tfct1,sep="_")] <- (tf.df.y.scale[paste(tfct2,tfct1,sep="_")]+aupr.mat.scale.select[tfct1,tfct2])/2
				tf.df.y.scale.all[paste(tfct2,tfct1,sep="_")] <- (tf.df.y.scale.all[paste(tfct2,tfct1,sep="_")]+aupr.mat.scale.all[tfct1,tfct2])/2
                	}else if (tfct1 > tfct2 & paste(tfct1,tfct2,sep="_") %in% names(tf.df.y.scale)){
                        	tf.df.y.scale[paste(tfct1,tfct2,sep="_")] <- (tf.df.y.scale[paste(tfct1,tfct2,sep="_")]+aupr.mat.scale.select[tfct1,tfct2])/2
				tf.df.y.scale.all[paste(tfct1,tfct2,sep="_")] <- (tf.df.y.scale.all[paste(tfct1,tfct2,sep="_")]+aupr.mat.scale.all[tfct1,tfct2])/2
                	}
                        #cat(paste(f,"\n",sep=""))
        	}
	}

	#combine x and y matrices
	tf.df.x <- tf.df.x[,sort(colnames(tf.df.x))]
	tf.df.y.scale <- tf.df.y.scale[sort(names(tf.df.y.scale))]
	tf.df <- rbind(tf.df.x,tf.df.y.scale)
	rownames(tf.df)[nrow(tf.df)] <- "y"
	tf.df.t <- as.data.frame(t(tf.df))
	tf.df.t.select <- tf.df.t[which(tf.df.t[,ncol(tf.df.t)]>-100),]

	tf.df.y.scale.all <- tf.df.y.scale.all[sort(names(tf.df.y.scale.all))]
	tf.df.all <- rbind(tf.df.x,tf.df.y.scale.all)
        rownames(tf.df.all)[nrow(tf.df.all)] <- "y"
        tf.df.t.all <- as.data.frame(t(tf.df.all))

	return(list(select.mat=tf.df.t.select,aupr.mat.scale.all=aupr.mat.scale.all,all.mat=tf.df.t.all))
}
