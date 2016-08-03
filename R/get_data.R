#'Function to load data
#'
#'\code{get.all.data} loads data for sparse logistic regression training and making predictions.
#'provides sampling of loaded data as an option.
#'@param path a directory specifying the path to all calculated scores 
#'@param chrom a character string, e.g. "chr12"
#'@param tf a character string, e.g. "CTCF"
#'@param ct a character string, e.g. "Hepg2"
#'@param n a numeric value
#'@param params a list of data frames returned by \code{em}
#'@param seqlen an integer specifying sequencing length, options are 24, 36, 40, 50, 75 and 100
#'@param coords logical
#'@param chip logical
#'@export 
#'@import fdrtool

get.all.data <- function(path,chrom,tf,ct,n,params,seqlen=24,coords=FALSE,chip=FALSE){

	setwd(path)
        map <- read.table(gzfile(paste("mapMotifs/",tf,chrom,"map",seqlen,".txt.gz",sep="")), header=F)
        motif <- read.table(gzfile(paste("bedMotifs/",tf,chrom,".bed.gz",sep="")), header=F)
        gc <- read.table(gzfile(paste("gcMotifs/",tf,chrom,"gc.txt.gz",sep="")), header=F)
        cons <- read.table(gzfile(paste("consMotifs/",tf,chrom,"cons.txt.gz",sep="")), header=F)
        fps <- read.table(gzfile(paste("fpsMotifs/",tf,ct,chrom,"fps.txt.gz",sep="")), header=F)
	bedCount <- cbind(motif[,c(1,2,3)],fps[,1])
	mocap <- run.Mocap(bedCount,params)
        target <- read.table(gzfile(paste("targetMotifs/",tf,chrom,"target.txt.gz",sep="")), header=F)

        if(coords){
                d <- cbind(motif[,c(1,2,3,5)],gc[,c(4,7,11,5,10)],mocap$acces,fps[,5],map,target[,2],cons,mocap$acces.score)
                colnames(d) <- c("chr","start","end","motif.score","repeats","gc.w1","gc.w4",
			"cpg.w1","cpg.wr4","acces","footprint.score","map","dtss","pc100v",
			"pc46v","pc46m","pc46p","pp100v","pp46v","pp46m","pp46p","acces.score")
        }else{
                d <- cbind(motif[,5],gc[,c(4,7,11,5,10)],mocap$acces,fps[,5],map,target[,2],cons,mocap$acces.score)
                colnames(d) <- c("motif.score","repeats","gc.w1","gc.w4","cpg.w1","cpg.wr4","acces",
                        "footprint.score","map","dtss","pc100v","pc46v","pc46m","pc46p","pp100v",
			"pp46v","pp46m","pp46p","acces.score")
        }

        d$tss <- as.numeric(as.factor(d$dtss < 1000))
        d$motif.score <- -log(d$motif.score)
        d$motif.score <- (d$motif.score - mean(d$motif.score))/mean(d$motif.score)
        d$repeats <- as.numeric(d$repeats)
        d$acces <- d$acces+1
        d$footprint.score <- d$footprint.score+1
        d$cpg.island <- as.numeric(d$cpg.wr4>0.6)*as.numeric(d$gc.w4>0.5)+1
        d <- subset(d, select = -c(cpg.wr4,gc.w4))
	if (chip){
		chip <- read.table(paste("chipMotifs/",tf,chrom,ct,"chip.txt",sep=""), header=F)
        	d$gs <- as.numeric(apply(chip, 1, function(x) all(x>0)))
	}
        map.fdr <- fdrtool(d$map,statistic="pvalue",plot=FALSE,color.figure=FALSE,verbose=FALSE,cutoff.method="fndr")
        d <- d[which(map.fdr$qval > map.fdr$param[1]),]
        d <- d[which(complete.cases(d)==TRUE),]
        samp <- sample(1:dim(d)[1],round(dim(d)[1]/n,0),replace=FALSE)
        d <- d[samp,]

        return(d)
}

#'function to get motif-associated accessibility cut count
#'
#'\code{get.cut.count}
#'@param path a directory path
#'@param chrom a character string, e.g. "chr12"
#'@param tf a character string, e.g. "CTCF"
#'@param ct a character string, e.g. "Hepg2"
#'@param n a numeric value
#'@param seqlen an integer, 24, 36, 50 or 100
#'@export
#'@import fdrtool

get.cut.count <- function(path,chrom,tf,ct,n,seqlen=24){

	setwd(path)
	bedCount <- read.table(gzfile(paste("countMotifs/",tf,chrom,ct,"count.txt.gz",sep="")), header=F)
	map <- read.table(gzfile(paste("mapMotifs/",tf,chrom,"map",seqlen,".txt.gz",sep="")), header=F)
	motif <- read.table(gzfile(paste("bedMotifs/",tf,chrom,".bed.gz",sep="")), header=F)
	#bedCount <- cbind(motif[,c(1,2,3)],count[,1])
	map.fdr <- fdrtool(map[,1],statistic="pvalue",plot=FALSE,color.figure=FALSE,verbose=FALSE,cutoff.method="fndr")
	bedCount <- bedCount[which(map.fdr$qval > map.fdr$param[1]),]
	colnames(bedCount) <- c("chr","start","end","count")

	samp <- sample(1:dim(bedCount)[1],round(dim(bedCount)[1]/n,0),replace=FALSE)
        bedCount <- bedCount[samp,]
	
	return(bedCount)
}
