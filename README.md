# Mocap
This document explains how to use the package Mocap to infer transcription factor binding sites (TFBS) from chromatin accessibility (DNase-Seq or ATAC-Seq) data.

## Installation

```{r echo=FALSE}
library(devtools)
install_github("xc406/Mocap")
library(Mocap)
```
## Example Usage
Below are some examples of using MocapG and MocapX to infer cell type-specific TFBS.

#### Simulate DNase I cut count data and run MocapG.
```{r echo=FALSE}
#simulate DNase I cut count as a mixture of zero-inflated negative binomial distributions
count.all <- sim.count(num=10000,prob0.true=0.01, prob1.true=0.89,
		mu1.true=10, theta1.true=3, mu2.true=98, theta2.true=5)

#specify a set of random initial guesses and the number of bootstrap runs
num.boots<-parallel::detectCores()

size<-1000
mu1.init <- 0.1
mu2.init <- 10000
theta1.init <- 5.1
theta2.init <- 2.1
prob0.init <- 0.1
prob1.init <- 0.8

#EM estimation of mixture model
count.resp <- vector("list", num.boots)
for (bootstrap in 1:num.boots) {
	count.resp[[bootstrap]] <- count.all[sample(length(count.all),size,replace=TRUE)]
}
params <- parallel::mclapply(1:num.boots,FUN = em, count = count.resp,
     mu1 = mu1.init,theta1 = theta1.init,mu2 = mu2.init,theta2 = theta2.init,
     prob0 = prob0.init,prob1=prob1.init,reltol=1e-10,
     mc.cores=num.boots,mc.set.seed=FALSE,mc.preschedule=FALSE)

#use the model to classify synthetic BED-formatted motif sites 
data(bedCount)##chr,start,end,count
synthetic.data <- cbind(bedCount[sample(1:nrow(bedCount),length(count.all),replace=FALSE),c("chr","start","end")],count.all)
res.mocapG <- run.MocapG(synthetic.data,params)
```

To use your own accessibility cut data, substitute the `simCount` function call with the following steps:
   1. sample DNase-Seq or ATAC-Seq cut counts of the human/mouse genome in 200bp windows, or use cut counts from a genome-wide motif +/- 100bp windows excluding the motif site itself
   2. filter out low mappability windows (< 0.8) using e.g. wgEncodeCrgMapabilityAlign50mer.bigWig from [ENCODE](ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/)

We provide the `get.cut.count` function in the Mocap package to create your own dataset. Please refer to the next section (score calculation scripts) for details on how to obtain the necessary files to make the function call.

#### Generate motif-associated accessibility cut count data and run MocapG.
```{r echo=FALSE}
#to use an example dataset
data(bedCount)

#to get count data directly from files
bedCount <- get.cut.count(path="/data/cgsb/bonneau/xchen/motifs/scores/",
		chrom="chr12",tf="STAT3",ct="K562",n=1,seqlen=24)

#specify a set of random initial guesses and the number of bootstrap runs
num.boots<-parallel::detectCores()

size<-10000
mu1.init <- 0.1
mu2.init <- 10000
theta1.init <- 5.1
theta2.init <- 2.1
prob0.init <- 0.1
prob1.init <- 0.8

#EM estimation of mixture model
count.resp <- vector("list", num.boots)
for (bootstrap in 1:num.boots) {
        count.resp[[bootstrap]] <- bedCount$count[sample(length(bedCount$count),size,replace=TRUE)]
}
params <- parallel::mclapply(1:num.boots,FUN = em, count = count.resp,
     mu1 = mu1.init,theta1 = theta1.init,mu2 = mu2.init,theta2 = theta2.init,
     prob0 = prob0.init,prob1=prob1.init,reltol=1e-10,
     mc.cores=num.boots,mc.set.seed=FALSE,mc.preschedule=FALSE)

#classify BED-formatted motif sites
res.mocapG <- run.MocapG(bedCount,params)
```
*Note that larger sample size yields more consistent parameter estimations, but also requires more iterations to converge.*

---

#### To run MocapX

 1. preprocess an accessibility (ATAC-Seq or DNase-Seq) sample

    ```
    ##create a temporary directory to store temp files
    export MYTMP=/state/partition1/$USER/$$
    mkdir -p $MYTMP
    python extractCutSite.py example.bam chrom_size_file output_dir 0 paired
    rm -rf $MYTMP
    ```
    This will return one processed bam file if stranded==0, and two strand-specific bam files if stranded==1.
    Use "paired" for paired-end sequencing, "single" for single-end sequencing.

    To convert to a wiggle file using the bam2wig program in RSeQC 
    ```
    #install or load rseqc
    #module load rseqc

    bam2wig.py -i file_p_cut.bam -s chrom_size_file -o file_p -t 0
    ```

 2. run score calculation scripts to generate data matrices for TFBS predictions

    Specify output directory, e.g. outpath="/scratch/$USER/scores/" to store calculated scores.
    Motif bed files (for hg19 build) can be downloaded from [here](http://whisper.bio.nyu.edu/~xc406/mocap/data).
    Move /bedMotifs folder to your outpath before running the below scripts.

    ```
    ##cut counts
    python countWig.py motif-bed-file wig-file out-path ctName

    ##footprint scores
    python getFps.py wig-file out-path motifChrom ctName cutoff tfName expName

    ##conservation scores
    python getCons.py path-to-cons-files out-path motifChrom tfName

    ##mapability scores
    ##install or load kent to convert downloaded bigWig files into per chromosome bedGraph files
    #module load kent
    bash makeMap.sh /path/wgEncodeCrgMapabilityAlign100mer.bigWig
    python getMap.py path-to-map-files map-file-name out-path motifChrom seqlen window

    ##sequence features, GC/CpG
    python getGC.py path-to-fasta-files out-path motifChrom

    ##gene features (closest target gene), dist-to-tss
    python getDtss.py refseq-file motifChrom out-path
    ```
    `getMap.py`, `getGC.py` and `getDtss.py` will run through all TF motif files under the out-path.

    ```{r echo=FALSE}
    #to load the data for TFBS predictions in a data frame using files from path=outpath
    testData <- get.all.data(path="/data/cgsb/bonneau/xchen/motifs/scores/",
	chrom="chr12",tf="STAT3",ct="K562",n=1,params=params,seqlen=24,coords=FALSE,chip=FALSE)

    #to use example datasets
    testData <- data.CTCFK562$data
    ``` 

 3. get the feature similarity X matrices and cross-prediction performance Y matrices

    ```
    #to generate pbs scripts for running rscript genXYmat.R to calculate the X and Y matrices on a cluster
    #generate all pair-wise X and Y scores using existing samples
    bash makeXYmat0.sh

    #generate pair-wise scores between existing samples and a new input TFCelltype
    bash makeXYmat.sh TFCelltype

    #example for running pbs rscripts on the cluster
    Rscript genXYmat.R --tf_name_cell_type TFCelltype \
                       --pred_tf_name_cell_type CTCFHepg2 \
                       --xdir /data/cgsb/bonneau/xchen/mocap/xdir \
                       --ydir /data/cgsb/bonneau/xchen/mocap/ydir

    #Files to generate X and Y matrices will be stored in xdir and ydir.
    ```

    A set of precalculated X and Y files used in the paper can be downloaded from [here](http://whisper.bio.nyu.edu/~xc406/mocap/data/).
  
 4. run sample mapping script to find the MocapX model for the input sample

    ```
    #generate pbs scripts and run sample_mapping.R on a cluster to find model-sample mappings
    bash makeSampleMapping.sh

    Rscript sample_mapping.R --xdir "/data/cgsb/bonneau/xchen/mocap/xdir/" \
                             --ydir "/data/cgsb/bonneau/xchen/mocap/ydir/" \
                             --gamma 1.31 \
                             --ind.test 1 \
                             --cores 20 \
                             --exp DNase-Seq \
                             --verbose TRUE
    #To find model mapping for a new input sample set ind.test to 0
    ```

 5. use the model trained in TFCT1 for TFBS predictions in TFCT2
    ```{r echo=FALSE}
    ##example using sparse logistic regression model trained in YY1K562 to predict binding sites in ETS1K562
    tfct1 <- "YY1K562"
    tfct2 <- "ETS1K562"

    ##get data for trained model in tfct1
    d1 <- get(paste0("data.",tfct1))

    ##get data for predicting binding sites in tfct2
    d2 <- get(paste0("data.",tfct2))
    dx <- add.int(d2$data[,d1$fv+3])[,d1$fv.i]
    s <- scale(dx[,c(1:ncol(dx)-1)],center=TRUE,scale=TRUE)

    ##making predictions
    p <- predict(d1$slr.model,s,proba=FALSE,decisionValues=TRUE)
    res.mocapX <- (-1)*p$decisionValues[,1]
    z.scores <- (res.mocapX-mean(res.mocapX))/sd(res.mocapX)

    ##save predictions in a bed file, option to convert posterior scores to z scores
    write.table(cbind(d2$data[,c(1,2,3)],z.scores),file=paste0("~/results/",tfct1,tfct2,"MocapX.bed"),
	row.names=F,col.names=F,quote=F,sep="\t")   

    ##plot AUPR curves and save the plots to file
    aupr.x <- get.pr.roc(dx,paste(tfct1,tfct2,sep="_"),d1$slr.model,plot=TRUE,dir="~/plots/")
    ```

