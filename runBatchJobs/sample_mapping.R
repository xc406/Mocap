library(Mocap)

##command line prompts
print.usage <- function(){
        cat("
        USAGE:

        sample_mapping.R --xdir /home/user/xdir 
        --ydir /home/user/ydir 
        --gamma 1.41 
        --ind.test 1 
	#--adjust 0.17 
	--cores 20 
	--exp \"ATAC-Seq\" 
	--verbose TRUE

        ")
}

cmd.args <- commandArgs(trailingOnly = T)

args <- sapply(strsplit(cmd.args," "), function(i) i)

if (length(args) != 14){
        ##prompt usage
        stop(print.usage())
}else{
        xdir <- args[2]
        ydir <- args[4]
        gamma <- as.numeric(args[6])
        ind.test <- as.numeric(args[8])
	#adjust <- as.numeric(args[10])
	cores <- as.numeric(args[10])
	exp <- args[12]
	verbose <- as.logical(args[14])
}


for (gamma in c(1.29,1.3,1.31,1.39,1.4,1.41,1.42)){
	res <- sample.mapping(xdir,ydir,gamma,ind.test,cores,exp,verbose)
	print(res)
}
