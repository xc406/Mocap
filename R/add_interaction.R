#'Function to add interactions
#'
#'\code{add.int} adds interaction features to a data frame
#'@param d.train a data frame with feature columns and gold standard column
#'@export
#'@examples
#'add.int(d.train)
add.int <- function(d.train){

        y <- d.train$gs
        d.train <- subset(d.train, select = -c(gs))
        n <- ncol(d.train)
        for (i in 1:n){
                for (j in i:n){
                        comb <- apply(d.train,1,function(x) x[i]*x[j])
                        name <- paste(colnames(d.train)[i],"x",colnames(d.train)[j],sep=" ")
                        d.train[,name] <- comb
                }
        }
        d.train[,"y"] <- y
        return(d.train)
}
