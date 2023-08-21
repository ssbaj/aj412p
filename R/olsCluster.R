#' olsCluster() Function
#' This function Clustered SE
#' Origin : https://www.fionaburlig.com/teaching/are212 
#' 

olsCluster <- function(data, y, X, clustervar) {

if (base::missing(data)) {
cat("  library(dplyr) ", "\n")
cat("  Const<-rep(1, nrow(df) )  ", "\n")
return(cat("  olsCluster(df, 'expend', c('Const', 'nfamily'), 'cluster변수') " )  )  }


tbldfGrabber <- function(data, varnames) {
matrixObject <- data %>% select(all_of(varnames)) %>% as.matrix()
return(matrixObject)
}

# usual setup
n <- nrow(data)
k <- length(X)

ydata <- tbldfGrabber(data, y)
xdata <- tbldfGrabber(data, X)
clusterdata <- tbldfGrabber(data, clustervar)

# grab each cluster identifier
clusters <- unique(clusterdata)
# number of clusters
G <- length(clusters)
betahat <- solve(t(xdata) %*% xdata) %*% t(xdata) %*% ydata
e <- ydata - xdata %*% betahat

# loop over all of the clusters
clusterMeat <- lapply(1:G, function(g) {
gindex <- which(clusterdata == clusters[g])
Xg <- xdata[gindex, ]
eg <- matrix(e[gindex, ])
XeeXg <- t(Xg) %*% eg %*% t(eg) %*% Xg
return(XeeXg)
})


# sandwich as per usual
meat <- (Reduce("+", clusterMeat)) / n
bread <- solve(t(xdata) %*% xdata)
# apply a DoF correction
dfc <- (G/(G-1))*((n-1)/(n-k))
sandwich <- dfc * n * t(bread) %*% meat %*% bread
# se as per usual
se <- sqrt(diag(sandwich))
output <- list(betahat, se)
names(output) <- c("betahat", "ClusteredSE")
return(output)

}

