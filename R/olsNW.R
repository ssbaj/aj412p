# Origin : https://www.fionaburlig.com/teaching/are212 

olsNW <- function(data, y, X, L) {


if (base::missing(data)) {
	    cat(" library(dplyr); library(HistData); WheatData<-Wheat ", "\n")
		cat(" Const<-rep(1, nrow(WheatData) ) ", "\n")
		cat(" WheatData<-cbind(WheatData, Const) ", "\n" )
		cat(" optimalLags<- nrow(WheatData)^(1/4) ", "\n" )
        return( cat(" myNW <- olsNW(WheatData, 'Wheat', c('Const', 'Wages'), optimalLags) " ) ) }


tbldfGrabber <- function(data, varnames) {
matrixObject <- data %>% select(all_of(varnames)) %>% as.matrix()
return(matrixObject)
}


# usual setup (note now we're calling n t instead)
# hashtag timeseries mindset

t <- nrow(data)
k <- length(X)

ydata <- tbldfGrabber(data, y)
xdata <- tbldfGrabber(data, X)
betahat <- solve(t(xdata) %*% xdata) %*% t(xdata) %*% ydata

# note that we now need our residuals to be a vector not a mat
# so that we can scalar-multiply them by our xmatrix

  e <- as.vector(ydata - xdata %*% betahat)

# emat is the scalar product of our residuals with the x mat

  emat <- e * xdata

# sets up our Newey-West weights -
# they will start at 1, and max out at 1 - l/(L + 1)

  weights <- seq(1, 0, by = -(1/(L + 1)))

# initialize our eMat'eMat matrix - note NOT = e'e
# scale by 1/2 bc we're only estimating half of the matrix

  eMatpeMat <- 0.5 * crossprod(emat) * weights[1]

# create the lag-weighted terms we're going to add to our main mat
# this spits out a list (note: lapply not sapply) of matrices

  weightStep <- lapply(2:length(weights), function(l) {

# crossprod(x,y) is just like t(x) %*% y, but faster

  weights[l] * crossprod(emat[1:(t - l + 1), ], emat[l:t,])})
  
  
# now actually add up all of the weighted summation terms
# Reduce() consecutively applies a fxn to each element of a list

  eMatpeMat <- eMatpeMat + Reduce("+", weightStep)

# fully populate the matrix - so far we've just estimated the lower diagonal

  eMatpeMat <- eMatpeMat + t(eMatpeMat)

# create the "meat" of the sandwich by multiplying this big matrix by 1/t

  meat <- eMatpeMat/t

# create the bread (as usual)

  bread <- solve(t(xdata) %*% xdata)

# finally, the full sandwich: multiply by t to get everything to work out

  sandwich <- t * (t(bread) %*% meat %*% bread)

# as usual, se's are the sqrt of the diagonal of the sandwich

  se <- sqrt(diag(sandwich))

# finally done!

  output <- list(betahat, sandwich, se)
  names(output) <- c("betahat", "sandwich", "se")
  return(output)
}

