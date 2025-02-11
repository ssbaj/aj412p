\name{olsNW} 
\alias{olsNW}
\title{Newey-West Standard Error}
\usage{
olsNW()
}
\description{
Newey-West Standard Error
}
\examples{

optimalLags <- (nrow(Adata))^(1/4)
olsNW(Adata, "종속변수", c("ones", "독립변수"), optimalLags)

library(sandwich)
cannedNW <- NeweyWest(wheatModel, lag = optimalLags, prewhite = FALSE, sandwich = TRUE)
cannedNWSE <- sqrt(diag(cannedNW))

}