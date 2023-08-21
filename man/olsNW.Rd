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
olsNW(Adata, "���Ӻ���", c("ones", "��������"), optimalLags)

library(sandwich)
cannedNW <- NeweyWest(wheatModel, lag = optimalLags, prewhite = FALSE, sandwich = TRUE)
cannedNWSE <- sqrt(diag(cannedNW))

}