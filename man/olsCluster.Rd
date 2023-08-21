\name{olsCluster}
\alias{olsCluster}
\title{Making STATA cluster SE}
\usage{
olsCluster()
}
\description{
Making STATA cluster SE
}
\examples{

library(dplyr); library(readxl)
Adata<-read_excel('Edata.xlsx')
olsCluster(Adata, '종속변수', c('ones', '독립변수'), 'Cluster표시변수' )

}