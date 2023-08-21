\name{pcse}
\alias{pcse}
\title{Panel Corrected Standard Error}
\usage{
pcse()
}
\description{
PCSE
}
\examples{

# 출처: https://rdrr.io/cran/pcse/src/R/pcse.R
library(dplyr); library(readxl); library(aj412)
Adata<-read_excel('data40_paneldata.xlsx')
re<-lm(사회복지~교부세+보조금+인구밀도, data=Adata)
pcse(re, Adata$Area, Adata$YEAR)
}