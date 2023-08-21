e_pcse<-function(explaining=0){
if(explaining==0) {
cat(" library(readxl) ", '\n')
cat(" Adata <- read_excel('Data40_PanelData.xlsx') ", '\n')
cat(" Re <- lm(사회복지~교부세+보조금, data=Adata) ", '\n')
cat(" ## PCSE 회귀분석 -------------- ", '\n')
cat(" pcse0 <- pcse(Re, Adata$Area, Adata$YEAR, pairwise = FALSE) ", '\n')
cat(" Re_pcse0 <- as.data.frame( cbind(pcse0$b, pcse0$pcse, pcse0$tstats ,pcse0$pval) ) ", '\n')
cat(" Re_pcse0 <- round(Re_pcse0, 4) ", '\n')
cat(" colnames(Re_pcse0)<-c('coef', 'pcse', 't-val', 'p-val') ", '\n')
cat(" Re_pcse0 ", '\n')
}}
