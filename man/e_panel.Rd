\name{e_panel}
\alias{e_panel}
\title{Calc: Panel test}
\usage{
e_panel()
}
\description{
Calculating Panel test
}
\examples{

 library(readxl); library(plm)  
 Adata<-read_excel('Data40_paneldata.xlsx')  
 colnames(Adata)<-c('index','year','area','sb','gyobu','bojo','pdens','age65')  
 panel.set<-pdata.frame(Adata, index=c('area', 'year'))  
   
 ## 자료에 자연로그를 취함  
 Adata$lnsb<-log(Adata$sb)  
 Adata$lngyobu<-log(Adata$gyobu)  
 Adata$lnbojo<-log(Adata$bojo)  
 Adata$lnpdens<-log(Adata$pdens)  
   
 ## 지역별 이분산성  
 library(gplots)  
 plotmeans(Adata$lnsb~Adata$area)  
   
 ## 시간별 이분산성  
 plotmeans(Adata$lnsb~Adata$year)  
   
 *** R의 plm패키지에서 각 모형을 나타내는 명령문 ***  
 random.method = 'swar' / 'amemiya' / 'walhus' / 'nerlove'  
 pooled OLS (model = 'pooling') / fixed effects (model = 'within') / random effects (model = 'random')  
 first differences (model = 'fd') / between (model = 'between')  
   
 ## 모형의 설정  
 panel.model01 <- lnsb~lngyobu+lnbojo+lnpdens+age65  
 panel.model02 <- lnsb~lngyobu+lnbojo+lnpdens+age65+factor(year)  
 
 #지역더미만 포함  
 fixed<-plm(panel.model01, data=Adata, index=c('area', 'year'), model='within') 
 
 #지역더미+연도더미 포함  
 fixed.time<-plm(panel.model02, data=Adata, index=c('area', 'year'), model='within') 
 
 #plm함수 사용 후, 지역더미변수의 회귀계수 출력  
 fixef(fixed)  
   
 *** R의 plm패키지에서 각 모형을 나타내는 명령문 ***  
 결과<-plm(lwage~educ+black+south+union, data=Adata, model='pooling')  
   
 *** Panel effect가 있는가를 검사 ***   
 -- 대립가설: Pooled OLS모형이 Fixed effect model(individual)보다 더 우수하다  
 plmtest(pooled, effect='individual',type=c('bp'))  
 -- 대립가설: Pooled OLS모형이 Fixed effect model(time)보다 더 우수하다  
 plmtest(pooled, effect='time',type=c('bp'))  
   
 *** Hausman Tets ***  
 -- 대립가설: Random effect model이 Fixed effect model 보다 우수하다  
 Fixed <- plm(inv ~ value + capital, data = Grunfeld, model='within')  
 Random <- plm(inv ~ value + capital, data = Grunfeld, model='random')  
 phtest(Fixed, Random)  
   
 *** Cross-sectional dependence 점검 *** --  
 -- 대립가설: p-value가 0.05보다 작으면 NULL을 기각하고 cross-sectional dependence가 존재  
 pcdtest(fe, test=c('lm'))  
 pcdtest(fe, test=c('cd'))  
   
 *** 자기상관의 점검 ***  
 Gdata<-pdata.frame(Grunfeld, index=c('firm', 'year'))  
 Result <- plm(inv ~ value + capital, data = Gdata, model='within')  
 pdwtest(Result, alternative='two.sided')  
 pbgtest(Result)  # Breusch–Godfrey Test for Panel Models  
   
 *** unit root가 존재하면 변수들을 차분 ***
 library(plm)  
 write.csv(Grunfeld, 'Gdata.csv')  
 Gdata<-read.csv('Gdata.csv', header=T)  
 Gdata<-Gdata[,-1]  
 Gdata<-pdata.frame(Gdata, index=c('firm', 'year'))  
 uroot.inv <- data.frame(split(Gdata$inv, Gdata$firm))  
 purtest(uroot.inv , pmax = 4, exo = 'intercept', test = 'madwu')  
 # inv변수가 non-stationary이므로 차분해야 함  
   
 ** 이분산성 점검 ***
 library(plm)  
 write.csv(Grunfeld, 'Gdata.csv')  
 Gdata<-read.csv('Gdata.csv', header=T)  
 Gdata<-Gdata[,-1]  
 Gdata<-pdata.frame(Gdata, index=c('firm', 'year'))  
 bptest(inv ~ value + capital, data = Gdata, model='within')  
   
 ***  Individual and/or Time Effects 검정을 위한 F-test ***
 library(plm)  
 write.csv(Grunfeld, 'Gdata.csv')  
 Gdata<-read.csv('Gdata.csv', header=T)  
 Gdata<-Gdata[,-1]  
 Gdata<-pdata.frame(Gdata, index=c('firm', 'year'))  
 gp <- plm(inv ~ value + capital, data = Gdata, model = 'pooling')  
 gi <- plm(inv ~ value + capital, data = Gdata, effect = 'individual', model = 'within')  
 gt <- plm(inv ~ value + capital, data = Gdata, effect = 'time', model = 'within')  
 gd <- plm(inv ~ value + capital, data = Gdata, effect = 'twoways', model = 'within')  
   
 -- 대립가설: gp(Pooled)모델이 gi(Fixed모델의 Individual모형)모델보다 우수하다  
 pFtest(gi, gp)  
 
 -- 대립가설: gp(Pooled)모델이 gt(Fixed모델의 Time모형)모델보다 우수하다  
 pFtest(gt, gp)  
 
 -- 대립가설: gp(Pooled)모델이 gd모델(Fixed모델의 Individual+Teim모형)보다 우수하다  
 pFtest(gd, gp)  
   
 *** Heteroskedasticity consistent coefficients (Arellano) ***  
 vcovHC(): three heteroskedasticity-consistent covariance estimators  
 white1 - for general heteroskedasticity but no serial correlation. Recommended for random effects.  
 white2 - is 'white1' restricted to a common variance within groups. Recommended for random effects.  
 arellano - both heteroskedasticity and serial correlation. Recommended for fixed effects.  
 예) coeftest(fixed, vcovHC(fixed, method = 'arellano'))  

}