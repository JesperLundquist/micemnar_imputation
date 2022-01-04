#### Bortfall 10%####
#library
library(MASS)
library(mice)
library(miceMNAR)
library("GJRM")
library(mvtnorm)
library(tidyverse)
library(plyr)
#Inställningar
set.seed(198) #För att kunna ta fram samma randomnummer igen
ni = 2000 #2000 observationer
mean.x=0
sd.x=0.3
rho0=0
rho03=0.3
rho06=0.6
rho09=0.9
mu.e=c(0,0)
b0=0
b1=1
b2=1
b0s=1.4
b1s=1
b2s=-0.5
b3s=1
true=1
sim = 1000     # Antal dataset som skapas

## Storage, sparar alla skattningar i matriser 
missing <- matrix(NA,sim,5) #andel bortfall
betaestimate <- matrix(NA,sim,25) #betaskattningar
betastd <- matrix(NA,sim,25) #beta standardfel
CIl <- matrix(NA,sim,25) #Confidens interval low
CIh <- matrix(NA,sim,25) #Confidens interval High

# Tar tiden
ptm <- proc.time()

for (i in 1:sim){ #  sim är antal iterationer
  #Generera feltermer rho=0
  sigma.e0 <- matrix(c(1 ,  rho0,
                       rho0 , 1),2,2)
  df0 <- mvrnorm(n=ni,mu=mu.e,Sigma=sigma.e0)  #simulering av data
  colnames(df0) <- c("e0","e.s0")
  
  #Generera feltermer rho=0.3
  sigma.e03 <- matrix(c(1 ,  rho03,
                        rho03 , 1),2,2)
  df03 <- mvrnorm(n=ni,mu=mu.e,Sigma=sigma.e03)  #simulera datan
  colnames(df03) <- c("e03","e.s03")
  
  #Generera feltermer rho=0.6
  sigma.e06 <- matrix(c(1 ,  rho06,
                        rho06 , 1),2,2)
  df06 <- mvrnorm(n=ni,mu=mu.e,Sigma=sigma.e06)  #simulera datan
  colnames(df06) <- c("e06","e.s06")
  
  #Generera feltermer rho=0.9
  sigma.e09 <- matrix(c(1 ,  rho09,
                        rho09 , 1),2,2)
  df09 <- mvrnorm(n=ni,mu=mu.e,Sigma=sigma.e09)  #simulera datan
  colnames(df09) <- c("e09","e.s09")
  
  #Generera variabler
  x1 <- rnorm(ni, mean = mean.x, sd=sd.x)
  x2 <- rnorm(ni, mean = mean.x, sd=sd.x) 
  x3 <- rnorm(ni, mean = mean.x, sd=sd.x)
  x4 <- rnorm(ni, mean = mean.x, sd=sd.x)
  e2 <- rnorm(ni, mean= mean.x, sd=1)
  
  data0 <- as.data.frame(cbind(x1,x2,x3,e2, df0))
  data03 <- as.data.frame(cbind(x1,x2,x3, df03))
  data06 <- as.data.frame(cbind(x1,x2,x3, df06))
  data09 <- as.data.frame(cbind(x1,x2,x3, df09))
  datagan <- as.data.frame(cbind(x1,x2,e2,x3, df0))
  
  ###Generera y för varje bortfallsmekanism (rho=0, rho=0.3, rho=0.6, rho=0.9, “Ganjali”)
  
  data0$y <- b0 + b1*data0$x1 + b2*data0$x2 + data0$e0
  data03$y <- b0 + b1*data03$x1 + b2*data03$x2 + data03$e03
  data06$y <- b0 + b1*data06$x1 + b2*data06$x2 + data06$e06
  data09$y <- b0 + b1*data09$x1 + b2*data09$x2 + data09$e09
  
  datagan$y <- b0 + b1*data0$x1 + b2*data0$x2 + data0$e2
  
  ###skattningar innan bortfall###
  est0<-lm (y~x1+x2, data=data0) 
  betaestimate[i,1] <- coef(summary(est0))["x1","Estimate"] # rho=0 medelvärde för b1hat
  betastd[i,1] <- coef(summary(est0))["x1","Std. Error"] # rho=0 std error.
  CIl[i,1]<- confint(est0, 'x1', level=0.95)[1] # rho=0 CI low 
  CIh[i,1]<- confint(est0, 'x1', level=0.95)[2] # rho=0 CI high
  
  est03<-lm (y~x1+x2, data=data03) 
  betaestimate[i,2] <- coef(summary(est03))["x1","Estimate"] # rho=0.3 medelvärde för b1hat
  betastd[i,2] <- coef(summary(est03))["x1","Std. Error"] # rho=0.3 std error.
  CIl[i,2]<- confint(est03, 'x1', level=0.95)[1] # rho=0.3 CI low 
  CIh[i,2]<- confint(est03, 'x1', level=0.95)[2] # rho=0.3 CI high
  
  est06<-lm (y~x1+x2, data=data06) 
  betaestimate[i,3] <- coef(summary(est06))["x1","Estimate"] # rho=0.6 medelvärde för b1hat
  betastd[i,3] <- coef(summary(est06))["x1","Std. Error"] # rho=0.6 std error.
  CIl[i,3]<- confint(est06, 'x1', level=0.95)[1] # rho=0.6 CI low 
  CIh[i,3]<- confint(est06, 'x1', level=0.95)[2] # rho=0.6 CI high
  
  est09<-lm (y~x1+x2, data=data09) 
  betaestimate[i,4] <- coef(summary(est09))["x1","Estimate"] # rho=0.9 medelvärde för b1hat
  betastd[i,4] <- coef(summary(est09))["x1","Std. Error"] # rho=0.9 std error.
  CIl[i,4]<- confint(est09, 'x1', level=0.95)[1] # rho=0.9 CI low 
  CIh[i,4]<- confint(est09, 'x1', level=0.95)[2] # rho=0.9 CI high
  
  estganjali<-lm (y~x1+x2, data=datagan) 
  betaestimate[i,5] <- coef(summary(estganjali))["x1","Estimate"] #ganjali medelvärde b1hat
  betastd[i,5] <- coef(summary(estganjali))["x1","Std. Error"] # ganjali std error.
  CIl[i,5]<- confint(estganjali, 'x1', level=0.95)[1] # ganjali CI low 
  CIh[i,5]<- confint(estganjali, 'x1', level=0.95)[2] #ganjali CI high
  
  ####Generera Bortfall#####
  #Bortfall rho=0
  #R=0 betyder bortfall
  data0$R <- ifelse(b0s + b1s*data0$x1 + b2s*data0$x2 + b3s*data0$x3 + data0$e.s < 0, 0, 1) #skapar bortfall enligt Heckman. värden på b från Galimard
  missing[i,1] <- sum(data0$R==0)/ni #kontrollera storlek på bortfallet
  yMNAR0 <- replace(data0$y, data0$R==0, NA)
  data0 <- as.data.frame(cbind(x1,x2,x3, yMNAR0))
  
  #Bortfall rho=0.3
  #R=0 betyder bortfall
  data03$R <- ifelse(b0s + b1s*data03$x1 + b2s*data03$x2 + b3s*data03$x3 + data03$e.s < 0, 0, 1) #skapar bortfall enligt Heckman. värden på b från Galimard
  missing[i,2] <- sum(data03$R==0)/ni #kontrollera storlek på bortfallet
  yMNAR03 <- replace(data03$y, data03$R==0, NA)
  data03 <- as.data.frame(cbind(x1,x2,x3, yMNAR03))
  
  #Bortfall rho=0.6
  #R=0 betyder bortfall
  data06$R <- ifelse(b0s + b1s*data06$x1 + b2s*data06$x2 + b3s*data06$x3 + data06$e.s < 0, 0, 1) #skapar bortfall enligt Heckman. värden på b från Galimard
  missing[i,3] <- sum(data06$R==0)/ni #kontrollera storlek på bortfallet
  yMNAR06 <- replace(data06$y, data06$R==0, NA)
  data06 <- as.data.frame(cbind(x1,x2,x3, yMNAR06))
  
  #Bortfall rho=0.9
  #R=0 betyder bortfall
  data09$R <- ifelse(b0s + b1s*data09$x1 + b2s*data09$x2 + b3s*data09$x3 + data09$e.s < 0, 0, 1) #skapar bortfall enligt Heckman. värden på b från Galimard
  missing[i,4] <-  sum(data09$R==0)/ni #kontrollera storlek på bortfallet
  yMNAR09 <- replace(data09$y, data09$R==0, NA)
  data09 <- as.data.frame(cbind(x1,x2,x3, yMNAR09))
  
  #Bortfall ganjali
  #kumulativ normalfördelning genererar Pr(Z<y+1.9)
  datagan$ganjali <- pnorm(1.9+datagan$y) 
  #skapar en binär variabel med sannolikheten för 1 = ganjali
  datagan$rber <- rbinom(ni,size=1, datagan$ganjali) 
  missing[i,5] <- sum(datagan$rber==0)/ni
  #visar y-variabeln i form av det andra bortfallet
  yMNARgan <- replace(datagan$y, datagan$rber==0, NA) 
  datagan <- as.data.frame(cbind(x1,x2,x3, yMNARgan))
  
  #CCA -Complete case analysis###
  #rho=0
  cca<-lm (yMNAR0~x1+x2, data=data0) 
  betaestimate[i,6] <- coef(summary(cca))["x1","Estimate"] # CCA medelvärde b1hat 
  betastd[i,6] <- coef(summary(cca))["x1","Std. Error"] # CCA std error.
  CIl[i,6]<- confint(cca, 'x1', level=0.95)[1]
  CIh[i,6]<- confint(cca, 'x1', level=0.95)[2]
  
  #rho=0.3
  cca<-lm (yMNAR03~x1+x2, data=data03) 
  betaestimate[i,7] <- coef(summary(cca))["x1","Estimate"] # CC medelvärde b1hat 
  betastd[i,7] <- coef(summary(cca))["x1","Std. Error"] # CCA std error.
  CIl[i,7]<- confint(cca, 'x1', level=0.95)[1]
  CIh[i,7]<- confint(cca, 'x1', level=0.95)[2]
  
  #rho=0.6  
  cca<-lm (yMNAR06~x1+x2, data=data06) 
  betaestimate[i,8] <- coef(summary(cca))["x1","Estimate"] # CCA medelvärde b1hat
  betastd[i,8] <- coef(summary(cca))["x1","Std. Error"] # CCA std error.
  CIl[i,8]<- confint(cca, 'x1', level=0.95)[1]
  CIh[i,8]<- confint(cca, 'x1', level=0.95)[2]
  
  #rho=0.9
  cca<-lm (yMNAR09~x1+x2, data=data09) 
  betaestimate[i,9] <- coef(summary(cca))["x1","Estimate"] # CCA medelvärde b1hat 
  betastd[i,9] <- coef(summary(cca))["x1","Std. Error"] # CCA std error.
  CIl[i,9]<- confint(cca, 'x1', level=0.95)[1]
  CIh[i,9]<- confint(cca, 'x1', level=0.95)[2]
  
  #Ganjali
  cca<-lm (yMNARgan~x1+x2, data=datagan) 
  betaestimate[i,10] <- coef(summary(cca))["x1","Estimate"] # CCA medelvärde b1hat
  betastd[i,10] <- coef(summary(cca))["x1","Std. Error"] # CCA std error.
  CIl[i,10]<- confint(cca, 'x1', level=0.95)[1]
  CIh[i,10]<- confint(cca, 'x1', level=0.95)[2]
  
  ###Imputering MInorm, stokastisk regressionsmodell### 
  #rho=0
  imp <- mice(data0, method = "norm",maxit=1, m = 10, print = FALSE)
  fit <- with(imp, lm(yMNAR0 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,11] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,11] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,11]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,11]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.3
  imp <- mice(data03, method = "norm",maxit=1, m = 10, print = FALSE)
  fit <- with(imp, lm(yMNAR03 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,12] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,12] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,12]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,12]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.6  
  imp <- mice(data06, method = "norm",maxit=1, m = 10, print = FALSE)
  fit <- with(imp, lm(yMNAR06 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,13] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,13] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,13]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,13]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.9 
  imp <- mice(data09, method = "norm",maxit=1, m = 10, print = FALSE)
  fit <- with(imp, lm(yMNAR09 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,14] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,14] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,14]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,14]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #Ganjali 
  imp <- mice(datagan, method = "norm",maxit=1, m = 10, print = FALSE)
  fit <- with(imp, lm(yMNARgan ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,15] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,15] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,15]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,15]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  ###Imputering MIHE2, Multipel imputering med Heckmans 2-stegsmodell###
  
  #rho=0
  JointModelEq <- generate_JointModelEq(data=data0,varMNAR = "yMNAR0")
  JointModelEq[,"yMNAR0_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR0_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data0,varMNAR="yMNAR0",JointModelEq=JointModelEq)
  arg$method["yMNAR0"] <- "hecknorm2step"
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR0 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,16] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,16] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,16]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,16]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.3
  JointModelEq <- generate_JointModelEq(data=data03,varMNAR = "yMNAR03")
  JointModelEq[,"yMNAR03_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR03_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data03,varMNAR="yMNAR03",JointModelEq=JointModelEq)
  arg$method["yMNAR03"] <- "hecknorm2step"
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR03 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,17] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,17] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,17]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,17]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.6 
  JointModelEq <- generate_JointModelEq(data=data06,varMNAR = "yMNAR06")
  JointModelEq[,"yMNAR06_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR06_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data06,varMNAR="yMNAR06",JointModelEq=JointModelEq)
  arg$method["yMNAR06"] <- "hecknorm2step"
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR06 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,18] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,18] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,18]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,18]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.9
  JointModelEq <- generate_JointModelEq(data=data09,varMNAR = "yMNAR09")
  JointModelEq[,"yMNAR09_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR09_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data09,varMNAR="yMNAR09",JointModelEq=JointModelEq)
  arg$method["yMNAR09"] <- "hecknorm2step"
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR09 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,19] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,19] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,19]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,19]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #Ganjali 
  JointModelEq <- generate_JointModelEq(data=datagan,varMNAR = "yMNARgan")
  JointModelEq[,"yMNARgan_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNARgan_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=datagan,varMNAR="yMNARgan",JointModelEq=JointModelEq)
  arg$method["yMNARgan"] <- "hecknorm2step"
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNARgan ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,20] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,20] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,20]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,20]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  ###Imputering MIHEml, Multipel imputering med Heckmans 1-stegsmodell/maximum likelihood###
  
  #rho=0
  JointModelEq <- generate_JointModelEq(data=data0,varMNAR = "yMNAR0")
  JointModelEq[,"yMNAR0_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR0_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data0,varMNAR="yMNAR0",JointModelEq=JointModelEq)
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR0 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,21] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,21] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,21]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,21]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.3
  JointModelEq <- generate_JointModelEq(data=data03,varMNAR = "yMNAR03")
  JointModelEq[,"yMNAR03_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR03_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data03,varMNAR="yMNAR03",JointModelEq=JointModelEq)
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR03 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,22] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,22] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,22]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,22]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.6
  JointModelEq <- generate_JointModelEq(data=data06,varMNAR = "yMNAR06")
  JointModelEq[,"yMNAR06_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR06_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data06,varMNAR="yMNAR06",JointModelEq=JointModelEq)
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR06 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,23] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,23] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,23]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,23]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #rho=0.9  
  JointModelEq <- generate_JointModelEq(data=data09,varMNAR = "yMNAR09")
  JointModelEq[,"yMNAR09_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNAR09_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=data09,varMNAR="yMNAR09",JointModelEq=JointModelEq)
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNAR09 ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,24] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,24] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,24]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,24]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
  
  #Ganjali
  JointModelEq <- generate_JointModelEq(data=datagan,varMNAR = "yMNARgan")
  JointModelEq[,"yMNARgan_var_sel"] <- c(1,1,1,0)
  JointModelEq[,"yMNARgan_var_out"] <- c(1,1,0,0)
  arg <- MNARargument(data=datagan,varMNAR="yMNARgan",JointModelEq=JointModelEq)
  imp <- mice(data = arg$data_mod,
              method = arg$method,
              predictorMatrix = arg$predictorMatrix,
              JointModelEq=arg$JointModelEq,
              control=arg$control, maxit=1, m=10, print = FALSE)
  fit <- with(imp, lm(yMNARgan ~ x1 + x2 ))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  betaestimate[i,25] <- as.numeric(tab[tab$term == "x1", "estimate"])
  betastd[i,25] <- as.numeric(tab[tab$term == "x1", "std.error"])
  CIl[i,25]<- as.numeric(tab[tab$term == "x1", "2.5 %"])
  CIh[i,25]<- as.numeric(tab[tab$term == "x1", "97.5 %"])
}
#Stop the loop
proc.time() - ptm # Stoppa tidtagning


#MÅTT#
beta <-colMeans(betaestimate) #skattning av beta
PB <- 100 * ((colMeans(betaestimate) - true)/ true) #percent bias
CR <- colMeans((CIl) < true & true < (CIh)) #Coverage rate
AW <- colMeans((CIh) - (CIl)) # “Average width”, Medelbredden på konfidensintervall 
St.e <- colMeans(betastd) #Standardfel
RMSE <- sqrt(colMeans(((betaestimate) - true)^2)) #Root mean square error


#Namnge metod och bortfallsmekanism
beta <- 
  as.data.frame(beta)%>%
  mutate(rho=rep(c("0","0.3","0.6","0.9","Non-Heckman"),5))%>%
  mutate(Metod=rep(c("Innan bortfall", "CCA", "MInorm", "MIHE2", "MIHEml"), each=5))
PB <- 
  as.data.frame(PB)%>%
  mutate(rho=rep(c("0","0.3","0.6","0.9","Non-Heckman"),5))%>%
  mutate(Metod=rep(c("Innan bortfall", "CCA", "MInorm", "MIHE2", "MIHEml"), each=5))
CR <- 
  as.data.frame(CR)%>%
  mutate(rho=rep(c("0","0.3","0.6","0.9","Non-Heckman"),5))%>%
  mutate(Metod=rep(c("Innan bortfall", "CCA", "MInorm", "MIHE2", "MIHEml"), each=5))
AW <- 
  as.data.frame(AW)%>%
  mutate(rho=rep(c("0","0.3","0.6","0.9","Non-Heckman"),5))%>%
  mutate(Metod=rep(c("Innan bortfall", "CCA", "MInorm", "MIHE2", "MIHEml"), each=5))
St.e <- 
  as.data.frame(St.e)%>%
  mutate(rho=rep(c("0","0.3","0.6","0.9","Non-Heckman"),5))%>%
  mutate(Metod=rep(c("Innan bortfall", "CCA", "MInorm", "MIHE2", "MIHEml"), each=5))
RMSE <- 
  as.data.frame(RMSE)%>%
  mutate(rho=rep(c("0","0.3","0.6","0.9","Non-Heckman"),5))%>%
  mutate(Metod=rep(c("Innan bortfall", "CCA", "MInorm", "MIHE2", "MIHEml"), each=5))

#Sammanställ resultat
result <- join_all(list(beta, PB, CR, AW, St.e, RMSE), by=c("rho","Metod"), type='left')

#Skapa tabell, avrunda, sätt i rätt ordning
tabell10<- result[,c("Metod","rho", "beta", "PB", "AW", "St.e", "RMSE","CR")]%>%
  mutate_if(is.numeric, round, digits=4)%>%
  mutate(CR=round(CR*100,1),
         PB=round(PB,1))%>%
  arrange(match(rho, c("0","0.3","0.6", "0.9", "Non-Heckman")))%>%
  arrange(match(Metod, c("Innan bortfall", "CCA","MInorm","MIHE2","MIHEml")))


#Spara tabell
write.table(tabell10, file = "tabell10.txt", sep = ",", quote = FALSE, row.names = F)

#Sparar resultat i csv-filer#
write.csv(tabell10, file="tabell10.csv")

