gfevd = function(model, n.ahead=10,normalize=TRUE,standardize=TRUE) {
  if (class(model) != "varest") {
    return("The model class is not varest!")
  }
  A <- Phi(model, (n.ahead-1))
  epsilon <- residuals(model)
  Sigma <- t(epsilon)%*%epsilon / (model$obs)
  gi <- array(0, dim(A))
  sigmas <- sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[,,j] <- t(A[,,j]%*%Sigma%*%solve(diag(sqrt(diag(Sigma)))))
  }
  
  if (standardize==TRUE){
    girf=array(NA, c(dim(gi)[1],dim(gi)[2], (dim(gi)[3])))
    for (i in 1:dim(gi)[3]){
      girf[,,i]=((gi[,,i])%*%solve(diag(diag(gi[,,1]))))
    } 
    gi = girf
  }
  
  num <- apply(gi^2,1:2,sum)
  den <- c(apply(num,1,sum))
  fevd <- t(num)/den
  if (normalize==TRUE) {
    fevd=(fevd/apply(fevd, 1, sum))
  } else {
    fevd=(fevd)
  }
  return = list(fevd=fevd, girf=gi)
} 

VS = function(CV){
  k = dim(CV)[1]
  SOFM = apply(CV,1:2,mean)*100
  VSI = (sum(rowSums(SOFM-diag(diag(SOFM))))/k)
  INC = colSums(SOFM)
  TO = colSums(SOFM-diag(diag(SOFM)))
  FROM = rowSums(SOFM-diag(diag(SOFM)))
  NET = TO-FROM
  
  ALL = rbind(rbind(rbind(cbind(SOFM,FROM),c(TO,sum(TO))),c(INC,NA)),c(NET,VSI))
  
  colnames(ALL) = c(rownames(CV),"FROM")
  rownames(ALL) = c(rownames(CV),"Directional TO Others","Directional Inlcuding Own","NET Directional Connectedness")
  ALL = format(round(ALL,2),nsmall=2)
  ALL[nrow(ALL)-1,ncol(ALL)] = "TCI"
  return = list(SOFM=SOFM,VSI=VSI,TO=TO,FROM=FROM,NET=NET,ALL=ALL,NPSO=INC)  
} 

####### Read the data
DATA = read.csv('/Users/tracy/文稿/SUIBE/2022大创/代码/对数收益率_百分数.csv')
DATE = as.Date(DATA$X)
Y = DATA[,-1]
A = Y[,6]
Y = Y[,1:5]
k = ncol(Y)

#######1 Description data
library(moments,tseries)
library(fUnitRoots)
library(fBasics)
ds1 <- rbind(apply(Y,2,mean),
               apply(Y,2,median),
               apply(Y,2,max),
               apply(Y,2,min),
               apply(Y,2,sd),
               apply(Y,2,skewness),
               apply(Y,2,kurtosis))
rownames(ds1) <- c("Mean","Median","Maximum","Minimum","Std.Deviation","Skewness","Kurtosis")
apply(Y,2,adfTest)
apply(Y,2,jarqueberaTest)
apply(Y,2,PP.test)



#######2 static spillover table
library(vars)
nlag <- 4   
nfore <- 10 
var_full <- VAR(Y, p=nlag, type="const")
CV_full <- gfevd(var_full, n.ahead=nfore)$fevd
rownames(CV_full)=colnames(CV_full)=colnames(Y)
TABLE2 = VS(CV_full)
TABLE2

########3 dynamic spillover table
t <- nrow(Y) 
space <- 100 
CV = array(NA, c(k, k, (t-space)))
colnames(CV) = rownames(CV) = colnames(Y)
for (i in 1:dim(CV)[3]) {
  var1 = VAR(Y[i:(space+i-1),], p=nlag, type="const")
  if(any(roots(var1)>1) & i>1){ 
    CV[,,i] = CV[,,(i-1)] 
  } else {
    CV[,,i] = gfevd(var1, n.ahead=nfore)$fevd
  }
  if (i%%500==0) {print(i)}
}

to = matrix(NA, ncol=k, nrow=(t-space))
from = matrix(NA, ncol=k, nrow=(t-space))
net = matrix(NA, ncol=k, nrow=(t-space))
total = matrix(NA,ncol=1,nrow=(t-space))
for (i in 1:dim(CV)[3]){
  vd = VS(CV[,,i])
  to[i,] = vd$TO/k
  from[i,] = vd$FROM/k
  net[i,] = vd$NET/k
  total[i,] = vd$VSI
}

########4 data visualization
date = DATE[-c(1:space)]
par(mfrow = c(1,1), oma = c(0,0,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
plot(date,total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=0.01)
polygon(c(date,rev(date)),c(c(rep(0,nrow(total))),rev(total)),col="grey20", border="grey20")
box()


###  Net Soillovers (ESG-Energy/ESG-Market)
t = 4
par(mfrow = c(2,2), oma = c(0,0,0,0) + 0.1, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))

#(a)*4 lightblue_esg, steelblue_energy
for (i in 1:t){
  plot(date,net[,5], xlab="",ylab="",type="l",xaxs="i",col="steelblue", las=1, main=paste(colnames(Y)[i], "/Energy"),ylim=c(floor(min(net)),ceiling(max(net))),tck=0.01,yaxs="i")
  polygon(c(date,rev(date)),c(c(rep(0,nrow(net))),rev(net[,5])),col="steelblue", border="steelblue")
  par(new=TRUE)
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="lightblue", las=1, main=paste(colnames(Y)[i], "/Energy"),ylim=c(floor(min(net)),ceiling(max(net))),tck=0.01,yaxs="i")
  polygon(c(date,rev(date)),c(c(rep(0,nrow(net))),rev(net[,i])),col="lightblue", border="lightblue")
  legend("topleft", inset=.00, c("ESG", "Energy"), col = c("lightblue", "steelblue" ), lty=c(1, 1),lwd=3)
  box()
}
box()

#(b)*4 lightblue_esg, grey40_market
for (i in 1:t){
  plot(date,net[,6], xlab="",ylab="",type="l",xaxs="i",col="grey40", las=1, main=paste(colnames(Y)[i], "/Market"),ylim=c(floor(min(net)),ceiling(max(net))),tck=0.01,yaxs="i")
  polygon(c(date,rev(date)),c(c(rep(0,nrow(net))),rev(net[,6])),col="grey40", border="grey40")
  par(new=TRUE)
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="lightblue", las=1, main=paste(colnames(Y)[i], "/Market"),ylim=c(floor(min(net)),ceiling(max(net))),tck=0.01,yaxs="i")
  polygon(c(date,rev(date)),c(c(rep(0,nrow(net))),rev(net[,i])),col="lightblue", border="lightblue")
  legend("topleft", inset=.00, c("ESG", "Market"), col = c("lightblue", "grey40" ), lty=c(1, 1),lwd=3)
  box()
}
box()

## Total Spillover
t <- nrow(Y) 
space1 <- 200
date1 <- DATE[-c(1:space1)]
CV1 = array(NA, c(k, k, (t-space1)))
colnames(CV1) = rownames(CV1) = colnames(Y)
for (i in 1:dim(CV1)[3]) {
  var1 = VAR(Y[i:(space1+i-1),], p=nlag, type="const")
  if(any(roots(var1)>1) & i>1){ 
    CV1[,,i] = CV1[,,(i-1)]
  } else {
    CV1[,,i] = gfevd(var1, n.ahead=nfore)$fevd
  }
  if (i%%500==0) {print(i)}
}

total1 = matrix(NA,ncol=1,nrow=(t-space1))
for (i in 1:dim(CV1)[3]){
  vd1 = VS(CV1[,,i])
  total1[i,] = vd1$VSI
}

win.graph()
par(mfrow  = c(1,1),oma = c(0,0,0,0) + 0.05, mar = c(8,1,8,1), mgp = c(0, 0.1, 0))
plot(date,total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=0.02, cex.text=0.5)
par(new=TRUE)
plot(date1,total1, type="l",xaxs="i",col="red", las=1, main="",ylab="",ylim=c(floor(min(total1)),ceiling(max(total1))),xlab="",tck=0.02,xaxt="n",yaxt="n")
legend("bottomleft", inset=.00, c("H=100", "H=200"), col=c("grey20", "red"), lty=c(1, 1), lwd=2, cex=0.6)


########5 ADCC-GARCH to obtain dynamic correlation
library(rmgarch)
library(rugarch)

uspec1 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1, 1), submodel=NULL, external.regressors=NULL), mean.model = list(armaOrder=c(1, 1), include.mean=TRUE, external.regressors=NULL), distribution.model = "snorm")
uspec2 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1, 1), submodel=NULL, external.regressors=NULL), mean.model = list(armaOrder=c(1, 1), include.mean=TRUE, external.regressors=NULL), distribution.model = "snorm")
uspec3 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1, 1), submodel=NULL, external.regressors=NULL), mean.model = list(armaOrder=c(1, 1), include.mean=TRUE, external.regressors=NULL), distribution.model = "snorm")
uspec4 <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1, 1), submodel=NULL, external.regressors=NULL), mean.model = list(armaOrder=c(1, 1), include.mean=TRUE, external.regressors=NULL), distribution.model = "snorm")
uspec5_energy <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1, 1), submodel=NULL, external.regressors=NULL), mean.model = list(armaOrder=c(1, 1), include.mean=TRUE, external.regressors=NULL), distribution.model = "snorm")
uspec6_market <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1, 1), submodel=NULL, external.regressors=NULL), mean.model = list(armaOrder=c(1, 1), include.mean=TRUE, external.regressors=NULL), distribution.model = "snorm")
mspec6 <- multispec(c(uspec1, uspec2, uspec3, uspec4, uspec5_energy, uspec6_market))

dcc_spec <- dccspec(mspec6, dccOrder=c(1, 1), model = "aDCC", distribution = "mvnorm")

fdcc <- dccfit(dcc_spec, Y)
fdcc

raw_dcov = rcov(fdcc)

########6 Hedge Effectiveness
t1 = dim(raw_dcov)[3]

dcc_cov <- matrix(NA, nrow=6, ncol=3)
rownames(dcc_cov) = colnames(Y)
colnames(dcc_cov) = c("h_self", "h_esg/energy","h_esg/market")
for (i in 1:t1){
  dcc_esg <- round(mean(raw_dcov[i, i, ]), digits = 3)
  dcc_esg2energy <- round(mean(raw_dcov[i, 5, ]), digits = 3)
  dcc_esg2market <- round(mean(raw_dcov[i, 6, ]), digits = 3)
  dcc_cov[i, ] <- c(dcc_esg, dcc_esg2energy, dcc_esg2market)}

esg_h = matrix(NA, nrow=4, ncol=5)
colnames(esg_h) = c("w_esg2energy", "w_esg2market","HR_esg2energy", "HR_esg2market","HE")

### Calculate Hedge Ratio (HR)
H_energy = dcc_cov[5,1]
H_market = dcc_cov[6,1]

#HR esg2energy
for (i in 1:4){
  h_i = dcc_cov[i,1]
  h_m = dcc_cov[i,2]
  w_i= (H_energy-h_m)/(h_i-2*h_m+H_energy)
  hr = h_m / h_i
  
  esg_h[i,1] <- w_i
  esg_h[i,3] <- hr
}
#HR esg2market
for (i in 1:4){
  h_i = dcc_cov[i,1]
  h_m = dcc_cov[i,3]
  w_i= (H_market-h_m)/(h_i-2*h_m+H_market)
  hr = h_m / h_i
  
  esg_h[i,2] <- w_i
  esg_h[i,4] <- hr
}

### calculate Hedge Effectiveness(HE)
for (i in 1:4){
  var_uh = var(Y[,i])
  var_h = var((Y[,i]+Y[,5])*esg_h[,i])
  HE = 1- var_h/var_uh
  esg_h[i,5] <- HE
}

esg_h = round(esg_h, digits=3)
esg_h


par(mfrow = c(2,2), oma = c(0,0,0,0) + 0.1, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:5){
  plot(DATE,Y[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1,ylim=c(-100,100),tck=0.01,yaxs="i",main=colnames(Y)[i])
}
plot(DATE,A, xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1,ylim=c(-100,100),tck=0.01,yaxs="i",main=colnames(A))

### END

