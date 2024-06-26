#########################
###### Simulation #######
#########################

library(tea) #for mindist function
library(VGAM) #for rpareto function
library(fda)



#the first three orthogonal basis functions
c=seq(0,1,by=0.01)
v1=sqrt(2)*sin((1-0.5)*pi*c)
v2=sqrt(2)*sin((2-0.5)*pi*c)
v3=sqrt(2)*sin((3-0.5)*pi*c)

#tail index
alpha=2.1

# Obtain rho to make rhoxy=(-1, -0.9, ....., 0.9, 1.0)
rhoxy<-seq(-1,1,.1)

crho=rep(NA,21)
for (m in 1:21) {
  f<- function(x) x/sqrt(x^2+(1-x^2)^{alpha/2})-rhoxy[m]
  crho[m]=uniroot(f, c(-1, 1))$root  
}


#DGP
J = length(c)
N = 100
R = 1000
z1 <- matrix(NA, ncol = R, nrow = N)
z2 <- matrix(NA, ncol = R, nrow = N)
n1 <- matrix(NA, ncol = R, nrow = N)
n2 <- matrix(NA, ncol = R, nrow = N)
n3 <- matrix(NA, ncol = R, nrow = N)
X <-array(data=NA, dim = c(J, N, R))
Y <-array(data=NA, dim = c(J, N, R))

f_res<-matrix(NA, 21,3)
#simulation
set.seed(2024)
for (m in 1:21) {
  
  rho = crho[m]
  
  for (r in 1:R) {
    z1[,r] = (rbinom(1,1,.5)*2-1)*rpareto(N, scale = 1, shape=alpha)
    z2[,r] = (rbinom(1,1,.5)*2-1)*rpareto(N, scale = 1, shape=alpha)
    n1[,r] = rnorm(N,0,1)
    n2[,r] = rnorm(N,0,1)
    n3[,r] = rnorm(N,0,1)
  }
  
  for (r in 1:R) {
    for (i in 1:N) {
      X[,i,r] =  z1[i,r]*v1+n1[i,r]*v2+n2[i,r]*v3 
      Y[,i,r] =  rho*(z1[i,r])*v1+((1-rho^{2})^{1/2})*(z2[i,r])*v2+n3[i,r]*v3  
    }
  }
  
  res <- rep(NA,R)
  resk <- rep(NA,R)
  for (r in 1:R) {
    
    rx = apply(X[,,r], 2, function (x) sqrt(sum(x^2)))
    ry = apply(Y[,,r], 2, function (x) sqrt(sum(x^2)))
    
    rxy = apply(cbind(rx,ry), 1, max)
    rrank = length(rxy)-rank(rxy)+1
    
    # find an optimal k
    k=mindist(rxy, method = "ks")$k0
    #k=2*floor(sqrt(N))
    
    idx = which(rrank <= k)
    rk=rxy[rrank==k]
    
    # rho_xy
    res[r]=ecc(X[,,r], Y[,,r], idx)
    
    # chosen k
    resk[r]=k
    
  }
  
  f_res[m,1]=mean(res)
  f_res[m,2]=sd(res)
  f_res[m,3]=mean(resk)
  print(m)
  
}

data_5_500 = data.frame(rho= seq(-1,1,.1), Bias=round(f_res[,1]-rhoxy, 3), SE=round(f_res[,2],3), k=round(f_res[,3],1))

data=data_3_500
data=data_3_100

data=data_4_500
data=data_4_100

data=data_5_100
data=data_5_500


for (m in 1:21) {
  print(paste(data[m,1],data[m,2], "(",data[m,3],")"))
}


# extremal correlation coefficient function
ecc<-function (x, y, idx){
  temp=0
  tempx=0
  tempy=0
  for (i in 1:length(idx)) {
    temp=temp + sum(x[,idx[i]]*y[,idx[i]])
    tempx=tempx + sum(x[,idx[i]]*x[,idx[i]])
    tempy=tempy + sum(y[,idx[i]]*y[,idx[i]])
  }
  temp/sqrt(tempx*tempy)
}   

##########################
### Data Example - ETF ###
##########################

# 1=Xly
# 2=Xlp
# 3=Xle
# 4=Xlf
# 5=Xlv
# 6=Xli
# 7=Xlb
# 8=Xlk
# 9=Xlu

# Tail index estimation
for (i in 1:9) {
  
  dt=Y[i,,]-rowMeans(Y[i,,])
  
  dtx = apply(dt, 2, function (x) sqrt(sum(x^2)))
  print(round(c(i, mindist(dtx, method = "ks")$tail.index),1))
}

# Pairwise ECC 

res<-matrix(0,9,9)

for (i in 1:8) {
  
  for (j in (i+1):9) {
    
    Xdt=Y[i,,]-rowMeans(Y[i,,])
    Ydt=Y[j,,]-rowMeans(Y[j,,])
    
    N=dim(Xdt)[2]
    J=dim(Xdt)[1]
    
    # Marginal tail estimates
    Xalpha=mindist(apply(Xdt, 2, function (x) sqrt(sum(x^2))), method = "ks")$tail.index
    Yalpha=mindist(apply(Ydt, 2, function (x) sqrt(sum(x^2))), method = "ks")$tail.index
    
    new_Xdt=NA
    new_Ydt=NA
    
    Xr = apply(Xdt, 2, function (x) sqrt(sum(x^2)))
    Yr = apply(Ydt, 2, function (x) sqrt(sum(x^2)))
    
    new_Xdt=Xdt/matrix(rep(Xr^{1-Xalpha/2.5},J),J,N,byrow=TRUE)
    new_Ydt=Ydt/matrix(rep(Yr^{1-Yalpha/2.5},J),J,N,byrow=TRUE)
    # Compute rhohat_{XY}
    rx = apply(new_Xdt, 2, function (x) sqrt(sum(x^2)))
    ry = apply(new_Ydt, 2, function (x) sqrt(sum(x^2)))
    
    rxy = apply(cbind(rx,ry),1, max)
    rrank = length(rxy)-rank(rxy)+1
    
    # find an optimal k
    k=mindist(rxy, method = "ks")$k0
    
    idx = which(rrank <= k)
    rk=rxy[rrank==k]
    
    temp=0
    tempx=0
    tempy=0
    for (k in 1:length(idx)) {
      temp=temp + sum(new_Xdt[,idx[k]]*new_Ydt[,idx[k]])
      tempx=tempx + sum(new_Xdt[,idx[k]]*new_Xdt[,idx[k]])
      tempy=tempy + sum(new_Ydt[,idx[k]]*new_Ydt[,idx[k]])
    }
    
    res[i,j]<-round(temp/sqrt(tempx*tempy),2)
    
  }
  
}

colnames(res)<-c("XLY","XLP","XLE","XLF","XLV","XLI","XLB","XLK","XLU")

library(corrplot)
corrplot(res, method = 'color', type="upper", col = COL2(n=200), cl.length = 21, tl.pos = 'd', tl.col = "black",
         addCoef.col = 'white')
