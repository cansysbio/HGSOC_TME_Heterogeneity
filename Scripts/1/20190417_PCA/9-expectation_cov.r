
rm(list=ls())

id.analysis = "20180606_PCA/1/"


# paths
if(.id=="dlc"){
    setwd("~/cruk/redmine/5218")
    path.data = "data/20180604_PCA/"
    path.data = paste0(path.data,'OVCTp_NES_Class_Clinical_data_ajs_dlc.csv')
    path.rd   = paste0("rd/",id.analysis)
}


# load import
print(load(paste0(path.rd,"0-import")))


# simluation parameters from AJS's data:
MU = apply(data_pca,2,mean)
SIGMA = cov(data_pca)
N = 5e+6
library(MASS)
DATA = mvrnorm(n = N, MU, SIGMA)
    

#
temp = c(n.patient,25,n.obs,2000)
n.n = length(temp)
id.n = .nf(data.frame(pos=1:n.n,id=temp))

# simulation
n.r = 2500
temp = 1:N
k_r = lapply(1:n.r,function(x){
    sample(temp,max(id.n$id))
    })

    
################################################################################
# covariance elements:

# object to fill:
temp = col(SIGMA)<row(SIGMA)
n.cov = sum(temp)
id.cov = .nf(data.frame(pos=1:n.cov,id=(1:prod(dim(temp)))[temp],sigma=SIGMA[temp]))
ar.estimate.rpn = array(NA,dim=c(n.r,n.cov,n.n),dimnames=list(1:n.r,id.cov$id,id.n$id))
for(rw in 1:n.r){
    for(nw in 1:n.n){
        ar.estimate.rpn[rw,,nw] = cov(DATA[k_r[[rw]][1:id.n$id[nw]],])[id.cov$id]
        }
    .cat(rw,n.r)    
    }
    
##
## plot estimated covariances
##

# small n
nw = 2
boxplot(ar.estimate.rpn[,,nw] - rep(id.cov$sigma,each=n.r),main=.p("sample size = ",id.n$id[nw]),
        ylab="resid",xlab="covariance elements")
points(apply(ar.estimate.rpn[,,nw] - rep(id.cov$sigma,each=n.r),2,mean),col="red")
abline(h=0,col="blue")
    
# large n
nw = n.n
boxplot(ar.estimate.rpn[,,nw] - rep(id.cov$sigma,each=n.r),main=.p("sample size = ",id.n$id[nw]),
        ylab="resid",xlab="covariance elements")
points(apply(ar.estimate.rpn[,,nw] - rep(id.cov$sigma,each=n.r),2,mean),col="red")
abline(h=0,col="blue")

    # CONCLUSION:
    # - covariances unbiased

    
    
################################################################################
# 2 first dim loadings:

# object to fill:
PCA = princomp(cov=SIGMA)
n.loading = n.var
id.loading = .nf(data.frame(pos=1:n.loading,id=1:n.loading,
                 value=PCA$loadings[,1]))

ar.estimate.rln = array(NA,dim=c(n.r,n.loading,n.n),dimnames=list(1:n.r,id.loading$id,id.n$id))
for(rw in 1:n.r){
    for(nw in 1:n.n){
        covw = cov(DATA[k_r[[rw]][1:id.n$id[nw]],])
        comp1w = princomp(cov=covw)$loadings[,1]
        if(cor(comp1w,id.loading$value)<0){comp1w=-comp1w}
        ar.estimate.rln[rw,,nw] = comp1w
        }
    .cat(rw,n.r)    
    }
    
##
## plot estimated covariances
##

# small n
nw = 2
boxplot(ar.estimate.rln[,,nw] - rep(id.loading$value,each=n.r),main=.p("sample size = ",id.n$id[nw]),
        ylab="resid",xlab="covariance elements")
points(apply(ar.estimate.rln[,,nw] - rep(id.loading$value,each=n.r),2,mean),col="red")
abline(h=0,col="blue")
#
boxplot(ar.estimate.rln[,,nw],main=.p("sample size = ",id.n$id[nw]),
        ylab="resid",xlab="covariance elements",col=rainbow(n.type)[.an(factor(id.var$type))])
points(apply(ar.estimate.rln[,,nw],2,mean),col="red")
points(id.loading$value,col="blue",pch=3)

# large n
nw = n.n
boxplot(ar.estimate.rln[,,nw] - rep(id.loading$value,each=n.r),main=.p("sample size = ",id.n$id[nw]),
        ylab="resid",xlab="covariance elements",col="light gray")
points(apply(ar.estimate.rln[,,nw] - rep(id.loading$value,each=n.r),2,mean),col="red")
abline(h=0,col="blue")
#
boxplot(ar.estimate.rln[,,nw],main=.p("sample size = ",id.n$id[nw]),
        ylab="resid",xlab="covariance elements",col="light gray")
points(apply(ar.estimate.rln[,,nw],2,mean),col="red")
points(id.loading$value,col="blue",pch=3)


    # CONCLUSION:
    # - covariances unbiased
    