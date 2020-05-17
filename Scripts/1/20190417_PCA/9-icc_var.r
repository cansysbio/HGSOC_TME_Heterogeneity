


rm(list=ls())

id.analysis = "20180606_PCA/1/"


# paths
if(.id=="dlc"){
    setwd("~/cruk/redmine/5218")
    path.data = "data/20180604_PCA/"
    path.data = paste0(path.data,'OVCTp_NES_Class_Clinical_data_ajs_dlc.csv')
    path.rd   = paste0("rd/",id.analysis)
}


# load import:
print(load(paste0(path.rd,"0-import")))


# dependence:
library(lme4)
lmer_var = apply(data_pca,2,function(x){
      lmer(x~1+(1|id.obs$patient))
      })
icc_var = sapply(lmer_var,function(x){
    sigma_i = sqrt(unlist(VarCorr(x)))
    sigma_e = attr(VarCorr(x), "sc")
    sigma_i/(sigma_e+sigma_i)
    })
plot(icc_var)
    # COMMENT: 
    # - ICC close to 0.45 for most of the variables of interest
    # - may uggest performing analysis on residuals 

    
##
## dataset on shifted resid
##

data_pca_resid = matrix(unlist(lapply(lmer_var,function(x)resid(x)+coef(summary(x))[1,1])),ncol=n.var)
dimnames(data_pca_resid) = dimnames(data_pca)
#data_pca_resid = apply(data_pca_resid,1,scale)
#
par(mfrow=c(2,2))
#
boxplot(data_pca,col="light gray")
boxplot(t(data_pca),col="light gray")
#
boxplot(data_pca_resid,col="light gray")
boxplot(t(data_pca_resid),col="light gray")
    # CONCLUSION:
    # - less variance per variable but same picture
    

##
## pca and comparison:
##

pca = princomp(covmat=cov(data_pca))
summary(pca)

pca_resid = princomp(covmat=cov(data_pca_resid))
summary(pca_resid)

# compare loadings
par(mfrow=c(1,5))
for(i in 1:5){
    plot(pca$loadings[,i],pca_resid$loadings[,i])
    cat(cor(pca$loadings[,i],pca_resid$loadings[,i]),"\t",
        cor(pca$loadings[,i+1],pca_resid$loadings[,i]),"\n")
    abline(0,1,col="blue")
}
    # CONCLUSINON:
    # - fairly similar results sugesting that patients are NOT driving effects 
    #   (as they could have) 











