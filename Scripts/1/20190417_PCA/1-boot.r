


rm(list=ls())

id.analysis = "20190417_PCA/1/"


# paths
if(.id=="dlc"){
    setwd("~/cruk/redmine/5218")
    path.rd   = paste0("rd/",id.analysis)
}


# load import
print(load(paste0(path.rd,"0-import")))


# id
n.boot  = 3
id.boot = data.frame(pos=1:n.boot,id=paste0("boot",0:(n.boot-1)),
                     name=c("Case boostrap","One-stage cluster bootstrap","two-stage cluster bootstrap"),
                     stringsAsFactors=FALSE)

              
              
##
## original estimates
##
              
n.comp = 5
original_pca = princomp(covmat=cov(data_pca))
mx.loadings_original.pc = original_pca$loadings[,1:n.comp]           



##
## bootstrap observation coordinates:
##

set.seed(478)# aj478@cam.ac.uk

# boot:
n.R = 10000
k_r_boot = as.list(rep(NA,n.boot))
names(k_r_boot) = id.boot$id

## Case boostrap
k_r_boot[[1]] = lapply(1:n.R,function(x){
    sample(id.obs$pos,n.obs,replace=TRUE)
    })
## One-stage cluster bootstrap
k_patient = split(id.obs$pos,id.obs$patient)
k_r_boot[[2]] = lapply(1:n.R,function(x){
    unlist(k_patient[sample(id.patient$pos,n.patient,replace=TRUE)])
    })
## two-stage cluster bootstrap
k_patient = split(id.obs$pos,id.obs$patient)
k_r_boot[[3]] = lapply(1:n.R,function(x){
    unlist(lapply(k_patient[sample(id.patient$pos,n.patient,replace=TRUE)],function(x)sample(x,length(x),replace=TRUE)))
    })
lapply(k_r_boot,length)

##
## PCA's loading per simulated sample and bootstrap type 
##

loadings1toc_r_boot = lapply(k_r_boot,function(x){
    lapply(x,function(y){
        princomp(covmat=cov(data_pca[y,]))$loadings[,1:n.comp]
    })
    })

ar.loadings.pcr_boot = lapply(loadings1toc_r_boot,function(x)array(unlist(x),dim=c(nrow(x[[1]]),ncol(x[[1]]),n.R),
    dimnames=list(rownames(x[[1]]),colnames(x[[2]]),1:n.R)))     
#  
pairs(t(ar.loadings.pcr_boot[[3]][1:10,1,1:1000]),col=paste0(substr(rainbow(10)[c(2,7)],1,7),50))
pairs(t(ar.loadings.pcr_boot[[3]][1:10,2,1:1000]),col=paste0(substr(rainbow(10)[c(2,7)],1,7),50))

    
##
## order the different dimensions...
##
    

# correction 1
manipulation1_r_boot = lapply(ar.loadings.pcr_boot,function(x,original){
    apply(x,3,function(y,original){# original = mx.loadings_original.pc;y = ar.loadings.pcr_boot[[1]][,,1]
    out = matrix(t(apply(y,2,function(x,original){# x=y[,1]
                       score = rbind(apply(original-x,2,function(x)sum(x^2)),
                                     apply(original+x,2,function(x)sum(x^2)))
                       score = score==min(score)
                       pos   = which(apply(score,2,sum)==1)
                       sign  = c(1,-1)[which(apply(score,1,sum)==1)]
                       c(pos,sign)},original=original)),
                  nrow=ncol(original),ncol=2,dimnames=list(paste0("comp",1:ncol(original)),c("position","sign")))
    as.data.frame(out)              
    },original = original)
    },original = mx.loadings_original.pc)
lapply(manipulation1_r_boot,function(x)table(unlist(lapply(x,function(x)length(tabulate(x[,1]))))))
    # COMMENT: duplicated dimensions frequent not for dimensions 1 and 2
 

# correction 2: NOT USED
if(FALSE){
    .coord = function(n.dk,n.tk){
        n.k=length(n.tk)
        if(n.k!=length(n.dk)){stop("coord and dimarray are different length")}
        for(k in 1:n.k){if(n.dk[k]>n.tk[k]){stop(cat("in dim", k, "n.dk > n.tk  :"))}}
        # calc
        temp=rep(n.dk[1],n.k)
        for(k in 1:(n.k-1)){temp[k]=prod(n.tk[1:(n.k-k)])*(n.dk[(n.k-k+1)]-1)}
        sum(temp)
        }
    permutation = expand.grid(1:n.comp,1:n.comp,1:n.comp,1:n.comp,1:n.comp)
    permutation = permutation[apply(permutation,1,function(x)length(unique(x))==n.comp),]
    permutation.k = permutation
    for(cw in 1:ncol(permutation)){
        permutation.k[,cw] = apply(cbind(permutation[,cw],cw),1,.coord,n.tk=c(ncol(permutation),ncol(permutation)))
        }
    manipulation2_r_boot = lapply(ar.loadings.pcr_boot,function(x,original,p,k,w){
        apply(x,3,function(y,original,p,k,w){# original = mx.loadings.pc;y = ar.loadings.pcr_boot[[1]][,,1]
            n.cw = ncol(original)
            # cor:
            mx.cor.cc = matrix(NA,n.cw,n.cw)
            for(i in 1:n.cw){
                x = y[,i]
                for(j in 1:n.cw){mx.cor.cc[i,j] = cor(x,original[,j],method="pearson")}
                }
            # solution:
            # for all permutation, pick the highest weighted correlation mean 
            # (weighted by variance of components)
            crit.perm = apply(matrix(abs(mx.cor.cc)[unlist(k)],ncol=n.cw)*matrix(rep(w,each=nrow(k)),ncol=n.cw),1,sum)
            # out
            k.best = which(crit.perm==max(crit.perm))
            out = data.frame(position=unlist(p[k.best,]),
                             sign = sign(mx.cor.cc[unlist(k[k.best,])]),
                             corr = mx.cor.cc[unlist(k[k.best,])])
            out                                                                                                      
            },original = original,p=p,k=k,w=w)
            },original = mx.loadings_original.pc,p=permutation,k=permutation.k,w=original_pca$sdev[1:ncol(permutation)]^1/sum(original_pca$sdev[1:ncol(permutation)]^1))
    }   
              
# use correction 1 and then delete outliers
ar.loadings_corr.pcr_boot = ar.loadings.pcr_boot
for(bw in 1:n.boot){
    for(rw in 1:n.R){
        ar.loadings_corr.pcr_boot[[bw]][,,rw] = ar.loadings.pcr_boot[[bw]][,manipulation1_r_boot[[bw]][[rw]]$position,rw]*rep(manipulation1_r_boot[[bw]][[rw]]$sign[manipulation1_r_boot[[bw]][[rw]]$position],each=n.var) 
        }
}              
# visualise
pairs(t(ar.loadings_corr.pcr_boot[[3]][1:10,1,1:1000]),col=paste0(gray(.5),50))
pairs(t(ar.loadings_corr.pcr_boot[[3]][1:10,2,1:1000]),col=paste0(gray(.5),50))
    # COMMENT:
    # -> way cleaner but some outliers remain -> detect and remove in senstitivity analysis
    
 
# outlier detection via robust covariance matrix estimator per boostrap type:
path.sc = "~/robustness/robustbiostat/sc/20160330/"
source(paste(path.sc,"covREM.r",sep=""))
source(paste(path.sc,"AttrPsyFunc.r",sep=""))
notoutliers.r_boot = lapply(ar.loadings_corr.pcr_boot,function(x){
    apply(apply(x[,1:2,],2,function(x){covREM(apply(t(x),2,scale))$w}),
          1,sum)>1
})
#
path.sc = "~/robustness/robustbiostat/sc/20160330/"
source(paste(path.sc,"covREM.r",sep=""))
source(paste(path.sc,"AttrPsyFunc.r",sep=""))
notoutliers.r_boot = lapply(ar.loadings_corr.pcr_boot,function(x){
    apply(apply(x[,1:2,],2,function(x){covREM(apply(t(x),2,scale))$w>0.5}),
          1,all)
})
    ## by hand:
    # x  = ar.loadings_corr.pcr_boot[[3]]
    # d1 = covREM(apply(t(x[,1,]),2,scale))
    # d2 = covREM(apply(t(x[,2,]),2,scale))
    # plot(d1$w,d2$w,col=paste0(substr(rainbow(10)[c(2,7)],1,7),50)[(d1$w>0.25&d2$w>0.25)+1])
    # notoutliers.r_boot[[3]] = d1$w>0.5&d2$w>0.5

    # COMMENT:
    # -> input of s estimator singular (mcd) but final ok
        
## check manipulation results: for dim1 and 2, compare ORIGINAL loadings versus CORRECTED ONES
pairs(t(ar.loadings_corr.pcr_boot[[3]][1:10,1,1:1000]),col=paste0(substr(rainbow(10)[c(2,7)],1,7),50)[as.numeric(notoutliers.r_boot[[3]])+1])
pairs(t(ar.loadings_corr.pcr_boot[[3]][11:20,1,1:1000]),col=paste0(substr(rainbow(10)[c(2,7)],1,7),50)[as.numeric(notoutliers.r_boot[[3]])+1])
pairs(t(ar.loadings_corr.pcr_boot[[3]][1:10,2,1:1000]),col=paste0(substr(rainbow(10)[c(2,7)],1,7),50)[as.numeric(notoutliers.r_boot[[3]])+1])
    # COMMENT:
    # - seems ok
    
    
save(id.boot,n.boot,ar.loadings_corr.pcr_boot,mx.loadings_original.pc,notoutliers.r_boot,n.R,file=paste0(path.rd,"1-boot-noscaling"))    
    

