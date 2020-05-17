

rm(list=ls())

id.analysis  = "20190417_PCA/1/"

##
## paths
##
if(.id=="dlc"){
    setwd("~/cruk/redmine/5218")
    path.data = paste0("data/20190417_PCA/OVCTp_NES_Class_Clinical_data_ajs_REBUTTAL.txt")
    path.id   = paste0("data/20190417_PCA/TreatmentNaive_SampleLabels_WESTumourCellularity.txt")
    path.rd   = paste0("rd/",id.analysis)
}




##
## import
##

data0 <- read.table(path.data,row.names='Term',header=T,
                    sep='\t',stringsAsFactors=FALSE)
col0 = colnames(data0)[-ncol(data0)]
row0 = rownames(data0)[-(nrow(data0)-0:1)]




##
## PCA with and without scaling:
## comparison of the loadings
##


# PCA with scaling by observations/patients
data_pca_scale = t(apply(apply(data0[row0,col0],2,as.numeric),2,scale))
colnames(data_pca_scale) = row0

# PCA without scaling 
data_pca_noscale = t(apply(data0[row0,col0],2,as.numeric))
colnames(data_pca_noscale) = row0
    
# Comparison of the loadings for dimensions 1 and 2
n.comp = 2
type_var = factor(data0[-(nrow(data0)-c(0:1)),ncol(data0)])
mx.load_scale.t2 = apply(princomp(covmat=cov(data_pca_scale))$loadings[,1:n.comp],2,
                         function(x,id){tapply(x,id,mean)},
                         id = type_var)
mx.load_noscale.t2 = apply(princomp(covmat=cov(data_pca_noscale))$loadings[,1:n.comp],2,
                           function(x,id){tapply(x,id,mean)},
                           id = type_var)
par(mfrow=c(1,2))
ylimw = max(abs(c(mx.load_scale.t2,mx.load_noscale.t2)))*c(-1,1)
xlimw = c(1,nlevels(type_var))+c(-.5,.5)
for(dw in 1:2){
    plot(1,1,ylim=ylimw,xlim=xlimw,axes=FALSE,xlab="",ylab="loadings",
         main=paste0("PCA dim ",dw))
    abline(h=0,col="light gray",lwd=2)
    points(mx.load_scale.t2[,dw],col="red")
    lines(mx.load_scale.t2[,dw],col="red")
    points(mx.load_noscale.t2[,dw],col="blue")
    lines(mx.load_noscale.t2[,dw],col="blue")
    axis(2,las=2,cex.axis=.7)
    axis(1,1:nlevels(type_var),levels(type_var),las=2,cex.axis=.7)
    legend("top",lty=c(1,1),col=c("red","blue"),legend=c("with scaling","without scaling"),ncol=2)
    }    


##
## PCA on chosen dataset (% var explained by method compared with AJS: OK)
##

data_pca = data_pca_noscale

# boxplots (compared with AJS: OK)
boxplot(data_pca,col="light gray")
boxplot(t(data_pca),col="light gray")

# princomp on covmatrix:
pca1 = princomp(covmat=cov(data_pca))
summary(pca1)
# prcomp on data (with n<p)
pca2 = prcomp(data_pca)
summary(pca2)
names(pca2)
# same results (with sometimes different sign)
for(pw in 1:nrow(data_pca)){cat("dim",pw,":\t",cor(pca2$rotation[,pw],pca1$loadings[,pw]),"\n")}




##
## ids
##

#
n.obs  = nrow(data_pca)
id.obs = data.frame(pos=1:n.obs,id=rownames(data_pca),
             patient = unlist(data0["Case_mRNA",-ncol(data0)]),
             stringsAsFactors=FALSE)

#
temp       = table(id.obs$patient)    
n.patient  = length(temp)
id.patient = data.frame(pos=1:n.patient,id=names(temp),n=c(temp),
                        stringsAsFactors=FALSE)

#
n.var  = ncol(data_pca)
id.var = data.frame(pos=1:n.var,id=colnames(data_pca),
                    type=data0[-(nrow(data0)-c(0:1)),ncol(data0)],
                    stringsAsFactors=FALSE)
    # all(rownames(data0)[-(nrow(data0)-c(0:1))]==id.var$id)

#
temp    = table(id.var$type)    
n.type  = length(temp)
id.type = data.frame(pos=1:n.type,id=names(temp),n=c(temp),
                     stringsAsFactors=FALSE)





##
## save
##
n = nrow(data_pca)
p = ncol(data_pca)
save(data_pca,n.obs,id.obs,n.patient,id.patient,n.var,id.var,n.type,id.type,file=paste0(path.rd,"0-import"))



