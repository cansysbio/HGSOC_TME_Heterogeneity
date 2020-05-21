covREM = function(data,est="s",start="mcd",boundfunc="biweight",bdp=0.4,arp=NULL,maxiter=500,eps=10^(-5),
                  quietly=F){
	##############################################################################	
	# useful functions
	find_start=function(full=y,complete=yd,Par=Par,mp.i=mp.i){
		#full=y;complete=yd;Par=Par;mp.i=mp.i
		small=nrow(complete)<=(ncol(complete)+1)
		# classic
		if(!Par$robust|Par$start=="mle"){
			Par$start="mle"
			if(small){out=list(center=apply(full,2,mean,na.rm=T),cov=cov(full,use="pairwise.complete.ob"),Par=Par)
			}else{out=list(center=apply(complete,2,mean,na.rm=T),cov=cov(complete),Par=Par)}
			# for compatibility reason
			if(!Par$na){
				out$cor=sweep(sweep(out$cov,1,sqrt(diag(out$cov)),"/"),2,sqrt(diag(out$cov)),"/")
				out$data=complete
				out$dist=sqrt(mahalanobis(complete,out$center,out$cov))
				out$stand_dist=((out$dist^2/Dim$p)^(1/3)-1+2/(9*Dim$p))/sqrt(2/(9*Dim$p))
				out$dim=Dim
				out$w=rep(1,Dim$n)
				class(out)="covREM"
				}
		out  
		# robust
		}else{
		# if small
		if(small&Par$start!="em-rogk"){
			if(!Par$quietly){warning(paste("\n'n < (p+1)' in the complete dataset",
				     "\nThe ROGK with imputation is used as starting point.\n",sep=""),call.=F,immediate.=T)}
				     Par$start="em-rogk";find_start(full,complete,Par,mp.i)
		}else{
		# mcd
		if(Par$start=="mcd"){
			require(robustbase)
			out=try(covMcd(complete),silent=T)
			if(class(out)=="try-error"){
				if(!Par$quietly){warning(paste("\nFAST MCD impossible on the complete dataset.",
				"\nThe ROGK is used as starting point.\n",sep=""),call.=F,immediate.=T)}
				Par$start="rogk";find_start(full,complete,Par)
			}else{list(center=out$center,cov=out$cov,Par=Par)}
		# rogk
		}else{if(Par$start=="rogk"){RobEM(complete,1,1,Par,rep(1,n.id),est="rogk")
		# s
		}else{if(Par$start=="s"){
			Par$start="mcd"
			start=find_start(full,complete,Par)
			RobEM(complete,start$center,start$cov,Par,rep(1,n.id),est="s")
		# em-rogk
		}else{if(Par$start=="em-rogk"){
			start=RobEM(complete,1,1,Par,rep(1,n.id),est="rogk")
			RobEM(full,start$center,start$cov,Par,mp.i,est="rogk")
		}}}}}}
		}
	robpar=function(est,boundfunc="t-biweight",eff=NULL,bdp=NULL,arp=NULL,
		m,n,eps,maxiter,mc,r,quietly,start){
		###
		Par=list(mc=mc,est=est,maxiter=maxiter,eps=eps,na=sum(r)>0,robust=est!="mle",start=start,
		         quietly=quietly)
		if(est=="mle"){Par=c(Par,c=1e+6,M=0,eff=1,bdp=0,arp=0)
		}else{  # possibly adapt bdp to m and n
			if(is.null(bdp)){bdp=(n-m)/(2*n)
			}else{if(bdp>(n-m)/(2*n)){bdp=(n-m)/(2*n);if(!quietly){cat("-> S-Step: bdp set to",round(bdp,3),"[chosen bdp was to high given n and m]\n")}}}
			### s-estimator 
			Par$s=list()
			Par$s$boundfunc=ifelse(is.null(names(boundfunc)),boundfunc[1],boundfunc[["s"]])
			Par$s$c=AttrPsyFunc("c",Par$s$boundfunc,p=m,bdp=bdp,arp=arp,eff=NULL)
			Par$s$M=AttrPsyFunc("M",Par$s$boundfunc,p=m,bdp=bdp,c=Par$s$c,arp=arp,eff=eff)
			# if solution has bdp>bpd*
			if((AttrPsyFunc("bdp",Par$s$boundfunc,p=m,c=Par$s$c,M=Par$s$M)-bdp)>0.00001){
				if(!quietly){cat("-> S-Step: bdp set to",round(bdp,3),"[chosen arp/eff induced a too high bdp given n and m]\n")}
				Par$s$c=AttrPsyFunc("c",Par$s$boundfunc,p=m,bdp=bdp,arp=NULL,eff=NULL)
				Par$s$M=AttrPsyFunc("M",Par$s$boundfunc,p=m,c=Par$s$c,bdp=bdp,arp=NULL,eff=NULL)
				}
			Par$s$arp=AttrPsyFunc("arp",Par$s$boundfunc,p=m,c=Par$s$c,M=Par$s$M)
			Par$s$bdp=AttrPsyFunc("bdp",Par$s$boundfunc,p=m,c=Par$s$c,M=Par$s$M)
			Par$s$eff=AttrPsyFunc("eff",Par$s$boundfunc,p=m,c=Par$s$c,M=Par$s$M)
			}
		Par
		}
	RobEM=function(y,mu,sigma,Par,mp.i,est){
		# Robust EM [uses sweep.mp() for imputation and rogk() or s_est() for Estimation]
		# useful function
		sweep.mp=function(y_mp,G){
			# Sweep algorithm [refer to Little & Rubin (1987, p. 112)] 
			# SC/DLC june2010
			n.i=dim(y_mp)[1]	
			n.p=dim(y_mp)[2]
			idp.miss=seq(1:n.p)[is.na(y_mp[1,])]	
			idp.nmiss=seq(1:n.p)[!is.na(y_mp[1,])]
			np.nmiss=length(idp.nmiss)
			y_imp=y_mp	
			## if missing values
			if(n.p-np.nmiss>0){
			for(pw in 1:np.nmiss){
				idpw=idp.nmiss[pw]+1
				G1=G-kronecker(G[idpw,,drop=F],matrix(rep(1,n.p+1),n.p+1))*
				   kronecker(matrix(rep(1,n.p+1),1),G[,idpw,drop=F])/G[idpw,idpw]
				G1[idpw,]=G1[,idpw]=G[idpw,]/G[idpw,idpw]
				G1[idpw,idpw]=-1/G[idpw,idpw]
				G=G1
				}
			# imputation
			y_imp[is.na(y_imp)]=
				(apply(array(t(kronecker(matrix(1,n.i),G[-1,-1])*
				kronecker(as.matrix(y_mp),matrix(1,n.p))),
				dim=c(n.p,n.p,n.i)),2,function(x)apply(x,2,sum,na.rm=T))+
				kronecker(G[1,-1,drop=F],matrix(1,n.i)))[is.na(y_imp)]
			# correction
			cov_cor=matrix(0,n.p,n.p)
			cov_cor[idp.miss,idp.miss]=G[-1,-1][idp.miss,idp.miss]
			## else
			}else{cov_cor=matrix(0,n.p,n.p)}
			## output
			list(y_imp=t(y_imp),cov_cor=cov_cor)
			}
		positive_eigen=function(x,limit=10^(-5),value=10^(-5)){
			LDtL=eigen(x)
			LDtL$values[LDtL$values<=limit]=value
			LDtL$vectors%*%diag(LDtL$values)%*%t(LDtL$vectors)
			}
		s_est=function(y_imp,mu,sigma,Par){	
			# weights 
			p=ncol(y_imp)
			n=nrow(y_imp)
			d=try(sqrt(mahalanobis(y_imp,mu,sigma)),silent=T)
			if(class(d)=="try-error"){
				sigma=positive_eigen(sigma)
				d=sqrt(mahalanobis(y_imp,mu,sigma))
				}
			# correction in order to respect the constraint of S-Estimator
			h<-floor((n+p+1)/2)
			d<-(d*sqrt(qchisq(h/(n+1),p)))/(sort(d)[h])
			# update estimates
			wbs=AttrPsyFunc("weights",boundfunc=Par$s$boundfunc,d=d,c=Par$s$c,M=Par$s$M)
 			vbs=d*AttrPsyFunc("psy",boundfunc=Par$s$boundfunc,d=d,c=Par$s$c,M=Par$s$M)
			mu=as.numeric((wbs%*%y_imp)/sum(wbs))
			sigma=p*matrix(apply(t(apply(sweep(y_imp,2,mu),1,function(x)x%*%t(x)))*wbs,2,sum),ncol=p)/sum(vbs)
			# out
			list(center=mu,cov=sigma,w=wbs,Par=Par)
			}
		rogk=function(data, c1 = 4.5, c2 = 3){
			# Robust OGK (similar to covOGK of 'robustbase' when sigmamu=s_mad)
			# useful function
			gnan0 = function(x, c1,c2){
				d <- apply(x, 2, mad)
				w <- sweep(x, 2, apply(x, 2, median), FUN = "-")
				w.c <- sweep(w, 2, d, FUN = "/")
				I <- ifelse(abs(w.c) <= c1, 1, 0)
				wc <- (1 - (w.c/c1)^2)^2 * I
				vec.mu <- apply(x * wc, 2, sum)/apply(wc, 2, sum)
				mat <- sweep(x, 2, vec.mu, FUN = "-")
				mat.c <- sweep(mat, 2, d, FUN = "/")
				rhoc <- ifelse(mat.c^2 <= c2^2, mat.c^2, c2^2)
				sig2 <- apply(rhoc, 2, mean) * d^2
				return(list(vec.mu=vec.mu,sig2=sig2))
				}
			# sizes
			n <- nrow(data)
			p <- ncol(data)
			k <- ((p - 1) * p)/2
			# compute D and standardize observations y
			D <- diag(sqrt(gnan0(data,c1,c2)$sig2))
			y <- sweep(data, 2, diag(D), FUN = "/")
			# compute matrix U, eigenvectors and eigenvalues
			pp <- matrix(1:(p^2),p,p)
			yp <- ym <- pos <- NULL
			for(i in 1:(p-1)) {
				for(j in (i+1):p){
					yp <- cbind(yp, y[, i] + y[, j], deparse.level = 1)
					ym <- cbind(ym, y[, i] - y[, j], deparse.level = 1)
					pos <- c(pos,pp[j,i])
				}
			}
			U <- diag(rep(1/2,p))
			U[pos] = 0.25*(gnan0(yp,c1,c2)$sig2-gnan0(ym,c1,c2)$sig2)
			U <- U + t(U)
			E <- eigen(U)$vectors
			# compute mu and cov
			A <- D%*%E
			final <- gnan0(y%*%E,c1,c2)
			mu <- as.numeric(A%*%final$vec.mu)
			V <- A%*%diag(final$sig2)%*%t(A)
			names(mu)=colnames(V)=rownames(V)=colnames(data)
			# compute weights
			mal <- mahalanobis(data, mu, V)
			do <- (qchisq(0.9, p) * median(mal))/qchisq(0.5, p)
			w <- ifelse(mal <= do, 1, 0)
			# output
			list(center=mu,cov=V,w=w)
			}
		# get info
		n.i=dim(y)[1]	
		n.p=dim(y)[2]
		maxiter=if(est=="rogk"){50}else{Par$maxiter}
		# optimisation
		iter=1
		crit=1000
		while((iter <= maxiter) & (crit>Par$eps)){
			mu.old=mu
			sigma.old=sigma
			#### imputation
			y_mpw=lapply(split(y,mp.i),function(x)matrix(x,ncol=n.p))
			sweep_mpw=lapply(y_mpw,sweep.mp,G=rbind(c(-1,mu),cbind(mu,sigma)))
			y_imp=matrix(unlist(lapply(sweep_mpw,function(x)x[[1]])),nrow=n.i,byrow=T)
			#### estimate mu and sigma
			if(est=="mle"){fit=list(center=apply(y_imp,2,mean),cov=cov(y_imp),w=rep(1,n.i))
			}else{if(est=="rogk"){fit=rogk(y_imp)
			}else{fit=s_est(y_imp,mu,sigma,Par);Par=fit$Par}}
			mu=fit$center
			sigma=fit$cov
			w=fit$w
			# convergence criterion
			crit=max(abs(c(mu-mu.old,sigma-sigma.old)))
			iter=iter+1
			}
		# output
		names(mu)=rownames(sigma)=colnames(sigma)=colnames(y_imp)=colnames(y)
		dist=try(sqrt(mahalanobis(y_imp,mu,sigma)),silent=T)
		if(class(dist)=="try-error"){
			sigma=positive_eigen(sigma)
			dist=sqrt(mahalanobis(y_imp,mu,sigma))
			}
		# Wilson-Hilferty transformation
		p.i=apply(!is.na(y),1,sum)
		stand_dist=((dist^2/(p.i))^(1/3)-1+2/(9*p.i))/sqrt(2/(9*p.i))
		# output
		Par=c(Par,iter=iter-1)
		out=list(center=mu,cov=sigma,cor=sweep(sweep(sigma,1,sqrt(diag(sigma)),"/"),2,sqrt(diag(sigma)),"/"),
			data=y,data_imp=if(iter==2){NULL}else{y_imp},w=w,dist=dist,stand_dist=stand_dist,
			Par=Par,dim=list(n=n.i,p=n.p),misspattern=NULL)
		class(out)="covREM"
		out
		}
	##############################################################################
	# prepare									
	mc=match.call()
	if(is.na(match(est,c("s","rogk","mle")))){
		stop("Non available 'est' specified",call.=F)}
	if(is.na(match(start,c("s","rogk","mle","mcd","em-rogk")))){
		stop("Non available 'start' specified",call.=F)}
	if(is.na(match(boundfunc,c("biweight","t-biweight","lmws")))){
		stop("Non available 'bounded function' specified",call.=F)}		
	if(any(apply(data,2,class)=="factor"|apply(data,2,class)=="character")){
		stop("Only numerical data are allowed, check column classes.",call.=F)}
	data <- as.matrix(data)
	data <- data[apply(is.na(data),1,sum)!=ncol(data),]
	data <- data[apply(is.na(data),2,sum)!=nrow(data),]
	Dim  <- list(n=nrow(data),p=ncol(data))
	if(Dim$n <= Dim$p){stop("n <= p -- you can't be serious!",call.=F)}
	if(Dim$n <= 1|Dim$p <= 1){stop("covREM fits multivariate datasets, check the dataset size",call.=F)}
	r <- is.na(data)
	##############################################################################
	# Parameters									
	Par=robpar(est=est,boundfunc=boundfunc,eff=eff,bdp=bdp,arp=arp,
		   m=Dim$p,n=Dim$n,eps=eps,maxiter=maxiter,mc=mc,r=r,quietly=quietly,start=start)
	##################################################################################
	# Pattern of missing values (misspattern) in the dataset				
	if(Par$na){
		mdp <- as.integer((r%*%(2^((1:Dim$p)-1)))+1)
		ro <- order(mdp)
		y  <- data[ro,]
		m  <- sort(unique(mdp))
		ni.mp <- tabulate(mdp)[m]
		n.mp  <- length(m)
		mp.i  <- rep(1:n.mp,ni.mp)
		id.mp <- 1-r[ro,][cumsum(ni.mp),,drop=F]
		np.mp <- as.vector(apply(id.mp, 1, sum))
		yd <- na.omit(y)
		n.id <- dim(yd)[1]	
	}else{  Par$start <- Par$est; y <- yd <- data; n.id <- Dim$n; mp.i=NA} 
	##################################################################################
	# Starting point (= final point if no missing)					
	if(is.character(start)){
		mu0sigma0=find_start(y,yd,Par,mp.i)	
	}else{mu0sigma0=start}
	##################################################################################
	# Estimation					
	## if no missing values -> output 
	if(!Par$na){
		if(!Par$quietly){warning("No missing values -> EM not performed",call.=F)}
		mu0sigma0
	## if missing -> EM
	}else{  mu.0=mu0sigma0$center
		sigma.0=mu0sigma0$cov
		Par=mu0sigma0$Par
		# EM										
		fit=RobEM(y,mu.0,sigma.0,Par,mp.i,est=Par$est)	
		# reorder the dataset
		fit$data=data
		fit$data_imp=fit$data_imp[order(ro),]
		fit$dist=fit$dist[order(ro)]
		fit$stand_dist=fit$stand_dist[order(ro)]		
		# add information about missing pattern and starting values
		fit$misspattern=list(id=id.mp,n=ni.mp,p=np.mp)
		# this is it 
		fit
		}
	}
	
	
print.covREM=function(obj,cor=F,digits=3,full=F){
	# get info
	Par=obj$Par
	poss.est=c("MLE","Robust OGK","Robust OGK with EM-Algorithm","Biweight S-Estimator","FAST MCD")
	poss.start=c("complete","complete","all","complete","complete")
	Par.est=c("mle","rogk","rogk","s","")
	names(poss.est)=names(Par.est)=names(poss.start)=c("mle","rogk","em-rogk","s","mcd")
	# call, title, est and convergence
	cat(paste("\nCall:",deparse(Par$mc),"\n"))
	if(full){cat("\n---\n")}
	if(full){
		cat(paste("\n",if(Par$robust){"Robust "},"Multivariate Location And Scale Estimates\n",sep=""))
		cat(paste("\nFit by ",poss.est[Par$est],if(Par$na){" with EM-Algorithm"},
		if(!is.null(Par$bdp)){paste(" (",round(Par$bdp*100,digits),"% BDP)",sep="")},sep=""))
		if(Par$na|Par$est=="s"){
			cat(paste("\n\tStarting values = ",poss.est[Par$start]," on ",poss.start[Par$start]," observations",sep=""))
			cat(paste("\n\tConvergence =",Par[[Par.est[Par$est]]]$crit<Par[[Par.est[Par$est]]]$eps,"\n"))
		}else{cat("\n")}}
	# center
	cat("\nLocation estimates:\n");print(round(obj$center,digits))
	# covariance or correlation
	if(cor){cat("\nCorrelation matrix:\n");print(format(round(as.data.frame(obj$cor),digits)))
	}else{cat("\nCovariance matrix:\n");print(format(round(as.data.frame(obj$cov),digits)))}
	#  missing pattern
	if(full&Par$na){
		p=ncol(obj$misspattern$id)
		n.mp=nrow(obj$misspattern$id)
		matx=as.data.frame(matrix(NA,ncol=p+6,nrow=n.mp))
		matx[,c(1,p+c(2,4,6))]="|"
		colnames(matx)=c("|",colnames(obj$misspattern$id),"|","obs","|","p","|")
		matx$obs=obj$misspattern$n
		matx$p=obj$misspattern$p
		matx[,colnames(obj$misspattern$id)]=apply(obj$misspattern$id,2,function(x){x[x==1]="x";x[x==0]="";x})
		cat("\nMissing Pattern:\n")
		print(matx)
		}
	if(full){cat("\n---\n\n")}
	}

summary.covREM=function(obj,cor=F,digits=3){
	print.covREM(obj,cor=cor,digits=digits,full=T)
	}

plot.covREM = function(obj1,obj2=NULL,id.n=3,compare=T,ask=T,which=1:3,label.obj2=NULL,col.points=c("#0080ff","#ff0000","green"),
                       col.ellipse=c("dark blue","brown4"),plot.mu_sd=T,round.mu_sd=2,cex.mu_sd=.75,tol.ellipse=.95,cex.names=1,
                       plot.cor_ellipse=T,cex.cor=1,round.cor=2){
	## useful functions
	panel.dist <- function(x,y,qw=qw,id.n=id.n,h=h,v=v,col) {
		panel.xyplot(x, y, pch = 21,col=col)
		n <- length(y)
		if (id.n > 0) {
		    yp.idx <- order(y)[(n - id.n + 1):n]
		    if (any(yp.idx)) 
			panel.text(x[yp.idx], y[yp.idx],paste(" ", yp.idx,sep = ""), adj = 0,cex=.8)
		}
		if(h){panel.abline(h = c(-qw,qw), lty = 2,col="gray")}
		if(v){panel.abline(v = c(-qw,qw), lty = 2,col="gray")
		      panel.abline(0,1,lty=1,col="gray")}
		invisible()
	    }
	ellipse=function(i,j,obj,level){
		points<-c(0:180,0)*pi/90
		mu=list(i=obj$center[i],j=obj$center[j])
		sigma=list(i=sqrt(diag(obj$cov)[i]),j=sqrt(diag(obj$cov)[j]))
		corr=obj$cor[i,j]
		chi=sqrt(qchisq(level,2))
		coord_cor=list(i=chi*cos(points+acos(corr)/2),
			       j=chi*cos(points-acos(corr)/2))
		coord_cov=list(i=coord_cor$i*sigma$i+mu$i,
			       j=coord_cor$j*sigma$j+mu$j)
		list(mu=mu,sigma=sigma,cor=corr,coord_cor=coord_cor,coord_cov=coord_cov)	       
		}
	## checks
	require(lattice)
	if(compare){which.plot=rep(T,3)}else{which.plot=c(T,F,T)}
	if(!is.null(which)){which.plot[-which]=F}
	if(sum(which.plot)==0){stop("incompatible arguments (which, compare)",call.=F)}        	
	if(sum(which.plot)==1){ask=F}
	if(ask){oask <- devAskNewPage(TRUE)
		on.exit(devAskNewPage(oask))
		}
	# obj2
	if(!is.null(obj2)){
		compare = T
		if(!is.list(obj2)){stop("obj2 has to be a list with elements 'center' and 'cov'",call.=F)}
		if((length(obj2$center)!=obj1$dim$p)){stop("center of obj1 and 2 are not compatible",call.=F)}
		if(all(dim(obj2$cov)!=rep(obj1$dim$p,2))){stop("cov of obj1 and 2 are not compatible",call.=F)}
		if(is.null(obj2$cor)){obj2$cor=sweep(sweep(obj2$cov,1,sqrt(diag(obj2$cov)),"/"),2,
			sqrt(diag(obj2$cov)),"/")}
		if(length(obj2$dist)!=obj1$dim$n){obj2$dist=sqrt(mahalanobis(
			if(obj1$Par$na){obj1$data_imp}else{obj1$data},obj2$center,obj2$cov))}
		if(length(obj2$stand_dist)!=obj1$dim$n){
			p.i=apply(!is.na(obj1$data),1,sum)
			obj2$stand_dist=((obj2$dist^2/(p.i))^(1/3)-1+2/(9*p.i))/sqrt(2/(9*p.i))}
		if(is.null(obj2$Par)){obj2$Par=list(est="self")}
	}else{if(compare){
		if(obj1$Par$robust){obj2=covREM(obj1$data,est="mle",quietly=T)
		}else{obj2=covREM(obj1$data,est="s",quietly=T)}
		}
	}
	# name of the est(s)
	poss.est=c("MLE","ROGK","S-Est.","Self Defined")
	names(poss.est)=c("mle","rogk","s","self")
	obj1$name=paste(if(compare){"[1]"},poss.est[obj1$Par$est],if(obj1$Par$na){"with EM"},
		if(obj1$Par$est=="s"){paste("(",round(obj1$Par$bdp*100,2),"%)",sep="")})	
	if(compare){obj2$name=paste("[2]",
		if(obj2$Par$est=="self"&!is.null(label.obj2)){label.obj2
		}else{paste(poss.est[obj2$Par$est],if(obj1$Par$na&obj2$Par$est!="self"){"with EM"})},
		if(obj2$Par$est=="s"){paste("(",round(obj2$Par$bdp*100,2),"%)",sep="")})}
	# different elements
	Limit=qnorm(0.975)
	Col=rep(col.points[1],obj1$dim$n)
	Col[abs(obj1$stand_dist)>Limit]=col.points[2]
	Names=colnames(obj1$data)
	if(!obj1$Par$na){obj1$data_imp=obj1$data}	
	## plot 1
	if(which.plot[1]){
	if(!ask){dev.new()}
		data=data.frame(dist=c(obj1$stand_dist,if(compare){obj2$stand_dist}),index=rep(1:obj1$dim$n,compare+1),
			est=rep(c(obj1$name,if(compare){obj2$name}),each=obj1$dim$n))
		print(xyplot(dist~index|est,data=data,ylab="Wilson-Hilferty Transformed Mahalanobis Distances",xlab="Index",
			panel=panel.dist,id.n=id.n,qw=Limit,h=T,v=F,col=rep(Col,compare+1)))
		}	
	## plot 2
	if(which.plot[2]){
	if(!ask){dev.new()}
		data=data.frame(dist1=obj1$stand_dist,dist2=obj2$stand_dist,
			est="Comparison of Wilson-Hilferty Transformed Mahalanobis Distances")
		print(xyplot(dist1~dist2|est,data=data,ylab=obj1$name,xlab=obj2$name,
		    panel=panel.dist,id.n=id.n,qw=Limit,h=T,v=T,col=Col))
		}
	## plot 3
	if(which.plot[3]){
		if(!ask){dev.new()}
		par(mfrow=c(obj1$dim$p,obj1$dim$p),mar=rep.int(0,4),oma=c(4, 4, 7, 6))
		for(i in 1:obj1$dim$p){
		for(j in 1:obj1$dim$p){
			# if outside diagonal
			if(i!=j){
				# info 
				coord1=ellipse(i,j,obj1,tol.ellipse)
				if(compare){coord2=ellipse(i,j,obj2,tol.ellipse)}else{coord2=NULL}
				# COVARIANCE-ELLIPSES
				if(i<j){
					# color of imputed data
					Col2=Col
					Col2[apply(is.na(obj1$data[,c(i,j)]),1,sum)>0]=col.points[3]
					# plot points
					plot(obj1$data_imp[,j],obj1$data_imp[,i],pch=21,col=Col2,axes=F,cex=.4,
					     ylim=range(coord1$coord_cov$i,coord2$coord_cov$i,obj1$data_imp[,i]),
					     xlim=range(coord1$coord_cov$j,coord2$coord_cov$j,obj1$data_imp[,j]))		    
					# plot ellipses
					lines(coord1$coord_cov$j,coord1$coord_cov$i,col=col.ellipse[1],lwd=2)
					abline(h=coord1$mu$i,col=col.ellipse[1],lwd=.5,lty=2)	
					abline(v=coord1$mu$j,col=col.ellipse[1],lwd=.5,lty=2)	
					if(compare){
						lines(coord2$coord_cov$j,coord2$coord_cov$i,col=col.ellipse[2],lwd=2)
						abline(h=coord2$mu$i,col=col.ellipse[2],lwd=.5,lty=2)	
						abline(v=coord2$mu$j,col=col.ellipse[2],lwd=.5,lty=2)	
						}
					box()
					# axes
					if(i==1){axis(3,cex.axis=.75,padj=if(j%%2==0){.5}else{1.5})}
					if(j==obj1$dim$p){axis(4,cex.axis=.75,padj=if(i%%2==0){-.5}else{-1.5})}
				# CORRELATIONS
				}else{	plot(0,0,pch="",axes=F,ylim=c(-3,3),xlim=c(-3,3))
					box()		
					# ellipse
					if(plot.cor_ellipse){
					lines(coord1$coord_cor$j,coord1$coord_cor$i,col=col.ellipse[1],lwd=.3)
					if(compare){lines(coord2$coord_cor$j,coord2$coord_cor$i,col=col.ellipse[2],lwd=.3)}}		
					# text
					cor.ij=format(round(c(coord1$cor,coord2$cor),round.cor))
					text(0,0,cor.ij[1],pos=if(compare){3}else{NULL},col=col.ellipse[1],cex=cex.cor)
					if(compare){text(0,0,cor.ij[2],pos=1,col=col.ellipse[2],cex=cex.cor)}
					}					
			# if diagonal element
			}else{plot(0,0,pch="",axes=F)
				box()
				text(0,0,Names[i],cex=cex.names)
				if(plot.mu_sd){
				# mu
				pos1=legend("topleft","1",plot=F,cex=cex.mu_sd)
				text(pos1$rect$left,pos1$rect$top-pos1$rect$h/2,pos=4,round(obj1$center[i],round.mu_sd),col=col.ellipse[1],cex=cex.mu_sd)
				if(compare){pos1=legend("topright","1",plot=F,cex=cex.mu_sd)
				text(pos1$rect$left+pos1$rect$w,pos1$rect$top-pos1$rect$h/2,pos=2,round(obj2$center[i],round.mu_sd),col=col.ellipse[2],cex=cex.mu_sd)}
				# sd
				pos2=legend("bottomleft","1",plot=F,cex=.75)
				text(pos2$rect$left,pos2$rect$top-pos2$rect$h/2,pos=4,round(sqrt(obj1$cov[i,i]),round.mu_sd),col=col.ellipse[1],cex=cex.mu_sd)
				if(compare){pos2=legend("bottomright","1",plot=F,cex=cex.mu_sd)
				text(pos2$rect$left+pos2$rect$w,pos2$rect$top-pos2$rect$h/2,pos=2,round(sqrt(obj2$cov[i,i]),round.mu_sd),col=col.ellipse[2],cex=cex.mu_sd)}
				# axis
				if(i==1|i==obj1$dim$p){
				axis(if(i==1){2}else{4},at=c(pos1$rect$top-pos1$rect$h/2,pos2$rect$top-pos2$rect$h/2),
					labels=c(expression(hat(mu)),expression(hat(sigma))),tick=T,las=2,lwd.tick=.5)
					}
				}}
			}}
		# axis 1
		outliers=if(compare){"outliers for [1]"}else{"outliers"}
		imputed_name=if(obj1$Par$na){if(compare){"imputed data by [1]"}else{"imputed data"}}else{""}
		imputed_coma=if(obj1$Par$na){", "}else{""}
		mtext(bquote(phantom("observed data")*.(imputed_coma)*phantom(.(imputed_name))*","*phantom(.(outliers))),
		      side=1,outer=T,line=1)
		mtext(bquote("observed data"*phantom(.(imputed_coma))*phantom(.(imputed_name))*phantom(", ")*phantom(.(outliers))),
		      side=1,outer=T,line=1,col=col.points[1])
		mtext(bquote(phantom("observed data")*phantom(.(imputed_coma))*.(imputed_name)*phantom(", ")*phantom(.(outliers))),
		      side=1,outer=T,line=1,col=col.points[3])
		mtext(bquote(phantom("observed data")*phantom(.(imputed_coma))*phantom(.(imputed_name))*phantom(", ")*.(outliers)),
		      side=1,outer=T,line=1,col=col.points[2])
		# axis 2
		mtext("Correlations",side=2,outer=T,col="black",line=1)
		# axis 3
		mtext(if(compare){"Multivariate Normal Parameters Comparison"}else{"Multivariate Normal Parameters"},side=3,outer=T,col="black",line=4,cex=1.5)
		if(compare){
			mtext(bquote(phantom(.(obj1$name))*" vs "*phantom(.(obj2$name))),side=3,outer=T,col="black",line=2)
			mtext(bquote(.(obj1$name)*phantom(" vs ")*phantom(.(obj2$name))),side=3,outer=T,col=col.ellipse[1],line=2)
			mtext(bquote(phantom(.(obj1$name))*phantom(" vs ")*.(obj2$name)),side=3,outer=T,col=col.ellipse[2],line=2)
		}else{mtext(obj1$name,side=3,outer=T,col=col.ellipse[1],line=1.75)}
		# axis 4
		mtext("95% Tolerance ellipses",side=4,outer=T,col="black",line=3)
	}# end 3
	}# end plot


