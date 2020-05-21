AttrPsyFunc=function(attrib=NULL,boundfunc="biweight",eff=NULL,bdp=NULL,arp=NULL,p=NULL,d=NULL,c=NULL,M=NULL){
	#
	# Tuckey's biweight 
	#
	if(boundfunc=="biweight"){
		if(attrib=="weights"){ifelse(abs(d)<=c,(1-(d/c)^2)^2,0)
		}else{if(attrib=="c"){
			if(!is.null(eff)){
			optimize(f=function(c,eff,p){drop((eff-AttrPsyFunc("eff","biweight",p=p,c=c))^2)},
				 eff=eff,p=p,interval=c(0.01,1000),tol=.Machine$double.eps)[[1]]
			}else{if(!is.null(bdp)){
			optimize(f=function(c,bdp,p){drop((bdp-AttrPsyFunc("bdp","biweight",p=p,c=c))^2)},
				 bdp=bdp,p=p,interval=c(0.01,1000),tol=.Machine$double.eps)[[1]]		
			}else{warning("Specify 'eff' or 'bdp'",call.=F)}}
		}else{if(attrib=="M"){0		
		}else{if(attrib=="p1"|attrib=="a"){
			(((945-420*c^2+90*c^4-12*c^6+c^8)+2*c*(105*c^2-13*c^4+c^6-945)*dnorm(c)/(2*pnorm(c)-1))/c^8)*pchisq(c^2,1)
		}else{if(attrib=="p2"|attrib=="b"){
			(((15-6*c^2+c^4)+2*c*(c^2-15)*dnorm(c)/(2*pnorm(c)-1))/c^4)*pchisq(c^2,1)
		}else{if(attrib=="p3"){
			(((675675-270270*c^2+51975*c^4-6300*c^6+525*c^8-30*c^10+c^12)+
			2*c*(45045*c^2-6930*c^4+558*c^6-31*c^8+c^10-675675)*dnorm(c)/(2*pnorm(c)-1))/c^12)*pchisq(c^2,1)
		}else{if(attrib=="p4"){
			(((525-240*c^2+54*c^4-8*c^6+c^8)+2*c*(65*c^2-9*c^4+c^6-525)*dnorm(c)/(2*pnorm(c)-1))/c^8)*pchisq(c^2,1)
		}else{if(attrib=="p5"){
			(((105-60*c^2+18*c^4-4*c^6+c^8)+2*c*(25*c^2-5*c^4+c^6-105)*dnorm(c)/(2*pnorm(c)-1))/c^8)*pchisq(c^2,1)
		}else{if(attrib=="p6"){
			(((2625-900*c^2+138*c^4-12*c^6+c^8)+2*c*(25*c^2-13*c^4+c^6-2625)*dnorm(c)/(2*pnorm(c)-1))/c^8)*pchisq(c^2,1)
		}else{if(attrib=="rho"){(d^2/2-d^4/(2*c^2)+d^6/(6*c^4))*(abs(d)<c)+c^2/6*(abs(d)>=c)
		}else{if(attrib=="psy"){ifelse(abs(d)<= c,d*(1-(d/c)^2)^2,0)
		}else{if(attrib=="e1"){
			(dchisq(c^2,p)*2*c^4*(c^6-3*c^4*p-p*(40+p*(10+p))+c^2*(32+p*(10+3*p))))/(p*(2+p)*(4+p)*(6+p)*(8+p))+
			(c^8-4*c^6*(2+p)+6*c^4*(2+p)*(4+p)-4*c^2*(2+p)*(4+p)*(6+p)+(2+p)*(4+p)*(6+p)*(8+p))*
			 pgamma(c^2/2,shape=5+p/2)/c^8
		}else{if(attrib=="e2"){
			dchisq(c^2,p)*(c^2*(2+c^2-p)-(8+(c^2-p)^2+6*p))/(p*(p/2+1))+
			(c^4-2*c^2*(2+p)+(2+p)*(4+p))*(pgamma(c^2/2,shape=(2+p)/2))/c^4
		}else{if(attrib=="eff"){
			AttrPsyFunc("e2","biweight",p=p,c=c)^2/AttrPsyFunc("e1","biweight",p=p,c=c)
		}else{if(attrib=="bdp"){
			(1/(c^6))*(p*(3*c^4-3*c^2*(2+p)+(2+p)*(4+p))-
			 3*c^4*p*(1-pgamma(c^2/2,shape=1+p/2))+3*c^2*p*(2+p)*(1-pgamma(c^2/2,shape=2+p/2))- 
			(p*(2+p)*(4+p))*(1-pgamma(c^2/2,shape=3+p/2))+c^6*(1-pgamma(c^2/2,shape=p/2)))
		}else{if(attrib=="arp"){1-pchisq(c^2,p)
		}}}}}}}}}}}}}}}}
	}else{stop("Specify a valid/possible bounded function [boundfunc]",call=T)}
	}

		
	
	
	

