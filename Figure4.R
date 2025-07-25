# Figure 4: After running the program in 'Simulations.R', we saved the simulation and estimation output in the image file
#            named "Plot-testing-Feb12-2014.RData" (also included as a Supplementary file).
#            Here we load the simulations and draw Figure 3 in the manuscript

 
#########  STEP 1:  The 4 following functions were used to draw the first column in Figure 3 
Belows.traj<- function(lam,K,beta,no,len){
	
	n.vec <- rep(0,len);
	n.vec[1] <- no
	
	for(i in 2:len){
		
		n.vec[i] <- n.vec[(i-1)]*(lam/(1+ (lam-1)*(n.vec[(i-1)]/K)^(beta) ) );
	}
	
	return(n.vec)
	
}


Belows.map<- function(lam,K,beta,nt.vec){
	
	len<- length(nt.vec);
	ntp1.vec <- rep(0,len);
	
	for(i in 1:len){
		
		ntp1.vec[i] <- nt.vec[i]*(lam/(1+ (lam-1)*(nt.vec[i]/K)^(beta) ) );
	}
	
	return(ntp1.vec)
	
}

Belows.gr<- function(lam,K,beta,nt.vec){
	
	len<- length(nt.vec);
	gr.vec <- rep(0,len);
	
	for(i in 1:len){
		
		gr.vec[i] <- lam/(1+ (lam-1)*(nt.vec[i]/K)^(beta));
	}
	
	return(gr.vec)
	
}


####- Below's Cobweb:

Belows.cobweb<-function(lam,K,beta,start.val,len){

	# len has to be an odd number	
	x.web <- rep(0,len);
	y.web <- x.web;
	x.web[1]<- start.val;
	y.web[1]<- 0;
	
	for(i in seq(2,(len-1),by=2)){
		
		next.pred <- x.web[(i-1)]*(lam/(1+ (lam-1)*(x.web[(i-1)]/K)^(beta) ) );
		x.web[i]  <- x.web[(i-1)];
		y.web[i]  <- next.pred;
		x.web[(i+1)] <- next.pred;
		y.web[(i+1)] <- next.pred;
		}		
	
	return(cbind(x.web,y.web));
}


###### STEP 2:  loading the R image file and drawing Figure 3 in the manuscript

load("Plot-testing-Feb12-2014.RData")
mlesmat.uc <- read.table("mlesmatuc1000b.txt", header=FALSE)
mlesmat.c <- read.table("mlesmatc1000b.txt", header=FALSE)
mlesmat.oc <- read.table("mlesmatoc1000b.txt", header=FALSE)


postscript("Fig4.eps", family="Helvetica", paper="letter")
### Undercompensatory density dependence
K <- 100;
lam <- 4.5;
beta<- 0.5;
Belowsmap <-Belows.map(lam=lam,K=K,beta=beta,nt.vec=0:500)
cobweb.plot<-Belows.cobweb(lam=lam,K=K,beta=beta,start.val=15,len=21)

#par(mfrow=c(3,1),font=1,mai=c(0.55,0.65,0.60,0.40),oma=c(1.45,1.25,0.5,0.25))
par(mfrow=c(3,3),mai=c(0.55,0.85,0.45,0.70),oma=c(1.15,1.15,0.5,0.25))
plot(0:500,Belowsmap,type="l",xlim=c(0,500),col="red", main=expression("Undercompensation: "*beta==0.5),ylab="N(t+1) for Below's map",xlab="N(t)",cex.lab=1.5,lwd=2, cex.main=1.5)
lines(cobweb.plot[,1],cobweb.plot[,2],col=4,type="l");abline(a=0,b=1,lty=2)
#abline(h=K, lty=1);


matplot(1:len2, t(undercomp[1:5,]), type="l", col=rep("blue",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3),xlab="Time", ylab="Log-abundances",cex.lab=1.5, main=expression("Case 1: Under-compensation"), cex.main=1.5)
abline(v=mean(mlesmat.uc[,10]), col="darkgrey", lty=2)
abline(v=thalf.uc, lty=3);
boxplot((mlesmat.uc[,10])/thalf.uc,ylab="Estimated/True",cex.lab=1.5, main=expression("Half-life, case 1"),cex.main=1.5, ylim=c(0,2))
abline(h=1)


#### Compensatory density dependence
K <- 100;
lam <- 4.5;
beta<- 1;
Belowsmap <-Belows.map(lam=lam,K=K,beta=beta,nt.vec=0:500)
cobweb.plot<-Belows.cobweb(lam=lam,K=K,beta=beta,start.val=15,len=21)
plot(0:500,Belowsmap,type="l",xlim=c(0,500),col="red", main=expression("Compensatory dynamics: "*beta==1),ylab="N(t+1) for Below's map",xlab="N(t)",cex.lab=1.5,cex.main=1.5, ylim=c(0,120),lwd=2)
lines(cobweb.plot[,1],cobweb.plot[,2],col=4,type="l");abline(a=0,b=1,lty=2)
#abline(h=K, lty=1);

matplot(1:75, t(comp[1:5,1:75]), type="l", col=rep("blue",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3) ,xlab="Time", ylab="Log-abundances",cex.lab=1.5, main=expression("Case2: Compensatory dynamics"), cex.main=1.5)
abline(v=mean(mlesmat.c[,10]), col="darkgrey", lty=2)
abline(v=thalf.c, lty=3);
boxplot((mlesmat.c[,10])/thalf.c,ylab="Estimated/True",cex.lab=1.5, main=expression("Half-life, case 2"),cex.main=1.5, ylim=c(0,2))
abline(h=1)


#### Overcompensatory density dependence
K <- 100;
lam <- 4.5;
beta<- 1.2;
Belowsmap <-Belows.map(lam=lam,K=K,beta=beta,nt.vec=0:500)
cobweb.plot<-Belows.cobweb(lam=lam,K=K,beta=beta,start.val=15,len=21)
plot(0:500,Belowsmap,type="l",xlim=c(0,500),col="red", main=expression("Overcompensation: "*beta==1.2),ylab="N(t+1) for Below's map",xlab="N(t)",cex.lab=1.5,lwd=2, cex.main=1.5, ylim=c(0,120))
lines(cobweb.plot[,1],cobweb.plot[,2],col=4,type="l");abline(a=0,b=1,lty=2)

matplot(1:75, t(overcomp[6:10,1:75]), type="l", col=rep("blue",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3),xlab="Time", ylab="Log-abundances",cex.lab=1.5, main=expression("Case3: Over-compensation"), cex.main=1.5)
abline(v=mean(mlesmat.oc[,10]), col="darkgrey", lty=2)
abline(v=thalf.c, lty=1);
boxplot((mlesmat.oc[,10])/thalf.oc,ylab="Estimated/True",cex.lab=1.5, main=expression("Half-life, case 3"),cex.main=1.5, ylim=c(0,2))
abline(h=1)

dev.off()
