### Simulations matching result 1 and result 4.3
###  Knowing how the variance changes implies that one can characterize how the probabilities of pseudeo-extinction change.
###  Not only that, but we can now quantify how the measures of stochastic stability specified by Ives et al change after the 
###  change point, because we know the covariances...

########## The functions used for this plot are saved in BPGSSToolkit.R ##########
source("BPGSSToolkit.R")

########## Drawing the figure with the mean and the variance decay #########


len1 <- 45;
tau1 <- 16;
a1 <- 0.3929;
c1 <- 0.9;
sigmasq1 <- sigmasq2 <- 0.07#0.09725;
a2 <- 0.39;
c2 <- 0.75;

mu1 <- a1/(1-c1);
mu2 <- a2/(1-c2);
t.half <- log(1/2)/log(abs(c2));
mid.thalf <- tau1+t.half/2;
mean.mean <- (mu1+mu2)/2;
nu1 <- sigmasq1/(1-c1*c1);
nu2 <- sigmasq2/(1-c2*c2);


lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1)
mean.t <- lowleft.m$mu.vec;
vars.t <- diag(lowleft.m$Sigma);
Vcov.t <- lowleft.m$Sigma; 

# 5 Simulations:
nsims <- 5
ts.reps <- my.rmvn(n=nsims, mu.vec=mean.t, cov.mat = Vcov.t)
ts.reps4var <- my.rmvn(n=10000, mu.vec=mean.t, cov.mat = Vcov.t) 
emp.covs <- var(t(ts.reps4var))
emp.vars <- diag(emp.covs)
emp.means <- apply(t(ts.reps4var),2,mean)

dim(Vcov.t)
dim(emp.covs)


plot(Vcov.t, emp.covs, pch=16)
abline(a=0,b=1, lwd=2)

# Set
# i = 9, j= 17 = tau1 + 1 = 16 +1
# look for matching quantities in theoretical matrix 
# position 9+1 and 17+1
Vcov.t[18,10]
nu1*(c1^(tau1-9))*(c2^(17-tau1))
emp.covs[18,10]




#pdf("~/Documents/UFL/JM/BPGSSmodel/TPB/R-Programs/AppendixFig.pdf", width=7,height=7, family="Helvetica", paper="special")

postscript("Fig3.eps", family="Helvetica", paper="letter")

par(mfrow=c(2,2),mai=c(0.45,0.75,0.48,0.58),oma=c(1.35,1.25,0.5,0.685) );

# First plot: The Expected value of the process over time
plot(1:len1, mean.t , type="l", col="darkgrey", lwd=2, xlab="", ylab="Log-abundance",ylim=c(mu2*0.2,mu1*1.5), xlim=c(-2,len1+5),cex.lab=1.5, main=expression(E(X(t+tau))==mu[1]*c[2]^t + mu[2]*(1-c[2]^t)),cex.main=1.5,bty="n", xaxt="n", yaxt="n");
arrows(x0=c(0,0), y0=c(mu2*0.35, mu2*0.35), x1=c(0, len1+5), y1=c(mu1*1.5,mu2*0.35),lwd=2, code=2,length=0.1)
text(x=-2,y=mu1, labels=expression(mu[1]),cex=1.25, srt=90);
text(x=-2,y=mu2, labels=expression(mu[2]),cex=1.25, srt=90);
segments(x0=tau1+1, y0=mu2*0.35, x1=tau1+1,y1=mu1*1.25, lwd=2,col="red")
text(x=tau1+1, y= mu2*0.35,labels=expression(0),cex=1.05,pos=1)
text(x=tau1+5, y= mu2*0.35,labels=expression(5),cex=1.05,pos=1)
text(x=tau1+10, y= mu2*0.35,labels=expression(10),cex=1.05,pos=1)
text(x=tau1+15, y= mu2*0.35,labels=expression(15),cex=1.05,pos=1)
text(x=tau1+20, y= mu2*0.35,labels=expression(20),cex=1.05,pos=1)
text(x=tau1+25, y= mu2*0.35,labels=expression(25),cex=1.05,pos=1)
text(x=tau1+30, y= mu2*0.35,labels=expression(t),cex=1.25,pos=1)

# Adding the simulations on top of the plot
matpoints(1:len1, ts.reps, type="l", col=rep("darkgrey", nsims), lty=rep(1,nsims))
points(1:len1,emp.means, type="l", lty=1, col="blue", lwd=2)

# Second plot:  theoretical weights for the weighted average of the expected value
tvec <- seq(0,20,by=0.1);
c2.t <- c2^(tvec);
onemc2.t <- 1- c2^(tvec);
plot(tvec,c2.t, type="l", col="blue", lwd=2, ylab=expression(c[2]^t), xlab=expression(t), cex.lab=1.15,main="Weights for the mean of X(t)",cex.main=1.5); 
points(tvec,onemc2.t, type="l", col="red", lwd=2); 
axis(side=4, at=c(0.0,0.2,0.4,0.6,0.8,1.0));
mtext(expression(1-c[2]^t),side=4,line=3, cex=1.15)
legend(x=10,y=0.95,legend=c(expression(c[2]^t), expression(1-c[2]^t) ), lwd=c(2,2), col=c("blue","red"),lty=c(1,1),bty="n")


# Third plot:  theoretical and empirical variances

plot(1:len1, vars.t , type="l", col="darkgrey", lwd=2, xlab="", ylab="Variance of log-abundance",ylim=c(nu2*0.45,nu1*1.5), xlim=c(-2,len1+5),cex.lab=1.5, main=expression(V(X(t+tau))==nu[1]*c[2]^{2*t} + nu[2]*(1-c[2]^{2*t})),cex.main=1.5,bty="n", xaxt="n", yaxt="n")
arrows(x0=c(0,0), y0=c(nu2*0.57, nu2*0.57), x1=c(0, len1+5), y1=c(nu1*1.5,nu2*0.57),lwd=2, code=2,length=0.1)
segments(x0=tau1+1, y0=nu2*0.57, x1=tau1+1,y1=mu1*1.25, lwd=2,col="red")
text(x=-2.5,y=nu1, labels=expression(nu[1]),cex=1.25, srt=90);
text(x=-2.5,y=nu2, labels=expression(nu[2]),cex=1.25, srt=90);
text(x=tau1+1, y= nu2*0.57,labels=expression(0),cex=1.05,pos=1)
text(x=tau1+5, y= nu2*0.57,labels=expression(5),cex=1.05,pos=1)
text(x=tau1+10, y= nu2*0.57,labels=expression(10),cex=1.05,pos=1)
text(x=tau1+15, y= nu2*0.57,labels=expression(15),cex=1.05,pos=1)
text(x=tau1+20, y= nu2*0.57,labels=expression(20),cex=1.05,pos=1)
text(x=tau1+25, y= nu2*0.57,labels=expression(25),cex=1.05,pos=1)
text(x=tau1+30, y= nu2*0.57,labels=expression(t),cex=1.25,pos=1)
# Adding the empirical estimates of the variance
points(1:len1, emp.vars, col="blue", lty=1, type="l", lwd=2)


tvec <- seq(0,20,by=0.1);
c2.2t <- c2^(2*tvec);
onemc2.2t <- 1- c2^(2*tvec);
plot(tvec,c2.2t, type="l", col="blue", lwd=2, ylab=expression(c[2]^{2*t}), xlab=expression(t), cex.lab=1.15,main="Weights for the variance of X(t)",cex.main=1.5); 
points(tvec,onemc2.2t, type="l", col="red", lwd=2); 
axis(side=4, at=c(0.0,0.2,0.4,0.6,0.8,1.0));
mtext(expression(1-c[2]^{2*t}),side=4,line=3, cex=1.15)
legend(x=10,y=0.95,legend=c(expression(c[2]^{2*t}), expression(1-c[2]^{2*t}) ), lwd=c(2,2), col=c("blue","red"),lty=c(1,1),bty="n")

dev.off()
