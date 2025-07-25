# Figure 1 in the manuscript:

########## Functions used saved in BPGSSToolkit.R ##########
source("BPGSSToolkit.R")

# What is the half-life of the density dependent effect?


len1 <- 160;
tau1 <- 60;
a1 <- 3.5;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.30;
c2 <- (1/2)^(1/10);
mu1 <- a1/(1-c1);
mu2 <- a2/(1-c2);
t.half <- log(1/2)/log(abs(c2))
mid.thalf <- tau1+t.half/2
mean.mean <- (mu1+mu2)/2

lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1);
mean.t <- lowleft.m$mu.vec;
var.t <- lowleft.m$Sigma;
ts.reps <- my.rmvn(n=5, mu.vec=mean.t, cov.mat = var.t)

# Plot of 5 time series replicates stopping at the breakpoint
par(bty="n", xaxt="n", yaxt="n", oma=c(1,1,1,2))
matplot(1:tau1, ts.reps[1:tau1,], type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3), xlim=c(-5,len1+30), xlab="Time", ylab="Log population abundance", cex.lab=1.5)
arrows(x0=c(0,0),y0=c(mu2-3,mu2-3),x1=c(0,len1+30),y1=c(1.05*max(ts.reps),mu2-3),lwd=2, code=2,length=0.1)


# Dotted line showing mean of first statio pdf
points(1:tau1, rep(mu1,tau1), type="l", lty=2, lwd=2)
text(x=-6,y=mu1, labels=expression(mu[1]),cex=1.25, srt=90)

# pdf of the statio distrib before the break point
sd1 <- sqrt(sigmasq1/(1-c1^2));
supp1 <- seq(mu1-3*sd1, mu1+3*sd1, by=0.01);
dens1 <- 10*dnorm(supp1,mean=mu1,sd=sd1);
segments(x0=tau1,y0=mu1-3*sd1,x1=tau1,y1=mu1+3*sd1,col="black",lwd=1);
points(tau1+dens1,supp1,type="l",lty=1,lwd=2,col="red")

#plot of the rest of the time series: from the breakpoint onwards
matplot(1:len1, ts.reps, type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3), xlim=c(-5,len1+30), xlab="Time", ylab="Log population abundance", cex.lab=1.5, add=TRUE)

# Dotted line showing mean of second statio pdf
points(1:len1, rep(mu2,len1), type="l", lty=2, lwd=2)
text(x=-6,y=mu2, labels=expression(mu[2]),cex=1.25, srt=90)

# pdf of the statio distrib after the break point
sd2 <- sqrt(sigmasq2/(1-c2^2));
supp2 <- seq(mu2-3*sd2, mu2+3*sd2, by=0.01);
dens2 <- 10*dnorm(supp2,mean=mu2,sd=sd2);
segments(x0=len1,y0=mu2-3*sd2,x1=len1,y1=mu2+3*sd2,col="black",lwd=1);
points(len1+dens2,supp2,type="l",lty=1,lwd=2,col="red")

# Dotted line showing mean of means
points(1:(tau1+t.half), rep(mean.mean,length(1:(tau1+t.half))), type="l", lty=2, lwd=2)
text(x=-6,y=mean.mean, labels=expression(bar(mu)),cex=1.25, srt=90)
#text(x=-6,y=mean.mean, labels=expression((mu[1] + mu[2])/2),cex=1.25, srt=90)


# pdf of the distrib. at the half life
mu3 <- mean.mean
sd3 <- sqrt((sd1^2)*c2^(2*t.half) + (sd2^2)*(1-(c2^(2*t.half))) )
supp3 <- seq(mu3-3*sd3, mu3+3*sd3, by=0.01);
dens3 <- 10*dnorm(supp3,mean=mu3,sd=sd3);
segments(x0=tau1+t.half,y0=mu3-3*sd3,x1=tau1+t.half,y1=mu3+3*sd3,col="black",lwd=1);
points(tau1+t.half+dens3,supp3,type="l",lty=1,lwd=2,col="red")

# Curly bracket showing halflife span
aa<-CurlyBraces(x=0,y=0 , mrange=t.half, direction = 2,my.lwd=2 )
points(aa$b_sequence+mid.thalf,aa$a_sequence+mean.mean,type="l", lwd=2,col="black")

# Arrow pointing to curly bracket
yz <- min(aa$a_sequence+mean.mean);
arrows(x0=tau1-7,y0=yz-1.3,x1=mid.thalf,y1=yz*0.98,code=2, length=0.07,lwd=2)

# Half life label
text(x=tau1-7, y=yz-1.75, labels="Half life", cex=1.25)




##### Same figure with two panels:
postscript("Fig1.eps", horizontal=FALSE, paper="letter")
par(mar=c(1,1,1,1),oma=c(3,3,0.5,0.5),mgp=c(1.25,0.5,0))
################################################## Panel 1 ##################################
len1 <- 160;
tau1 <- 60;
a1 <- 3.5;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.30;
c2 <- (1/2)^(1/10);
mu1 <- a1/(1-c1);
mu2 <- a2/(1-c2);
t.half <- log(1/2)/log(abs(c2))
mid.thalf <- tau1+t.half/2
mean.mean <- (mu1+mu2)/2

lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1);
mean.t <- lowleft.m$mu.vec;
var.t <- lowleft.m$Sigma;
ts.reps <- my.rmvn(n=5, mu.vec=mean.t, cov.mat = var.t)


# Plot of 5 time series replicates stopping at the breakpoint
par(mfrow=c(2,1),bty="n", xaxt="n", yaxt="n",oma=c(1,1,1,2))
matplot(1:tau1, ts.reps[1:tau1,], type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3), xlim=c(-5,len1+30), xlab="", ylab="", cex.lab=1.5)
arrows(x0=c(0,0),y0=c(mu2-3,mu2-3),x1=c(0,len1+30),y1=c(1.05*max(ts.reps),mu2-3),lwd=2, code=2,length=0.1)


# Dotted line showing mean of first statio pdf
points(1:tau1, rep(mu1,tau1), type="l", lty=2, lwd=2)
text(x=-6,y=mu1, labels=expression(mu[1]),cex=1.25, srt=90)

# pdf of the statio distrib before the break point
sd1 <- sqrt(sigmasq1/(1-c1^2));
supp1 <- seq(mu1-3*sd1, mu1+3*sd1, by=0.01);
dens1 <- 10*dnorm(supp1,mean=mu1,sd=sd1);
segments(x0=tau1,y0=mu1-3*sd1,x1=tau1,y1=mu1+3*sd1,col="black",lwd=1);
points(tau1+dens1,supp1,type="l",lty=1,lwd=2,col="black")

#plot of the rest of the time series: from the breakpoint onwards
matplot(1:len1, ts.reps, type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu2-3, mu1+3), xlim=c(-5,len1+30), xlab="Time", ylab="Log population abundance", cex.lab=1.5, add=TRUE)

# Dotted line showing mean of second statio pdf
points(1:len1, rep(mu2,len1), type="l", lty=2, lwd=2)
text(x=-6,y=mu2, labels=expression(mu[2]),cex=1.25, srt=90)

# pdf of the statio distrib after the break point
sd2 <- sqrt(sigmasq2/(1-c2^2));
supp2 <- seq(mu2-3*sd2, mu2+3*sd2, by=0.01);
dens2 <- 10*dnorm(supp2,mean=mu2,sd=sd2);
segments(x0=len1,y0=mu2-3*sd2,x1=len1,y1=mu2+3*sd2,col="black",lwd=1);
points(len1+dens2,supp2,type="l",lty=1,lwd=2,col="black")

# Dotted line showing mean of means
points(1:(tau1+t.half), rep(mean.mean,length(1:(tau1+t.half))), type="l", lty=2, lwd=2)
text(x=-6,y=mean.mean, labels=expression(bar(mu)),cex=1.25, srt=90)
#text(x=-6,y=mean.mean, labels=expression((mu[1] + mu[2])/2),cex=1.25, srt=90)


# pdf of the distrib. at the half life
mu3 <- mean.mean
sd3 <- sqrt((sd1^2)*c2^(2*t.half) + (sd2^2)*(1-(c2^(2*t.half))) )
supp3 <- seq(mu3-3*sd3, mu3+3*sd3, by=0.01);
dens3 <- 10*dnorm(supp3,mean=mu3,sd=sd3);
segments(x0=tau1+t.half,y0=mu3-3*sd3,x1=tau1+t.half,y1=mu3+3*sd3,col="black",lwd=1);
points(tau1+t.half+dens3,supp3,type="l",lty=1,lwd=2,col="black")

# Curly bracket showing halflife span
aa<-CurlyBraces(x=0,y=0 , mrange=t.half, direction = 2,my.lwd=2 )
points(aa$b_sequence+mid.thalf,aa$a_sequence+mean.mean,type="l", lwd=2,col="red")

# Arrow pointing to curly bracket
yz <- min(aa$a_sequence+mean.mean);
arrows(x0=tau1-7,y0=yz-1.3,x1=mid.thalf,y1=yz*0.98,code=2, length=0.07,lwd=2)

# Half life label
text(x=tau1-7, y=yz-1.75, labels="Half life", cex=1.25)

############################################  Panel 2  #############################

len1 <- 160;
tau1 <- 60;
a2 <- 3.5;
c2 <- 0.88#0.75;
sigmasq1 <- sigmasq2 <- 0.12#0.09725;
a1 <- 0.30;
c1 <- (1/2)^(1/10);
mu1 <- a1/(1-c1);
mu2 <- a2/(1-c2);
t.half <- log(1/2)/log(abs(c2))
mid.thalf <- tau1+t.half/2
mean.mean <- (mu1+mu2)/2

lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1);
mean.t <- lowleft.m$mu.vec;
var.t <- lowleft.m$Sigma;
ts.reps <- my.rmvn(n=5, mu.vec=mean.t, cov.mat = var.t)

# Plot of 5 time series replicates stopping at the breakpoint
#par(bty="n", xaxt="n", yaxt="n", oma=c(1,1,1,2))
matplot(1:tau1, ts.reps[1:tau1,], type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu1-3, mu2+3), xlim=c(-5,len1+30), xlab="", ylab="", cex.lab=1.5)
arrows(x0=c(0,0),y0=c(mu1-3,mu1-3),x1=c(0,len1+30),y1=c(1.05*max(ts.reps),mu1-3),lwd=2, code=2,length=0.1)


# Dotted line showing mean of first statio pdf
points(1:tau1, rep(mu1,tau1), type="l", lty=2, lwd=2)
text(x=-6,y=mu1, labels=expression(mu[1]),cex=1.25, srt=90)

# pdf of the statio distrib before the break point
sd1 <- sqrt(sigmasq1/(1-c1^2));
supp1 <- seq(mu1-3*sd1, mu1+3*sd1, by=0.01);
dens1 <- 10*dnorm(supp1,mean=mu1,sd=sd1);
segments(x0=tau1,y0=mu1-3*sd1,x1=tau1,y1=mu1+3*sd1,col="black",lwd=1);
points(tau1+dens1,supp1,type="l",lty=1,lwd=2,col="black")

#plot of the rest of the time series: from the breakpoint onwards
matplot(1:len1, ts.reps, type="l", col=rep("darkgrey",5), lty=rep(1,5), add=TRUE)

# Dotted line showing mean of second statio pdf
points(1:len1, rep(mu2,len1), type="l", lty=2, lwd=2)
text(x=-6,y=mu2, labels=expression(mu[2]),cex=1.25, srt=90)

# pdf of the statio distrib after the break point
sd2 <- sqrt(sigmasq2/(1-c2^2));
supp2 <- seq(mu2-3*sd2, mu2+3*sd2, by=0.01);
dens2 <- 10*dnorm(supp2,mean=mu2,sd=sd2);
segments(x0=len1,y0=mu2-3*sd2,x1=len1,y1=mu2+3*sd2,col="black",lwd=1);
points(len1+dens2,supp2,type="l",lty=1,lwd=2,col="black")

# Dotted line showing mean of means
points(1:(tau1+t.half), rep(mean.mean,length(1:(tau1+t.half))), type="l", lty=2, lwd=2)
text(x=-6,y=mean.mean, labels=expression(bar(mu)),cex=1.25, srt=90)
#text(x=-6,y=mean.mean, labels=expression((mu[1] + mu[2])/2),cex=1.25, srt=90)


# pdf of the distrib. at the half life
mu3 <- mean.mean
sd3 <- sqrt((sd1^2)*c2^(2*t.half) + (sd2^2)*(1-(c2^(2*t.half))) )
supp3 <- seq(mu3-3*sd3, mu3+3*sd3, by=0.01);
dens3 <- 10*dnorm(supp3,mean=mu3,sd=sd3);
segments(x0=tau1+t.half,y0=mu3-3*sd3,x1=tau1+t.half,y1=mu3+3*sd3,col="black",lwd=1);
points(tau1+t.half+dens3,supp3,type="l",lty=1,lwd=2,col="black")

# Curly bracket showing halflife span
aa<-CurlyBraces(x=0,y=0 , mrange=t.half, direction = 1,my.lwd=2 )
points(aa$b_sequence+mid.thalf,aa$a_sequence+mean.mean+0.6,type="l", lwd=2,col="red")

# Arrow pointing to curly bracket
yz <- max(aa$a_sequence+mean.mean);
arrows(x0=tau1-7,y0=yz+1.7,x1=mid.thalf,y1=mean.mean+1.5,code=2, length=0.07,lwd=2)

# Half life label
text(x=tau1-10, y=yz+4, labels="Half life", cex=1.25)
mtext("Time",side=1,line=-0.75,outer=TRUE,cex=1.5)
mtext("Log population abundance",side=2,line=-0.5,outer=TRUE,cex=1.5)
dev.off()







