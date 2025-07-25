## Figure 2 in the manuscript:

source("BPGSSToolkit.R")


postscript("Fig2.eps", horizontal=TRUE, paper="letter")
par(mar=c(3,3,1,3),oma=c(0.5,0.5,0.5,0.5),mgp=c(1.75,0.75,0))

c2.vals.left <- seq(-0.92,-0.001,by=0.001);
c2.vals.right <- seq(0.001,0.92,by=0.001);
half.life.left <- log(1/2)/log(abs(c2.vals.left));
half.life.right <- log(1/2)/log(abs(c2.vals.right));
par(ps=12,bty="n", oma=c(1,0,1,2))
plot(c2.vals.left,half.life.left, type="l", col="red", lwd=2, xlab=expression(c[2]), ylab="", xlim=c(-0.99,1.3), xaxt="n",yaxt="n", tck=0, cex.lab=1.5,asp=-0.5)
axis(1, cex.axis=1.25, tck=NA)
axis(4, cex.axis=1.25, tck=NA, at=c(0,1,2,4,6,8))
mtext("Half Life",side=4,line=3, cex=1.5)
points(c2.vals.right, half.life.right, type="l", col="red", lwd=2)
segments(x0=-0.5,y0=0,x1=-0.5,y1=1, lty=2)
segments(x0=-0.5,y0=1,x1=0.5,y1=1, lty=2)
segments(x0=0.5,y0=1,x1=0.5,y1=0, lty=2)

# Low left corner: c2 = -0.35
par(mar=c(3,3,1,3), tcl=-0.25, mgp=c(1.75,0.5,0))
len1 <- 31;
tau1 <- 16;
a1 <- 0.3929;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.39;
c2 <- -0.35;
lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1)
ts.reps <- my.rmvn(n=5, mu.vec=lowleft.m$mu.vec, cov.mat = lowleft.m$Sigma)
par(plt=c(0.03,0.23,0.22,0.42), new=T, las=1, ps=9,bty="n", xaxt="s", yaxt="s", oma=c(0.5,0.5,0.75,1))
matplot(1:len1, y=ts.reps, type="l", xlab="", ylab="",ylim=c(-1,2.75))
points(1:len1, lowleft.m$mu.vec  , type="l", lwd=2, xlab="", ylab="")
mtext(expression(c[2]==-0.35), side=3,adj=0, cex=1.25)

# Top left corner: c2 = -0.75
len1 <- 31;
tau1 <- 16;
a1 <- 0.3929;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.39;
c2 <- -0.75;
lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1)
ts.reps <- my.rmvn(n=5, mu.vec=lowleft.m$mu.vec, cov.mat = lowleft.m$Sigma)
par(plt=c(0.20,0.40,0.7,0.9), new=T, las=1, ps=9,bty="n", xaxt="s", yaxt="s", oma=c(0.5,0.5,0.75,1))
matplot(1:len1, y=ts.reps, type="l", xlab="", ylab="", main= expression(c[2]==-0.75),ylim=c(-1,2.75))
points(1:len1, lowleft.m$mu.vec  , type="l", lwd=2, xlab="", ylab="")

# Low Middle: c2 = 0
len1 <- 31;
tau1 <- 16;
a1 <- 0.3929;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.39;
c2 <- 0;
lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1)
ts.reps <- my.rmvn(n=5, mu.vec=lowleft.m$mu.vec, cov.mat = lowleft.m$Sigma)
par(plt=c(0.38,0.58,0.35,0.55), new=T, las=1, ps=9,bty="n", xaxt="s", yaxt="s", oma=c(0.5,0.5,0.75,1))
matplot(1:len1, y=ts.reps, type="l", xlab="", ylab="", main= expression(c[2]==0),ylim=c(-1,2.75))
points(1:len1, lowleft.m$mu.vec  , type="l", lwd=2, xlab="", ylab="")

# Right middle: c2 = 0.35
len1 <- 31;
tau1 <- 16;
a1 <- 0.3929;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.39;
c2 <- 0.35;
lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1)
ts.reps <- my.rmvn(n=5, mu.vec=lowleft.m$mu.vec, cov.mat = lowleft.m$Sigma)
par(plt=c(0.7,0.9,0.22,0.42), new=T, las=1, ps=9,bty="n", xaxt="s", yaxt="s", oma=c(0.5,0.5,0.75,1))
matplot(1:len1, y=ts.reps, type="l", xlab="", ylab="", main= expression(c[2]==0.35),ylim=c(-1,2.75))
points(1:len1, lowleft.m$mu.vec  , type="l", lwd=2, xlab="", ylab="")

# Top Middle: c2 = 0.75
len1 <- 31;
tau1 <- 16;
a1 <- 0.3929;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.09725;
a2 <- 0.1;
c2 <- 0.75;
lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1)
ts.reps <- my.rmvn(n=5, mu.vec=lowleft.m$mu.vec, cov.mat = lowleft.m$Sigma)
par(plt=c(0.57,0.77,0.7,0.9), new=T, las=1, ps=9,bty="n", xaxt="s", yaxt="s", oma=c(0.5,0.5,0.75,1))
matplot(1:len1, y=ts.reps, type="l", xlab="", ylab="", main= expression(c[2]==0.75),ylim=c(-1,2.75))
points(1:len1, lowleft.m$mu.vec  , type="l", lwd=2, xlab="", ylab="")

dev.off()
