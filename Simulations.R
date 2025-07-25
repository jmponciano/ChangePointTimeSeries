source("BPGSSToolkit.R")

# Simulation and estimation leading to Figure 3:


####  Approximating different density-dependent models with Gompertz using as metric
####  the half life


######---- STEP 1:  Set up and sample plots of simulations under different 
######              compensatory dynamics
 
a1 <- 1.35;#3.5;
c1 <- 0.75;
a2 <- 0.25;
c2 <- (1/2)^(1/10);
mu1 <- a1/(1-c1);
mu2 <- a2/(1-c2);
print(mu1)
print(mu2)
t.half <- log(1/2)/log(abs(c2))
mid.thalf <- tau1+t.half/2
mean.mean <- (mu1+mu2)/2
 
len <- 500;#len1; # 160 above
tau1 <- 40;
qp1 <- len - tau1;
K1 <- exp(mu1);
K2 <- exp(mu2);
alpha <- 20#1000; #scale parameter of the Gamma
k     <- 4.5*alpha; # shape parameter of the Gamma
b1 <- -log(k/alpha)/K1;
b2 <- -log(k/alpha)/K2;
#X11()
short.len <- 150
par(mfrow=c(3,1))
# Undercompensatory
pop.sim1a <- negbinde.sim(no=K1,len=tau1,scalep=alpha,shapep=k,ddp=list("Below",K=K1, beta = 0.5))
pop.sim1b <- negbinde.sim(no=pop.sim1a[tau1],len=qp1,scalep=alpha,shapep=k,ddp=list("Below",K=K2, beta = 0.10))
pop.sim1 <- c(pop.sim1a,pop.sim1b)
plot(1:short.len, (pop.sim1)[1:short.len], type="l",col="red", lwd=2)


# Compensatory
pop.sim1a <- negbinde.sim(no=K1,len=tau1,scalep=alpha,shapep=k,ddp=list("Below",K=K1, beta = 0.5))
pop.sim1b <- negbinde.sim(no=pop.sim1a[tau1],len=qp1,scalep=alpha,shapep=k,ddp=list("Below",K=K2, beta = 1))
pop.sim1 <- c(pop.sim1a,pop.sim1b)
plot(1:short.len, (pop.sim1)[1:short.len], type="l",col="red", lwd=2)

# Overcompensatory
pop.sim1a <- negbinde.sim(no=K1,len=tau1,scalep=alpha,shapep=k,ddp=list("Below",K=K1, beta = 0.5))
pop.sim1b <- negbinde.sim(no=pop.sim1a[tau1],len=qp1,scalep=alpha,shapep=k,ddp=list("Below",K=K2, beta = 1.2))
pop.sim1 <- c(pop.sim1a,pop.sim1b)
plot(1:short.len, (pop.sim1)[1:short.len], type="l",col="red", lwd=2)

######---- End of STEP 1---######

######---- STEP 2: Simulate a large number of time series with a change point at time ----- ######
######  'tau1' defined above

nsims <- 1000;
len2  <- 150;
qp1   <- len2 - tau1;
undercomp <- matrix(0, nrow=nsims, ncol=len2);
comp <- matrix(0, nrow=nsims, ncol=len2);
overcomp <- matrix(0, nrow=nsims, ncol=len2);
tausq <- 0.2315


# simulating data
for(i in 1:nsims){
	
	# Undercompensatory
	pop.sim1a <- negbinde.sim(no=K1,len=tau1,scalep=alpha,shapep=k,ddp=list("Below",K=K1, beta = 0.5))
	pop.sim1b <- negbinde.sim(no=pop.sim1a[tau1],len=qp1,scalep=alpha,shapep=k,ddp=list("Below",K=K2, beta = 0.10))
	undercomp[i,] <- log(c(pop.sim1a,pop.sim1b)) + rnorm(n=len2, mean=0, sd=sqrt(tausq));

	# Compensatory
	pop.sim1a <- negbinde.sim(no=K1,len=tau1,scalep=alpha,shapep=k,ddp=list("Below",K=K1, beta = 0.5))
	pop.sim1b <- negbinde.sim(no=pop.sim1a[tau1],len=qp1,scalep=alpha,shapep=k,ddp=list("Below",K=K2, beta = 1))
	comp[i,] <- log(c(pop.sim1a,pop.sim1b)) + rnorm(n=len2, mean=0, sd=sqrt(tausq));

	# Overcompensatory
	pop.sim1a <- negbinde.sim(no=K1,len=tau1,scalep=alpha,shapep=k,ddp=list("Below",K=K1, beta = 0.5))
	pop.sim1b <- negbinde.sim(no=pop.sim1a[tau1],len=qp1,scalep=alpha,shapep=k,ddp=list("Below",K=K2, beta = 1.2))
	overcomp[i,] <- log(c(pop.sim1a,pop.sim1b)) + rnorm(n=len2, mean=0, sd=sqrt(tausq));	
}

######---- STEP 3:  Obtain the empirical estimates of the half life 

mean.undercomp <- apply(undercomp,2,mean);
mean.comp <- apply(comp,2,mean);
mean.overcomp <- apply(overcomp,2,mean);
mean.mean.sim <- (mu1+mu2)/2;
delta.ucomp <- mean.undercomp-mean.mean.sim;
delta.comp <- mean.comp-mean.mean.sim;
delta.ocomp <- mean.overcomp-mean.mean.sim;
x1.uc <- (which(delta.ucomp<0, arr.ind=TRUE))[1];
x0.uc <- x1.uc-1; 
y1.uc <- mean.undercomp[x1.uc];y0.uc <- mean.undercomp[x0.uc];
thalf.uc <- x0.uc + (mean.mean.sim-y0.uc)/(y1.uc-y0.uc);

x1.c <- (which(delta.comp<0, arr.ind=TRUE))[1];
x0.c <- x1.c-1;
y1.c <- mean.comp[x1.c];y0.c <- mean.comp[x0.c];
thalf.c <- x0.c + (mean.mean.sim-y0.c)/(y1.c-y0.c);

x1.oc <- (which(delta.ocomp<0, arr.ind=TRUE))[1];
x0.oc <- x1.oc -1;
y1.oc <- mean.overcomp[x1.oc];y0.oc <- mean.overcomp[x0.oc];
thalf.oc <- x0.oc + (mean.mean.sim-y0.oc)/(y1.oc-y0.oc);

######---- STEP 4: Fitting the Gompertz Breakpoint model to each of the 1000 time series simulated above, for each
#####              type of compensatory dynamics (uc=under-compensatory, c=compensatory, oc=overcompensatory)
#####              Beware, this step will take hours to run (depending of course on how big 'B' is)!!!

B <- 2; #nsims;
# Store mles, two means and mles of half-life
mlesmat.uc <- matrix(0,nrow=B,ncol=10);
mlesmat.c <- matrix(0,nrow=B, ncol=10);
mlesmat.oc <- matrix(0,nrow=B, ncol=10);

# Now doing the estimation on many time series
# On purpose, I will do the loops for each case (uc, c, oc) SEPARATELY!

# Specify initial parameter values again...
guess8 <- c(log(0.3929), log(0.3929), atanh(0.7934),atanh(0.7934), log(0.09725), log(0.09725), log(0.2315));

# Undercompensatory dynamics bootstrap
for(i in 1:B){
	
	lSetophaga <- undercomp[i,]
	early.stop <- is.infinite(lSetophaga)
	yes.no <- sum(early.stop)
	if(yes.no>0){index <- which(early.stop==T,arr.ind=T)[1]; lSetophaga <- lSetophaga[1:(index-1)]}
	
	trial.optim <- optim(par=guess8, fn=model.profmle, method="Nelder-Mead", data.vec=lSetophaga, M=tau1, model.flag="m8")
	mles <- c(exp(trial.optim$par[1:2]), tanh(trial.optim$par[3:4]), exp(trial.optim$par[5:7]))
	mu1.hat <- mles[1]/(1-mles[3])
	mu2.hat <- mles[2]/(1-mles[4])
	t.half.hat <- log(0.5)/log(abs(mles[4]))
	mlesmat.uc[i,1:7] <- mles
	mlesmat.uc[i,8] <- mu1.hat
	mlesmat.uc[i,9] <- mu2.hat
	mlesmat.uc[i,10] <- tau1+t.half.hat
	
}


#Compensatory dynamics boostrap
for(i in 1:B){
	
	lSetophaga <- comp[i,]
	early.stop <- is.infinite(lSetophaga)
	yes.no <- sum(early.stop)
	if(yes.no>0){index <- which(early.stop==T,arr.ind=T)[1]; lSetophaga <- lSetophaga[1:(index-1)]}

	trial.optim <- optim(par=guess8, fn=model.profmle, method="Nelder-Mead", data.vec=lSetophaga, M=tau1, model.flag="m8")
	mles <- c(exp(trial.optim$par[1:2]), tanh(trial.optim$par[3:4]), exp(trial.optim$par[5:7]))
	mu1.hat <- mles[1]/(1-mles[3])
	mu2.hat <- mles[2]/(1-mles[4])
	t.half.hat <- log(0.5)/log(abs(mles[4]))
	mlesmat.c[i,1:7] <- mles
	mlesmat.c[i,8] <- mu1.hat
	mlesmat.c[i,9] <- mu2.hat
	mlesmat.c[i,10] <- tau1+t.half.hat
	
}


#Overcompensatory dynamics bootstrap
for(i in 1:B){
	
	lSetophaga <- overcomp[i,]
	early.stop <- is.infinite(lSetophaga)
	yes.no <- sum(early.stop)
	if(yes.no>0){index <- which(early.stop==T,arr.ind=T)[1]; lSetophaga <- lSetophaga[1:(index-1)]}

	trial.optim <- optim(par=guess8, fn=model.profmle, method="Nelder-Mead", data.vec=lSetophaga, M=tau1, model.flag="m8")
	mles <- c(exp(trial.optim$par[1:2]), tanh(trial.optim$par[3:4]), exp(trial.optim$par[5:7]))
	mu1.hat <- mles[1]/(1-mles[3])
	mu2.hat <- mles[2]/(1-mles[4])
	t.half.hat <- log(0.5)/log(abs(mles[4]))
	mlesmat.oc[i,1:7] <- mles
	mlesmat.oc[i,8] <- mu1.hat
	mlesmat.oc[i,9] <- mu2.hat
	mlesmat.oc[i,10] <- tau1+t.half.hat
	
}

