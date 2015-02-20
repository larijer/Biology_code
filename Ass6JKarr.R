######-----Jered Karr-----#######
######-----Created November 14 2012-----#######
######-----Created for Brett Melbourn's class EBIO 5460-----#######


#------Function definitions--------------





####### Bite Mass Model ########
# R = rate of processing food (grams per minute)
# S = bite mass (in grams)
# h = time required to crop a bite ( in minutes )

BMMo <- function (R,S,h){
	I <- ( R*S )/( R*h + S )
	return(I)
}



####### Density Model ########
# V = Velocity of travel when not cropping (meters per minute)
# D = density of plants (number per square meter)
# S = bite mass (in grams)
# h = time required to crop a bite ( in minutes )

DMo <- function(V,D,S,h,a){
	I <- (V * a * sqrt(D) * S) / (1 + a * h * V * sqrt(D))
	return(I)
}


##############   Function  ##################
nll_BMMo_norm <- function(p,trueS,trueI){
	R <- p[1]
	h <- exp(p[2])
	sd <- exp(p[3])  #the standard deviation of the normal
	S <- trueS
	I <- BMMo(R,S,h)
	nll <- -sum(dnorm(trueI,mean=I,sd=sd,log=TRUE))
	return(nll)
}

nll_DMo_norm <- function(p,trueS,trueD, a, trueI){
	V <- p[1]
	h <- exp(p[2])
	sd <- exp(p[3])  #the standard deviation of the normal
	a <- a
	D <- trueD
	S <- trueS
	I <- DMo(V,D,S,h,a)
	nll <- -sum(dnorm(trueI,mean=I,sd=sd,log=TRUE)) 
	return(nll)
}


#----nll_DMo_norm_V() -----------------------------
#A variant of nll_DMo_norm that computes the likelihood
#for a given value of V, while optimizing the other
#parameters.
#Notice that we now pass V as a fixed argument - it is not
#in the vector of parameters (p) to be optimized.
# Then corresponding functions for all the parameters


nll_DMo_norm_V <- function(p,V,trueS,trueD,a,trueI){
	h <- exp(p[1])
	sd <- exp(p[2])  #the standard deviation of the normal
	a <- a
	D <- trueD
	S <- trueS
	I <- DMo(V,D,S,h,a)
	nll <- -sum(dnorm(trueI,mean=I,sd=sd,log=TRUE)) 
	return(nll)
}
nll_DMo_norm_h <- function(p,h,trueS,trueD,a,trueI){
	V <- exp(p[1])
	sd <- exp(p[2])  #the standard deviation of the normal
	a <- a
	D <- trueD
	S <- trueS
	I <- DMo(V,D,S,h,a)
	nll <- -sum(dnorm(trueI,mean=I,sd=sd,log=TRUE)) 
	return(nll)
}

nll_DMo_norm_sd <- function(p,sd,trueS,trueD,a,trueI){
	V <- exp(p[1])
	h <- exp(p[2])
	a <- a
	D <- trueD
	S <- trueS
	I <- DMo(V,D,S,h,a)
	nll <- -sum(dnorm(trueI,mean=I,sd=sd,log=TRUE)) 
	return(nll)
}
#----Variance of the Negative log likelihood
# This is used for analysis of point deletion
# Only for density model because is best fit for my data
DMo_nll.var <- function(V,D,S,h,a,sd) {
	mu <- DMo(V,D,S,h,a)
	return( mu + mu*mu/sd )
}

realdata <- read.csv('student14_ass6.csv')
trueS <- realdata$S  #  plant size (grams)
trueD <- realdata$D  #  plant density (number plants m^-2)
trueI <- realdata$I  #  intake (grams min^-1)

plot(trueS, trueI)

#----Parameters------
R <- 48 		#the rate of processing of food in the mouth as Rmax (grams per minute)
h <- .02	#time required to crop a bite as h (in minutes)
V <- 64.8	#Vmax (meters per minute) as the velocity of travel in the absence of cropping
sd <- 2
a <- 1		# a=2 when pattern is random, a =1 when the pattern is uniform;

plot(log(trueS),log(trueD)) #Plot to see if distribuation (a) is uniform (yes it is)

#-------Fit the model to the data by Direct Search
R_range <- (seq(1000,2000,length.out=100))
h_range <- seq(.01,.1,length.out=100) 
sd_range <- seq(5,20,length.out=100)

#Initialize storage objects (matrix, vector, and scalar)
n <- length(R_range)*length(h_range)*length(sd_range)
results_BMMo <- matrix(NA,n,4)  #columns: R, h, sd, nll
colnames(results) <- c("R","h","sd","nll")
p_best_BMMo <- NA
nll_best_BMMo <- Inf

#Direct search over the three parameters: R, h, sd for BMMo
j <- 1 #Row index for results matrix
for (R in R_range) {
	for (h in h_range) {
		for (sd in sd_range){

				parm <- c(R, log(h),log(sd))
				nll <- nll_BMMo_norm(parm, trueS, trueI)
	
	#	Keep the results
				results_BMMo[j,] <- c(R, h, sd, nll)
				j <- j + 1

	#	Record current best solution
		if (nll < nll_best_BMMo) {
			nll_best_BMMo <- nll
			p_best_BMMo <- cbind(R, h, sd) #cbind automatically gives labels
			}
		}
	}
#	Monitor progress
	print( paste(round(100*j/n),"%",sep=""), quote=FALSE )
}

#Print the minimum sum of squares and best parameter values
nll_best_BMMo
p_best_BMMo


#-------Fit the model to the data by Direct Search
V_range <- (seq(1,100,length.out=100))
h_range <- seq(.01,.1,length.out=100) 
sd_range <- seq(.1,3,length.out=100)

#Initialize storage objects (matrix, vector, and scalar)
n <- length(V_range)*length(h_range)*length(sd_range)
results_DMo <- matrix(NA,n,4)  #columns: V, h, sd, nll
colnames(results) <- c("V","h","sd","nll")
p_best_DMo <- NA
nll_best_DMo <- Inf

#Direct search over the three parameters: V, h, sd for DMo
j <- 1 #Row index for results matrix
for (V in V_range) {
	for (h in h_range) {
		for (sd in sd_range){

				paro <- c(V, log(h),log(sd))
				nll <- nll_DMo_norm(paro,trueS,trueD, a, trueI)
	
	#	Keep the results
				results_DMo[j,] <- c(V, h, sd, nll)
				j <- j + 1

	#	Record current best solution
		if (nll < nll_best_DMo) {
			nll_best_DMo <- nll
			p_best_DMo <- cbind(V, h, sd) #cbind automatically gives labels
			}
		}
	}
#	Monitor progress
	print( paste(round(100*j/n),"%",sep=""), quote=FALSE )
}

#Print the minimum sum of squares and best parameter values
nll_best_DMo
p_best_DMo
max(results_DMo[,4])

#----Plot sum of squares profiles for BMMo

quartz(width=10,height=5)
par(mfrow=c(1,3))
scale <- 1.5  #This is a zoom setting. Must be > 1. Smaller is zooming in.
plot(results_BMMo[,1],(results_BMMo[,4]),xlab="R",ylab="Neg Log Lik",ylim=c((nll_best_BMMo),(scale*nll_best_BMMo)), col="blue")
plot(results_BMMo[,2],(results_BMMo[,4]),xlab="h",ylab="Neg Log Lik",ylim=c((nll_best_BMMo),(scale*nll_best_BMMo)), col="blue")
plot(results_BMMo[,3],(results_BMMo[,4]),xlab="sd",ylab="Neg Log Lik",ylim=c((nll_best_BMMo),(scale*nll_best_BMMo)), col="blue")
mtext("Neg Log Likelihood profiles - direct search for Bite Mass Model",outer=TRUE,line=-2.5)


#----Plot sum of squares profiles for DMo
quartz(width=10,height=5)
par(mfrow=c(1,3))
scale <- 1.05  #This is a zoom setting. Must be > 1. Smaller is zooming in.
plot(results_DMo[,1],(results_DMo[,4]),xlab="V",ylab="Neg Log Lik",ylim=c((nll_best_DMo),(scale*nll_best_DMo)), col="blue")
plot(results_DMo[,2],(results_DMo[,4]),xlab="h",ylab="Neg Log Lik",ylim=c((nll_best_DMo),(scale*nll_best_DMo)), col="blue")
plot(results_DMo[,3],(results_DMo[,4]),xlab="sd",ylab="Neg Log Lik",ylim=c((nll_best_DMo),(scale*nll_best_DMo)), col="blue")
mtext("Neg Log Likelihood profiles - direct search of Density Model",outer=TRUE,line=-2.5)


R <- p_best_BMMo[1]
h <- p_best_BMMo[2]
sd <- p_best_BMMo[3]


#Optimization: finding the minimum negative log likelihood
parm <- c(R, log(h),log(sd))
fit_BMMo <- optim( parm, nll_BMMo_norm, trueS=trueS,trueI=trueI)
fit_BMMo		#Note convergence
fit_BMMo$par	#The parameter estimates
fit_BMMo$value	#The negative log likelihood
Bitemasschop <- exp(fit_BMMo$par[2])    # correct for negative values of h           
Bitemassstd <- exp(fit_BMMo$par[3])		# correct for negative values of sd
Bitemassrate <- fit_BMMo$par[1]
  
             
parD <- c(V, log(h), log(sd))
fit_DMo <- optim( parD, nll_DMo_norm, trueD=trueD, trueS=trueS, , a=a, trueI=trueI)
fit_DMo		#Note convergence
fit_DMo$par	#The parameter estimates
fit_DMo$value	#The negative log likelihood  
optimVDMo <- fit_DMo$par[1]          
optimhDMo <- exp(fit_DMo$par[2])    # correct for negative values of h           
optimsdDMo <- exp(fit_DMo$par[3]) 	# correct for negative values of sd

#Calculate for I using models using caluculated parameters.

optimumBMMo <- BMMo(Bitemassrate, trueS, 	Bitemasschop)  # calcualte I using a opitmized parameter for BMM
DMO_I<-DMo(optimVDMo, trueD, trueS, optimhDMo, 1) # Calcuate I using optimized parameters for DM


#AIC
# Calculatng AIC for each model

AIC_DMo <- 2*(fit_DMo$value) + 2*3
AIC_BMMo <- 2*(fit_BMMo$value) + 2*3

deltaAIC2 <- -(AIC_DMo-AIC_BMMo)

AICc <- AIC_DMo + (2*3*(3 + 1))/(length(trueI)-3-1)
AICcBMM <- AIC_BMMo + (2*3*(3 + 1))/(length(trueI)-3-1)
AIC <- rbind(AIC_DMo, AIC_BMMo)
loglik <- rbind(fit_DMo$value,fit_BMMo$value)
deltaAIC <- rbind(0, deltaAIC2)
cbind(loglik, AIC, deltaAIC)


#Akaike weights
AKB <- exp(-0.5*deltaAIC2) 
AKD <- exp(-0.5*0)
AKDD<-AKD/(AKB+AKD) 
AKBB<-AKB/(AKB+AKD) 


#Plot the data and fit
	#Graph using optim fot BMM
plot(trueS, optimumBMMo, type="l", main = "Bite Mass Model", xlab = "Bite Mass (in grams)", ylab = "Intake (grams min^-1)")
points(trueS,trueI, col="red")
segments(trueS,optimumBMMo,trueS,trueI,col="green")

par(mfrow=c(1,2))
plot(trueI, optimumBMMo, main = "Bite Mass Model", xlab = "Intake (grams min^-1)", ylab = "Modeled Intake (grams min^-1)") #Graph using optim fot DM
abline(a= 0, b=1, col="red", lw=2)



plot(trueI, optimDMo, main = "Density Model", xlab = "Intake (grams min^-1)", ylab = "Modeled Intake (grams min^-1)") #Graph using optim fot DM
abline(a= 0, b=1, col="red", lw=2)
points(trueI, optimDMo)

#Confidence intervals
#confidence interval for V for Density Model

nll.V <- rep(NA,50)
V_range <- seq(50,65,length.out=length(nll.V))
par <- c(log(optimhDMo),log(optimsdDMo))  #starting values for h, sd
i <- 1
for (V in V_range){
	nll.V[i] <- optim( par, nll_DMo_norm_V, V=V, trueS=trueS, trueD=trueD, a=a, trueI=trueI )$value
	print(i) #Monitor progress
	i <- i + 1
}

#Plot the profile
plot(V_range,nll.V,xlab="V",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for V")
axis(1,at=seq(min(V_range),max(V_range),by=2),tcl=-0.25,labels=FALSE)


#Add MLE to graph
min_nll <-  fit_DMo$value
#First check if we find a slightly better minimum than optim() from profiling
if ( min(nll.V) < fit_DMo$value ){
	optimVDMo <- V_range[which.min(nll.V)]
	min_nll <- min(nll.V)
}
abline(h=min_nll,col="red")
text(optimVDMo,min_nll,signif(optimVDMo,4),pos=3)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- min(fit_DMo$value,min(nll.V)) #We might find a slightly better fit by profiling
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2
abline(h=conf_lim,col="red")

#Calculate CI by interpolation
nll_lo <- nll.V[1:which.min(nll.V)]
V_lo <- V_range[1:which.min(nll.V)]
conf_lo <- approx(nll_lo,V_lo,xout=conf_lim)$y
nll_hi <- nll.V[which.min(nll.V):length(nll.V)]
V_hi <- V_range[which.min(nll.V):length(nll.V)]
conf_hi <- approx(nll_hi,V_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=conf_lo,col="red",lty=2)
abline(v=conf_hi,col="red",lty=2)
text(conf_lo,min_nll+4,signif(conf_lo,4),pos=4)
text(conf_hi,min_nll+4,signif(conf_hi,4),pos=2)

#The confidence interval
cbind(conf_lo,conf_hi)

#confidence interval for h for Density Model

nll.h <- rep(NA,50)
h_range <- seq(.03,.04,length.out=length(nll.h))
par <- c(log(optimVDMo),log(optimsdDMo))  #starting values for V, sd
i <- 1
for (h in h_range){
	nll.h[i] <- optim( par, nll_DMo_norm_h, h=h, trueS=trueS, trueD=trueD, a=a, trueI=trueI )$value
	print(i) #Monitor progress
	i <- i + 1
}

#Plot the profile
plot(h_range,nll.h,xlab="h",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for h")
axis(1,at=seq(min(h_range),max(h_range),by=2),tcl=-0.25,labels=FALSE)


#Add MLE to graph
min_nll <-  fit_DMo$value
#First check if we find a slightly better minimum than optim() from profiling
if ( min(nll.h) < fit_DMo$value ){
	optimhDMo <- h_range[which.min(nll.h)]
	min_nll <- min(nll.h)
}
abline(h=min_nll,col="red")
text(optimhDMo,min_nll,signif(optimhDMo,4),pos=3)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- min(fit_DMo$value,min(nll.h)) #We might find a slightly better fit by profiling
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2
abline(h=conf_lim,col="red")

#Calculate CI by interpolation
nll_lo <- nll.h[1:which.min(nll.h)]
h_lo <- h_range[1:which.min(nll.h)]
h_conf_lo <- approx(nll_lo,h_lo,xout=conf_lim)$y
nll_hi <- nll.h[which.min(nll.h):length(nll.h)]
h_hi <- h_range[which.min(nll.h):length(nll.h)]
h_conf_hi <- approx(nll_hi,h_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=h_conf_lo,col="red",lty=2)
abline(v=h_conf_hi,col="red",lty=2)


#The confidence interval
cbind(h_conf_lo,h_conf_hi)


#confidence interval for sd for Density Model

nll.sd <- rep(NA,50)
sd_range <- seq(1.5,2.5,length.out=length(nll.sd))
par <- c(log(optimVDMo),log(optimhDMo))  #starting values for V, h
i <- 1
for (sd in sd_range){
	nll.sd[i] <- optim( par, nll_DMo_norm_sd, sd=sd, trueS=trueS, trueD=trueD, a=a, trueI=trueI )$value
	print(i) #Monitor progress
	i <- i + 1
}

#Plot the profile
plot(sd_range,nll.sd,xlab="sd",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for sd")
axis(1,at=seq(min(sd_range),max(sd_range),by=2),tcl=-0.25,labels=FALSE)


#Add MLE to graph
min_nll <-  fit_DMo$value
#First check if we find a slightly better minimum than optim() from profiling
if ( min(nll.sd) < fit_DMo$value ){
	optimsdDMo <- sd_range[which.min(nll.sd)]
	min_nll <- min(nll.sd)
}
abline(h=min_nll,col="red")
text(optimsdDMo,min_nll,signif(optimsdDMo,4),pos=3)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- min(fit_DMo$value,min(nll.sd)) #We might find a slightly better fit by profiling
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2
abline(h=conf_lim,col="red")

#Calculate CI by interpolation
nll_lo <- nll.sd[1:which.min(nll.sd)]
sd_lo <- sd_range[1:which.min(nll.sd)]
sd_conf_lo <- approx(nll_lo,sd_lo,xout=conf_lim)$y
nll_hi <- nll.sd[which.min(nll.sd):length(nll.sd)]
sd_hi <- sd_range[which.min(nll.sd):length(nll.sd)]
sd_conf_hi <- approx(nll_hi,sd_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=sd_conf_lo,col="red",lty=2)
abline(v=sd_conf_hi,col="red",lty=2)


#The confidence interval
cbind(sd_conf_lo,sd_conf_hi)
nintyfive <- (rbind(cbind(h_conf_lo,h_conf_hi),cbind(sd_conf_lo,sd_conf_hi),cbind(conf_lo,conf_hi))) 
rownames(nintyfive)<- c("h","sd","V")
colnames(nintyfive)<-c("low",'hi')
# diagnostics


residuals <- trueI-DMO_I
#--3 Plots
par(mfrow=c(1,3))

#1. histogram of residuals assessed against theoretical distribution
hist(residuals, xlab = "Residuals", main = "Histogram of residuals",freq=FALSE)
rr <- seq(min(residuals),max(residuals),length.out=100)
lines(rr,dnorm(rr,mean=0,sd=sd(residuals)),col="red")
box()

#2. Residuals vs fitted
plot(DMO_I,residuals, ylab = "Residuals", xlab = "Fitted values",
     main = "Residuals vs fitted")
abline(h=0,col="red")

#3. Quantile-quantile plot
qqnorm(residuals)
qqline(residuals)

#Setup for point deletion

n <- nrow(realdata)
casedel <- matrix(nrow=n,ncol=4)
colnames(casedel) <- list("V","h","sd","nll")
case <- 1:n
V <- optimVDMo
h <- optimhDMo
sd <- optimsdDMo

#Likelihood fit for full dataset
parD <- c(V, log(h), log(sd))  #starting parameter values
fullfit <- optim( parD, nll_DMo_norm, trueD=trueD, trueS=trueS, , a=a, trueI=trueI,method="BFGS")
fullphat <- c(fullfit$par[1], exp(fullfit$par[2]), exp(fullfit$par[3]))
names(fullphat) <- c("V","h","sd")
fullnll <- fullfit$value

#Plot the data and fit


#2. Residuals vs fitted
plot(DMO_I,residuals, ylab = "Residuals", xlab = "Fitted values",
     main = "Residuals vs fitted")
abline(h=0,col="red")


#Standardized residuals are the residuals divided by their
#theoretically expected standard deviation according to the 
#negative log likelihood distribution. Standardized residuals vs
#fitted values plot suggests that smaller fitted values have
#larger residuals.
sr <- residuals / sqrt(DMo_nll.var(V,trueD,trueS,h,a,sd))
plot(DMO_I,sr,xlab="Fitted values",ylab="Standardized residuals", main="Standardized Residuals vs fitted")
abline(h=0,col="red")



#Case deletion fits.
#This is simply refitting the model with each data point left
#out in turn. We record the new parameter estimates and nll.
p <- c(fullfit$par[1], fullfit$par[2], fullfit$par[3])
for (i in case) {
	trueS_minus <- realdata$S[-i]
	trueD_minus <- realdata$D[-i]
	trueI_minus <- realdata$I[-i]
	fit <- optim(p,nll_DMo_norm, trueD_minus=trueD_minus, trueS_minus=trueS_minus, , a=a, trueI_minus=trueI_minus)
	casedel[i,1:3] <- c(fit$par[1], exp(fit$par[2]), exp(fit$par[3]))
	casedel[i,4] <- fit$value
	print(paste("Deleted",i,"of",n,sep=" ")) #Monitoring
}

#Likelihood displacement

LD <- 2 * ( casedel[,"nll"] - fullnll )

#Percent change in parameters

Vpc <- 100 * ( casedel[,"V"] - fullphat["V"] ) / fullphat["V"]
hpc <- 100 * ( casedel[,"h"] - fullphat["h"] ) / fullphat["h"]
sdpc <- 100 * ( casedel[,"sd"] - fullphat["sd"] ) / fullphat["sd"]


max(Vpc) # to figure out what datapoint had the largest effect on nll
max(hpc) # to figure out what datapoint had the largest effect on nll



#Diagnostic plots for point deletion


par(mfrow=c(2,2))

plot(case,abs(LD),xlab="Case",ylab="Likelihood displacement")


plot(case,Vpc,xlab="Case",ylab="Percent change in V")
abline(h=0,col="gray")
text(107,Vpc[107],"107",cex=0.8,pos=2)

plot(case,hpc,xlab="Case",ylab="Percent change in h")
abline(h=0,col="gray")
text(114,hpc[114],"114",cex=0.8, pos=2)

plot(case,sdpc,xlab="Case",ylab="Percent change in sd")
abline(h=0,col="gray")

mtext("Influence - case deletion diagnostics",3,-2.5,outer=TRUE)


