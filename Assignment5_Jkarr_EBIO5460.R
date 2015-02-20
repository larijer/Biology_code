#Jered Karr
#November 1 2012
#Assignment 5 for EBIO 5460 with Brett Melbourne



#----Function definitions----------------------------------------

#----SIR_ts() ------------------------------
#Returns the time series of the SIR model
#N_init: Initial population size
#alpha Probability that an infected individual infects another
#gamma Probability of recovery (1/mean duration of illness)
#maxt:   Calculated time points will be 0:maxt
#
SIR_ts <- function(N_init,I0,alpha,gamma,maxt) {

#	Initialize population vector
	S <- rep(NA,maxt+1)
	S[1] <- (N_init - I0)
	R <- rep(NA,maxt+1)
	R[1] <- 0
	I <- rep(NA,maxt+1)
	I[1] <- round(I0,0)
	
#	Iterate the model through time
	for (t in 1:maxt){
	S[t+1] <- S[t] - S[t] * ( 1 - (1-alpha)^I[t] )
	I[t+1] <- I[t] + S[t] * ( 1 - (1-alpha)^I[t] ) - gamma*I[t]
	R[t+1] <- R[t] + gamma*I[t]
}
	return(list(S=S,I=I,R=R))
}

#----ssq_SIR_obserr() -----------------------------
#Returns the sum of squares for the SIR model with
#observation error. This is set up for use with
#optim.
#par:    Vector of initial values for N_init,I0,alpha,gamma
#Robs: The data, time series of observed removed
#Iobs: The data, time series of observed infected
#maxt: Passed to SIR_ts().


ssq_SIR_obserr <- function(par,N_init,Iobs,Robs, maxt){
	 
	I0 <- round(par[1],0)      	#Here, I "unpack" the parameter vector 
	alpha <- par[2]		#into its separate parameters. This is not
	gamma <- par[3]     #necessary but the code is now self documenting.
	predict <- SIR_ts(N_init,I0,alpha,gamma,maxt)
	d <- Iobs - predict$I
	b <- Robs - predict$R
	ssq <- sum(d^2) + sum(b^2)
	return(ssq)
}

#######   Read in Data  ###########

islanddata<-read.csv("student14_ass5.csv") #Data from D2L
time <- islanddata$X
Iobs <- islanddata$I
Robs <- islanddata$R
Sobs <- islanddata$S

plot(time,Iobs,xlab="Days",ylab="Number of People (N)",ylim=c(0,6000),main="Infection Thru 40 Days")
points(time,Iobs,pch=8)
lines(time,Iobs)
points(time,Robs,col="red",,pch=8)
lines(time,Robs, col="red")
points(time,Sobs,col=3,,pch=8)
lines(time,Sobs, col=3)
legend("left", c("Susceptible","Infectious","Removed"),pch=rep(16,3),
       col=3:1)
#----Fit the model to the data by using optim()
      
#-------Fit the model to the data by Direct Search
alpha_range <- (seq(0,.1,length.out=100))
I0_range <- seq(1,4,length.out=100) # I don't know if this should only be an interger but maybe...
gamma_range <- seq(0.1,.2,length.out=100)


#Initialize storage objects (matrix, vector, and scalar)
n <- length(I0_range)*length(alpha_range)*length(gamma_range)
results <- matrix(NA,n,4)  #columns: i, s, ssq
colnames(results) <- c("I0","alpha","gamma","ssq")
p_best <- NA
ssq_best <- Inf
N_init<-6000
maxt <- 39
#Direct search over the two parameters: i and s
j <- 1 #Row index for results matrix
for (i in I0_range) {
	for (a in alpha_range) {
		for (g in gamma_range){

par <- c(i,a,g)
ssq <- ssq_SIR_obserr(par,N_init,Iobs,Robs, maxt)
	
	#	Keep the results
		results[j,] <- c(i,a,g,ssq)
		j <- j + 1

	#	Record current best solution
		if (ssq < ssq_best) {
			ssq_best <- ssq
			p_best <- cbind(i,a,g) #cbind automatically gives labels
			}
		}
	}
#	Monitor progress
	print( paste(round(100*j/n),"%",sep=""), quote=FALSE )
}

#Print the minimum sum of squares and best parameter values
ssq_best
p_best

#----Plot sum of squares profiles
quartz(width=10,height=5)
par(mfrow=c(1,3))
scale <- 1.003  #This is a zoom setting. Must be > 1. Smaller is zooming in.
plot(results[,1],log(results[,4]),xlab="i",ylab="log Sum of squares",ylim=c(log(ssq_best),log(scale*ssq_best)), col="blue")
plot(results[,2],log(results[,4]),xlab="a",ylab="log Sum of squares",ylim=c(log(ssq_best),log(scale*ssq_best)), col="blue")
plot(results[,3],log(results[,4]),xlab="g",ylab="log Sum of squares",ylim=c(log(ssq_best),log(scale*ssq_best)), col="blue")
mtext("Sum of squares profiles - direct search",outer=TRUE,line=-2.5)


#Initialize parameters: use values from direct search BUT TRY DIFFERENT STARTING VALUES TOO****

alpha <- p_best[1,2]  #Probability that an infected individual infects another
gamma <- p_best[1,3]  #Probability of recovery (1/mean duration of illness)
N_init <- 6000        #Total population size
I0 <- p_best[1,1]     #Number of initially infected individuals
maxt <- 200    		  #Maximum time to run the simulation start with 39 but can use larger number to make predictions

#Optimization: finding the minimum SSQ
par <- c(I0,alpha,gamma)  #Put the starting parameters in a vector
fit <- optim( par, ssq_SIR_obserr, N_init=N_init, Iobs=Iobs, Robs=Robs, maxt=maxt )
fit   #Note: check convergence code = 0.

#Calculate fitted model dynamics for best parameter values
#Also called "fitted values"
N <- SIR_ts(N_init,I=fit$par[1],alpha=fit$par[2],gamma=fit$par[3],maxt)


#Plot the variables through time
plot(NA,NA,type="n",xlim=c(0,maxt),ylim=c(0,6000),main="Modeled verus Measured",xlab="Days",ylab="N")

points(0:(maxt),N$S,col="green",lwd=2)
lines(0:(maxt),N$S,col="green",lwd=2)
lines(0:(maxt),N$R,col="black",lwd=2)
points(0:(maxt),N$R,col="black",lwd=2)
lines(0:(maxt),N$I,col="red",lwd=2)
points(0:(maxt),N$I,col="red",lwd=2)
legend("left", c("Model Susceptible","Model Infectious","Model Removed"),pch=c(16,16,16),
       col=3:1)
       
#plot predicted       
max(N$R)
max(N$I)
       
plot(NA,NA,type="n",xlim=c(0,maxt),ylim=c(0,6000),main="Infection Model",xlab="Days",ylab="N")

points(0:(maxt),N$S,col="green",lwd=2)
lines(0:(maxt),N$S,col="green",lwd=2)
lines(0:(maxt),N$R,col="black",lwd=2)
points(0:(maxt),N$R,col="black",lwd=2)
lines(0:(maxt),N$I,col="red",lwd=2)
points(0:(maxt),N$I,col="red",lwd=2)
abline(h=max(N$R),col="orange",lwd=3)
abline(h=max(N$I),col="blue",lwd=3)
legend("right", c("Susceptible","Infectious","Removed","Max Dead","Max infected"),lwd=rep(2,5)
       ,col=c(3:1,"orange","blue"))

