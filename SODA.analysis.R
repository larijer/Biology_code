#SODA data analysis
#written by: Broken wing, Feb 2014



SODA <-read.csv("~/Documents/HoloClimate Research/Data Analysis/SODA/SODA.temperature.profiles.soledad.csv")
depths <-read.csv("~/Documents/HoloClimate Research/Data Analysis/SODA/depths.csv",header=F)


#SODA <- read.csv("SODA.temperature.profiles.soledad.csv")
#depths <-read.csv("depths.csv",header=F)

#how many years do I have?
length(levels(as.factor(SODA$year)))
levels(as.factor(SODA$year))[10]
mon.col <-rainbow(12)
Mon.col <-rep(NA,12)
Mon.col[6] <- mon.col[1]
Mon.col[7] <- mon.col[2]
Mon.col[5] <- mon.col[3]
Mon.col[8] <- mon.col[4]
Mon.col[4] <- mon.col[5]
Mon.col[9] <- mon.col[6]
Mon.col[3] <- mon.col[7]
Mon.col[10] <- mon.col[8]
Mon.col[2] <- mon.col[9]
Mon.col[11] <- mon.col[10]
Mon.col[1] <- mon.col[11]
Mon.col[12] <- mon.col[12]

length(levels(as.factor(SODA$month)))

#y= months
#depths = duh


depth.profiles <- function(depths,y,x,P,cul,year){
	
	plot(NA, NA, ylim=(c(max(as.numeric(depths)),0)),xlim= c(min(x[, 3:22]),max(x[, 3:22])),xlab="Temperature",ylab="Depth",main=levels(as.factor(year)))

for(i in 1:length(levels(as.factor(y)))){
	
	lines(as.numeric(P[i,3:22]) ,as.numeric(depths), type ="l",col=cul[i])
	text(x[i,3] , y=0, levels(as.factor(y[i])), col=cul[i])
	
	
}

legend("bottomright",levels(as.factor(y)),col=cul,lty =1)
}


depth.profiles.one <- function(depths,y,x,P,cul,year){
	
	plot(NA, NA, ylim=(c(max(as.numeric(depths)),0)),xlim= c(min(x[, 3:22]),max(x[, 3:22])),xlab="Temperature",ylab="Depth",main=levels(as.factor(year)))

for(i in 1:length(levels(as.factor(y)))){
	
	lines(as.numeric(P[i,3:22]) ,as.numeric(depths), type ="l",col=cul[i])
	text(x[i,3] , y=0, levels(as.factor(y[i])), col=cul[i])
	
	
}
	lines(as.numeric(colMeans(P[,3:22])) ,as.numeric(depths), type ="l",col="black",lwd=2)
legend("bottomright",levels(as.factor(y)),col=cul,lty =1)
}





list.function <-function(y,x){
	name.list <-levels(as.factor(x))
	lis <- vector("list",length(levels(as.factor(x))))
	tit <- rep(NA,length(name.list))
	
	for(i in 1:length(levels(as.factor(x)))){
		lis[[i]] <-assign(name.list[i],subset(y,x==name.list[i]))
		
	}
	return(lis)
}



SO <- list.function(SODA,SODA$year)


a <- c(128, 118, 110)
quartz()
par(mfrow= c(1,3))

for (i in a){
	
	depth.profiles.one(depths,SO[[c(i,2)]],SO[[i]],SO[[i]],Mon.col,SO[[c(i,1)]])
	
	
}

par(mfrow= c(1,3))
xx <-rep(NA,138)
for (i in 1:138){
	xx[i] <-SO[[i]][3,3]-SO[[i]][3,8]


}
xx <-na.omit(xx)
quartz()
plot(levels(as.factor(SODA[,1])),as.numeric(xx),type="b",ylim=c(min(xx),max(xx)))
lines(levels(as.factor(SODA[,1])),as.numeric(f.xx),col="blue")
lines(levels(as.factor(SODA[,1])),as.numeric(xx),col="red")
lines(levels(as.factor(SODA[,1])),as.numeric(xx),col="green")

ma5 = c(1,1,1,1,1)/5
f.xx <-filter(xx,ma5)

ha <-numeric(length(xx))
ha[f.xx-xx>1.5] <-"el nino"

f.xx[128]-xx[128]

1998-128

plot(xx)


