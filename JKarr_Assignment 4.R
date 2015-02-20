# Jered Karr
# Friday October 19th 2012
# Assignment 4 for EBIO 5460



##############Pseudocode#######################

# This is a model for the Nitrogen cycle with 4 componets (plants, herbivores, Biounavailable Soil Nitrogen, Bioavailable Soil Nitrogen)
# controlled fires destroy 30% of the plant matter when Plant Biomass rises above 25 kg/ha
# There is a steady stream of Bio available Nitrogen from the Atmosphere adding 10 kg/ha/year.
# Total Nitrogen for system is 80 Kg/ha.
# Run for 1000 days
# Soil Bio Unavailable Nitrogen is slowly converted to Soil Bio-available Nitrogen at a rate of every month or (1/30 per day) Herbivores add 5% of N mass everyday and die very iregularly plants also add to this pool by dieing every 3 months (1/90 pre day of Plant Mass is converted to this pool)
# Soil Bio-available Nitrogen pool is being added to by Soil Bio-unavailable Nitrogen getting converted (add 1/30 0f Soil Bio-unavailable Nitrogen per day) Soil Bio-available Nitrogen is being used by plants at a rate of 0.002 kg N per Ha per day when nitrogen concentrion is 1 kg/ha N. Atmospheric deposition of N is adding 10/365 Kg/ha/day.
#Herbivore growth is dependenant on Plant Biomass(Relationship of .0028 Kg N). They lose N at a steady rate of 5% per day and a very small amount due to death.
#  Plant growth rate is dependant on Soil Bio-available Nitrogen. But is limited by Herbivore growth rate and steady death(1/90 of total plant Biomass is converted to Soil Biounavailable Nitrogen per day).
# Last due to acumalation of plant material controlled fires happen whenver the total Nitrogen of plants exceed 25 Kg/ha reducing it by 30%
# These values are converted to a vector and kept later to graph and examine patterns observed.


######CODE############


#Initialize variables
SUAN <- 20  # Soil Bio Unavailable Nitrogen
SAN <- 44	# Soil Bio-Available Nitrogen
PN <- 15	# Plant Biomass
HN <- 1	# Herbivore Biomass change to 0 to exclude herbivores from the sytem
t <- 0		# days
End <- 1000 # days

#Initializing values
excrete <- .05 # Percentage of Nitrogen Herbivores excrete every day
plt_grth <- .002 # rate of plant growth
herb_grth <- .0028 # rate of animal growth
death <- 1-.9^(1/365) #death rate of herbivores
at_dep <- (10/365) # rate of atmospheric deposition
fire <- .7	# Amount of Nitrogen is left after fires
#Initialize vectors to store the results
SUANv <- SUAN
SANv <- SAN
PNv <- PN
HNv <- HN
tv <- t


for(i in 0:End){	
	if(PN > 25){ #gatekeeper to decide if fire is needed
		SAND <- SAN # Making copies of original values so they stay the same throughout time steps 
		SUAND <- SUAN
		HND <- HN
		PND <- PN
		
		# Soil Avaiable Nitrogen (-Plant growth, + Soil N conversion and atmosphere Depostion)
		SAN <- SAND - (SAND * PND*plt_grth) + (SUAND/30) + at_dep #control at_dep by adding # in front of + at_dep
		
		#Soil UnAvailable Nitrogen (-Soil N conversion +Plant Death, Herbivore Death, Herbivore excretement)
		SUAN <- SUAND + excrete * HND - (SUAND/30) + (PND/90) + (HND * death)
		
		#Herbivore Nitrogen (-Herbivore excretement, Herbivore death +Herbivore growth due to plant consumption)
		HN <- HND - (excrete * HND) + (PND * HND * herb_grth) -(HND * death)
		
		#Plant Nitrogen(-Herbivore growth due to plant consumption,Plant Death + Plant growth due to uptake of SAN) Also fire destroys 30% of plant Biomass
		PN <- (PND + (SAND * PND * plt_grth) - (PND * HND * herb_grth) - (PND/90)) * fire #control adding fires by putting # in front of *fire
		
		t <- t+1 #time counter in case I want to use while
		SUANv <- c(SUANv, SUAN) # Vectors to store values
		SANv <- c(SANv, SAN)
		PNv <- c(PNv, PN)
		tv <- c(tv, t)
		HNv <- c(HNv, HN)
		
	}else { #if plant Biomass is less than 25 no fire!
		SAND <-SAN
		SUAND <- SUAN
		HND <- HN
		PND <- PN
		SAN <- SAND - (SAND * PND*plt_grth) + (SUAND/30) + at_dep #control at_dep by adding # in front of + at_dep
		SUAN <- SUAND + excrete * HND - (SUAND/30) + (PND/90) + (HND * death)
		HN <- HND - (excrete * HND) + (PND * HND * herb_grth) -(HND * death)
		PN <- PND + (SAND * PND * plt_grth) - (PND * HND * herb_grth) - (PND/90)
		t <- t+1
		SUANv <- c(SUANv, SUAN)
		SANv <- c(SANv, SAN)
		PNv <- c(PNv, PN)
		tv <- c(tv, t)
		HNv <- c(HNv, HN)	

	}
}

plot(NA,NA,type="n",xlim=c(0,tv[length(tv)]),ylim=c(0, 50),xlab="Time(days) ",ylab="N (kg/ha)", main="N cycle")
lines(tv, HNv, col="red",lwd=2)
lines(tv, SUANv, col="blue", lwd=2)
lines(tv, SANv, col="green", lwd=2)
lines(tv, PNv, lwd=2)
legend("topright", c("Soil Unavailable","Soil Available","Herbivore","Plant" ),pch=rep(16,3), col=4:1)

HN + PN + SUAN + SAN #check amount of Nitrogen