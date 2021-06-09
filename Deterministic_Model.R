#Deterministic Johnson Grass Model
#Programmed by K Harrison Holmes in 2018
#University of Nebraska - Lincoln

# We loop over all crop rotation strategies.
# Populations and allele frequencies are recorded at each time step for each strategy.

source("parameterfunc.r")
parameter.vector <- parameterfunc() #Calls parameters from external file

############
#Parameters#
############
### For genetic model make set.1=1, for same model with allele frequency plot make set.1=2, otherwise make set.1=0 for only demography model ###

set.1 = 1 

### Set simulation time ###

T = 25

if (set.1 == 0) {
  gamma = 0
}
if (set.1 > 0) {
  gamma = parameter.vector[36]			#gene flow
}

### Let h -> -inf if resistance is recessive, h -> inf if resistance is dominant, or h in reals if resistance is partially dominant. ###

h = parameter.vector[1]			#dominance
h.0 = (0/4)^(exp(-h))	#the value and number of the following h's depend on the ploidy
h.1 = (1/4)^(exp(-h))
h.2 = (2/4)^(exp(-h))
h.3 = (3/4)^(exp(-h))
h.4 = (4/4)^(exp(-h))

#MAX4plot = 0 					#Used for scaling the graph


predation = parameter.vector[2]		#Multiply this against s.1J for s.1P
s.1J = parameter.vector[3] 			#seed survival June - Oct
s.1N = parameter.vector[4] 			#seed survival Nov - March
s.1A = parameter.vector[5]			#seed survival April - May

s.2J = parameter.vector[6]			#bud survival June - Oct
s.2N1 = parameter.vector[7]			#Bud survival Oct - March
s.2A = parameter.vector[8] 			#Bud survival April - May
s.2N2 = parameter.vector[9]			#Decreased bud survival Oct - March

s.3J = parameter.vector[10]			#Seedling survival June - flowering

s.4J = parameter.vector[12]			#Sprout survival June - flowering

g.1temp = parameter.vector[14]		#Seed germination April - May
g.2temp = parameter.vector[15]		#Bud sprouting April - May

alpha.3A = parameter.vector[16]
alpha.4A = alpha.3A

delta.3A = parameter.vector[17]
delta.4A = parameter.vector[18]


#Seedling to seed transition
m.SS = parameter.vector[19]					#max production
k.w.SS = parameter.vector[20]				#weed competitiveness
k.c.SS = parameter.vector[21]				#crop competitiveness

#Sprout to seed transition
m.BS = parameter.vector[22]		 			#max production
k.w.BS = parameter.vector[23]				#weed competitiveness
k.c.BS = parameter.vector[24]				#crop competitiveness

#Seedling to bud transition
m.SB = parameter.vector[25]		 			#max production
k.w.SB = parameter.vector[26]				#weed competitiveness
k.c.SB = parameter.vector[27]				#crop competitiveness

#Sprout to bud transition
m.BB = parameter.vector[28]		 			#max production
k.w.BB = parameter.vector[29]				#weed competitiveness
k.c.BB = parameter.vector[30]				#crop competitiveness

d.c = c(parameter.vector[31],parameter.vector[32],parameter.vector[33],parameter.vector[34],parameter.vector[35])	#Density values for Inzen, Sorghum, Soybean, Wheat and Fallow respectively.

sigma = c(parameter.vector[37],parameter.vector[38],parameter.vector[39],parameter.vector[40],parameter.vector[41])	#Survival of herbicide by culm
tau = c(parameter.vector[42],parameter.vector[43],parameter.vector[44],parameter.vector[45],parameter.vector[46])	#Survival of herbicide by bud
mu = parameter.vector[47]		#limited fecundidty for resistant weed

#Create matrices to collect all rotation strategies over time.

AG.matrix = matrix(0,T+1,6)
Seeds.matrix = matrix(0,T+1,6)
Buds.matrix = matrix(0,T+1,6)
Rfreq.matrix = matrix(0,T+1,6)

#Note that 1 = Inzen, 2 = Sorghum, 3 = Soybean, 4 = Wheat, and 5 = Fallow

rotation = list()
rotation[[1]] = c(1)
rotation[[2]] = c(1,2)
rotation[[3]] = c(1,3)
rotation[[4]] = c(1,3,2,3)
rotation[[5]] = c(1,5,4)
rotation[[6]] = c(1,5,4,2,5,4)

for (plot.loop in 1:6) {
  
### Choose crop rotation below ###

set.2 = plot.loop               	#Choose one of the crop rotations below

#^ In this case, the loop collects all rotations for us.
rot.control = rep(rotation[[set.2]],times=(ceiling(T/length(rotation[[set.2]]))+1))
d.c.rot = d.c[rot.control[(1:(T))]]

#we set p.crop and gamma to mark the presence of inzen in rotation
p.crop = c(1)
for (i in 1:(T+1)) {
	if (rot.control[i] == 1) {
		p.crop[i] = 1
	}
	if (rot.control[i] > 1) {
		p.crop[i] = 0
	}
}

sigma.rot = sigma[rot.control[(1:(T+1))]]		
tau.rot = tau[rot.control[(1:(T+1))]]

#Vectors for collecting each stage over time.

S.0 = rep(0,T+1)
S.1 = rep(0,T+1)
S.2 = rep(0,T+1)
S.3 = rep(0,T+1)
S.4 = rep(0,T+1)
B.0 = rep(0,T+1)
B.1 = rep(0,T+1)
B.2 = rep(0,T+1)
B.3 = rep(0,T+1)
B.4 = rep(0,T+1)
P.S.0 = rep(0,T+1)
P.S.1 = rep(0,T+1)
P.S.2 = rep(0,T+1)
P.S.3 = rep(0,T+1)
P.S.4 = rep(0,T+1)
P.B.0 = rep(0,T+1)
P.B.1 = rep(0,T+1)
P.B.2 = rep(0,T+1)
P.B.3 = rep(0,T+1)
P.B.4 = rep(0,T+1)
p.Rcount = rep(0,T+1)
A.G = rep(0,T+1)
All.Seeds = rep(0,T+1)
All.Buds = rep(0,T+1)

#Set initial values

S.0[1] = parameter.vector[48]
All.Seeds[1] = S.0[1]
B.0[1] = parameter.vector[49]
All.Buds[1] = B.0[1]
P.S.0[1] = parameter.vector[50]
P.B.0[1] = parameter.vector[51]
A.G[1] = P.S.0[1] + P.B.0[1]

#Begin time loop

for (t in 1:T) {
	{if (p.crop[t]==1) {
		sigma.t = c(((sigma.rot[t])*(1-h.0)+h.0),((sigma.rot[t])*(1-h.1)+h.1),((sigma.rot[t])*(1-h.2)+h.2),((sigma.rot[t])*(1-h.3)+h.3),((sigma.rot[t])*(1-h.4)+h.4))
		tau.t = c(((tau.rot[t])*(1-h.0)+h.0),((tau.rot[t])*(1-h.1)+h.1),((tau.rot[t])*(1-h.2)+h.2),((tau.rot[t])*(1-h.3)+h.3),((tau.rot[t])*(1-h.4)+h.4))
	}
	else {
		sigma.t = c(sigma.rot[t],sigma.rot[t],sigma.rot[t],sigma.rot[t],sigma.rot[t])
		tau.t = c(tau.rot[t],tau.rot[t],tau.rot[t],tau.rot[t],tau.rot[t])
	}}
	mu.t = c(((1-h.0)+mu*h.0),((1-h.1)+mu*h.1),((1-h.2)+mu*h.2),((1-h.3)+mu*h.3),((1-h.4)+mu*h.4))
	P.surv = (sigma.t[1])*(s.3J*P.S.0[t]+s.4J*P.B.0[t]) + (sigma.t[2])*(s.3J*P.S.1[t]+s.4J*P.B.1[t])+(sigma.t[3])*(s.3J*P.S.2[t]+s.4J*P.B.2[t])+(sigma.t[4])*(s.3J*P.S.3[t]+s.4J*P.B.3[t])+(sigma.t[5])*(s.3J*P.S.4[t]+s.4J*P.B.4[t])
	f.SS = (m.SS*k.w.SS)/(1+k.w.SS*P.surv+k.c.SS*d.c.rot[t])
	f.BS = (m.BS*k.w.BS)/(1+k.w.BS*P.surv+k.c.BS*d.c.rot[t])
	f.SB = (m.SB*k.w.SB)/(1+k.w.SB*P.surv+k.c.SB*d.c.rot[t])
	f.BB = (m.BB*k.w.BB)/(1+k.w.BB*P.surv+k.c.BB*d.c.rot[t])
	P.fec.0 = predation*s.1J*(sigma.t[1])*(mu.t[1])*(s.3J*f.SS*P.S.0[t]+s.4J*f.BS*P.B.0[t])
	P.fec.1 = predation*s.1J*(sigma.t[2])*(mu.t[2])*(s.3J*f.SS*P.S.1[t]+s.4J*f.BS*P.B.1[t])
	P.fec.2 = predation*s.1J*(sigma.t[3])*(mu.t[3])*(s.3J*f.SS*P.S.2[t]+s.4J*f.BS*P.B.2[t])
	P.fec.3 = predation*s.1J*(sigma.t[4])*(mu.t[4])*(s.3J*f.SS*P.S.3[t]+s.4J*f.BS*P.B.3[t])
	P.fec.4 = predation*s.1J*(sigma.t[5])*(mu.t[5])*(s.3J*f.SS*P.S.4[t]+s.4J*f.BS*P.B.4[t])
#	p.tild = (1-gamma)*(((1/4)*P.fec.1+(1/2)*P.fec.2+(3/4)*P.fec.3+P.fec.4)/(P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4))+gamma*p.crop[t]
	if (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4<=0) {p.tild = 0} else 
	{p.tild = (1-gamma)*(((1/4)*P.fec.1+(1/2)*P.fec.2+(3/4)*P.fec.3+P.fec.4)/(P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4))+gamma*p.crop[t]}
	
	culm.spr = (s.1J*s.1N*g.1temp*S.0[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*((1-p.tild)^4)*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.1[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*p.tild*((1-p.tild)^3))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.2[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(6*(p.tild^2)*((1-p.tild)^2))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.3[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*(p.tild^3)*(1-p.tild))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.4[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(p.tild^4)*s.1N*g.1temp) + (s.2J*(tau.t[1])*s.2N2*g.2temp*B.0[t] + s.3J*(sigma.t[1])*f.SB*s.2N1*g.2temp*P.S.0[t] + s.4J*(sigma.t[1])*f.BB*s.2N1*g.2temp*P.B.0[t]) + (s.2J*(tau.t[2])*s.2N2*g.2temp*B.1[t] + s.3J*(sigma.t[2])*f.SB*s.2N1*g.2temp*P.S.1[t] + s.4J*(sigma.t[2])*f.BB*s.2N1*g.2temp*P.B.1[t]) + (s.2J*(tau.t[3])*s.2N2*g.2temp*B.2[t] + s.3J*(sigma.t[3])*f.SB*s.2N1*g.2temp*P.S.2[t] + s.4J*(sigma.t[3])*f.BB*s.2N1*g.2temp*P.B.2[t]) + (s.2J*(tau.t[4])*s.2N2*g.2temp*B.3[t] + s.3J*(sigma.t[4])*f.SB*s.2N1*g.2temp*P.S.3[t] + s.4J*(sigma.t[4])*f.BB*s.2N1*g.2temp*P.B.3[t]) + (s.2J*(tau.t[5])*s.2N2*g.2temp*B.4[t] + s.3J*(sigma.t[5])*f.SB*s.2N1*g.2temp*P.S.4[t] + s.4J*(sigma.t[5])*f.BB*s.2N1*g.2temp*P.B.4[t])
	
	culm.spr.seeds = (s.1J*s.1N*g.1temp*S.0[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*((1-p.tild)^4)*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.1[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*p.tild*((1-p.tild)^3))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.2[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(6*(p.tild^2)*((1-p.tild)^2))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.3[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*(p.tild^3)*(1-p.tild))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.4[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(p.tild^4)*s.1N*g.1temp) 
	
	culm.spr.buds = (s.2J*(tau.t[1])*s.2N2*g.2temp*B.0[t] + s.3J*(sigma.t[1])*f.SB*s.2N1*g.2temp*P.S.0[t] + s.4J*(sigma.t[1])*f.BB*s.2N1*g.2temp*P.B.0[t]) + (s.2J*(tau.t[2])*s.2N2*g.2temp*B.1[t] + s.3J*(sigma.t[2])*f.SB*s.2N1*g.2temp*P.S.1[t] + s.4J*(sigma.t[2])*f.BB*s.2N1*g.2temp*P.B.1[t]) + (s.2J*(tau.t[3])*s.2N2*g.2temp*B.2[t] + s.3J*(sigma.t[3])*f.SB*s.2N1*g.2temp*P.S.2[t] + s.4J*(sigma.t[3])*f.BB*s.2N1*g.2temp*P.B.2[t]) + (s.2J*(tau.t[4])*s.2N2*g.2temp*B.3[t] + s.3J*(sigma.t[4])*f.SB*s.2N1*g.2temp*P.S.3[t] + s.4J*(sigma.t[4])*f.BB*s.2N1*g.2temp*P.B.3[t]) + (s.2J*(tau.t[5])*s.2N2*g.2temp*B.4[t] + s.3J*(sigma.t[5])*f.SB*s.2N1*g.2temp*P.S.4[t] + s.4J*(sigma.t[5])*f.BB*s.2N1*g.2temp*P.B.4[t])
	
	s.3A = 1/(1+alpha.3A*(culm.spr.seeds + delta.3A*culm.spr.buds))	
	s.4A = 1/(1+alpha.4A*(culm.spr.seeds + delta.4A*culm.spr.buds)) 	
	g.1 = s.1N*g.1temp*s.3A
	g.2.1 = s.2N1*g.2temp*s.4A
	g.2.2 = s.2N2*g.2temp*s.4A
	g.1.hat = s.1N*(1-g.1temp)*s.1A
	g.2.1.hat = s.2N1*(1-g.2temp)*s.2A
	g.2.2.hat = s.2N2*(1-g.2temp)*s.2A
	

	S.0[t+1] = s.1J*g.1.hat*S.0[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*((1-p.tild)^4)*g.1.hat
	S.1[t+1] = s.1J*g.1.hat*S.1[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*p.tild*((1-p.tild)^3))*g.1.hat
	S.2[t+1] = s.1J*g.1.hat*S.2[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(6*(p.tild^2)*((1-p.tild)^2))*g.1.hat
	S.3[t+1] = s.1J*g.1.hat*S.3[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*(p.tild^3)*(1-p.tild))*g.1.hat
	S.4[t+1] = s.1J*g.1.hat*S.4[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(p.tild^4)*g.1.hat
	B.0[t+1] = s.2J*(tau.t[1])*g.2.2.hat*B.0[t] + s.3J*(sigma.t[1])*f.SB*g.2.1.hat*P.S.0[t] + s.4J*(sigma.t[1])*f.BB*g.2.1.hat*P.B.0[t]
	B.1[t+1] = s.2J*(tau.t[2])*g.2.2.hat*B.1[t] + s.3J*(sigma.t[2])*f.SB*g.2.1.hat*P.S.1[t] + s.4J*(sigma.t[2])*f.BB*g.2.1.hat*P.B.1[t]
	B.2[t+1] = s.2J*(tau.t[3])*g.2.2.hat*B.2[t] + s.3J*(sigma.t[3])*f.SB*g.2.1.hat*P.S.2[t] + s.4J*(sigma.t[3])*f.BB*g.2.1.hat*P.B.2[t]
	B.3[t+1] = s.2J*(tau.t[4])*g.2.2.hat*B.3[t] + s.3J*(sigma.t[4])*f.SB*g.2.1.hat*P.S.3[t] + s.4J*(sigma.t[4])*f.BB*g.2.1.hat*P.B.3[t]
	B.4[t+1] = s.2J*(tau.t[5])*g.2.2.hat*B.4[t] + s.3J*(sigma.t[5])*f.SB*g.2.1.hat*P.S.4[t] + s.4J*(sigma.t[5])*f.BB*g.2.1.hat*P.B.4[t]
	P.S.0[t+1] = s.1J*g.1*S.0[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*((1-p.tild)^4)*g.1
	P.S.1[t+1] = s.1J*g.1*S.1[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*p.tild*((1-p.tild)^3))*g.1
	P.S.2[t+1] = s.1J*g.1*S.2[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(6*(p.tild^2)*((1-p.tild)^2))*g.1
	P.S.3[t+1] = s.1J*g.1*S.3[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*(p.tild^3)*(1-p.tild))*g.1
	P.S.4[t+1] = s.1J*g.1*S.4[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(p.tild^4)*g.1
	P.B.0[t+1] = s.2J*(tau.t[1])*g.2.2*B.0[t] + s.3J*(sigma.t[1])*f.SB*g.2.1*P.S.0[t] + s.4J*(sigma.t[1])*f.BB*g.2.1*P.B.0[t]
	P.B.1[t+1] = s.2J*(tau.t[2])*g.2.2*B.1[t] + s.3J*(sigma.t[2])*f.SB*g.2.1*P.S.1[t] + s.4J*(sigma.t[2])*f.BB*g.2.1*P.B.1[t]
	P.B.2[t+1] = s.2J*(tau.t[3])*g.2.2*B.2[t] + s.3J*(sigma.t[3])*f.SB*g.2.1*P.S.2[t] + s.4J*(sigma.t[3])*f.BB*g.2.1*P.B.2[t]
	P.B.3[t+1] = s.2J*(tau.t[4])*g.2.2*B.3[t] + s.3J*(sigma.t[4])*f.SB*g.2.1*P.S.3[t] + s.4J*(sigma.t[4])*f.BB*g.2.1*P.B.3[t]
	P.B.4[t+1] = s.2J*(tau.t[5])*g.2.2*B.4[t] + s.3J*(sigma.t[5])*f.SB*g.2.1*P.S.4[t] + s.4J*(sigma.t[5])*f.BB*g.2.1*P.B.4[t]
	p.Rcount[t+1] = ((1/4)*(S.1[t+1]+B.1[t+1]+P.S.1[t+1]+P.B.1[t+1])+(1/2)*(S.2[t+1]+B.2[t+1]+P.S.2[t+1]+P.B.2[t+1])+(3/4)*(S.3[t+1]+B.3[t+1]+P.S.3[t+1]+P.B.3[t+1])+(S.4[t+1]+B.4[t+1]+P.S.4[t+1]+P.B.4[t+1]))/(S.0[t+1]+B.0[t+1]+P.S.0[t+1]+P.B.0[t+1]+S.1[t+1]+B.1[t+1]+P.S.1[t+1]+P.B.1[t+1]+S.2[t+1]+B.2[t+1]+P.S.2[t+1]+P.B.2[t+1]+S.3[t+1]+B.3[t+1]+P.S.3[t+1]+P.B.3[t+1]+S.4[t+1]+B.4[t+1]+P.S.4[t+1]+P.B.4[t+1])
	A.G[t+1] = P.S.0[t+1] + P.S.1[t+1] + P.S.2[t+1] + P.S.3[t+1] + P.S.4[t+1] + P.B.0[t+1] + P.B.1[t+1] + P.B.2[t+1] + P.B.3[t+1] + P.B.4[t+1]
	All.Seeds[t+1] = S.0[t+1] + S.1[t+1] + S.2[t+1] + S.3[t+1] + S.4[t+1]
	All.Buds[t+1] = B.0[t+1] +  B.1[t+1] + B.2[t+1] + B.3[t+1] + B.4[t+1]
}

AG.matrix[,plot.loop] = A.G
Seeds.matrix[,plot.loop] = All.Seeds
Buds.matrix[,plot.loop] = All.Buds
Rfreq.matrix[,plot.loop] = p.Rcount
#MAX4plot = max(MAX4plot,max(A.G,All.Seeds,All.Buds))

}


#Plotting results

library(dplyr)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggplot2)
library(scales)

U=T+1
D.results=data.frame(rot=c(rep(1,U),rep(2,U),rep(3,U),rep(4,U),
                         rep(5,U),rep(6,U)),time=rep(1:U,6),
                   Culms=rep(0,6*U), Seeds=rep(0,6*U),
                   Buds=rep(0,6*U), 
                   P=rep(0,6*U))
for (i in 1:6){
  dum = D.results$rot==i 
  D.results[dum,]$Culms[1:U] = AG.matrix[,i]
  D.results[dum,]$Seeds[1:U] = Seeds.matrix[,i]/10
  D.results[dum,]$Buds[1:U] = Buds.matrix[,i]
  D.results[dum,]$P[1:U] = Rfreq.matrix[,i]
  
}


mytheme=theme(plot.title = element_text(size = 18, face = "bold"),
              axis.text.y = element_text(size=18),
              axis.text.x = element_text(size=18),
              legend.position= c(0.73,0.7), #c(0.75,0.25),
              legend.title =  element_blank(),
              legend.text = element_text(size=14))

mythemeNL=theme(plot.title = element_text(size = 18, face = "bold"),
                axis.text.y = element_text(size=18),
                axis.text.x = element_text(size=18),
                legend.position="none",
                legend.title =  element_blank())
#Figure 4
plotF <-
  ggplot(data=D.results, mapping=aes(x=time, y=P, group = rot, color=factor(rot)))+
  geom_line(aes(y = P), size = 1.1)+
  labs(x = "Time (Years)", y = "Resistance Allele Frequency", size=3) +
  scale_color_brewer(palette="BrBG")+
  theme_bw() +
  theme(axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.title = element_text(size = 18),
        legend.position="none", 
        plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  annotate ("text",label=expression(paste(omega[1]," (IS,...)")),
            y=0.67,x=20, angle=15, size =3)+
  annotate ("text",label=expression(paste(omega[2]," (IS,CS,...)")),
            y=0.511,x=18,angle=20,size=3)+
  annotate ("text",label=expression(paste(omega[3]," (IS,B,...)")),
              y=0.497,x=13,angle=27,size = 3)+
  annotate ("text",label=expression(paste(omega[4]," (IS,B,CS,B,...)")),
              y=0.38,x=22, size =3)+
  annotate ("text",label=expression(paste(omega[5]," (IS,-,W,...)")),
              y=0.48,x=22.5, size = 3)+
  annotate ("text",label=expression(paste(omega[6]," (IS,-,W,CS,-,W,...)")),
              y=0.13,x=21, size =3)

# Saving Figure 4
ggsave("Det_p (Fig.4).pdf",plotF, dpi = 300, device = cairo_pdf)

#Producing subplots for figure 3
plot1 <-
  D.results[,1:5] %>% 
  gather(key = "stage", value = "value", Culms, Seeds, Buds) %>% 
  filter(rot == 1)  %>% 
  ggplot(aes(x=time, y=value, group = stage))+
  geom_line(aes(y = value, color=stage), size = 1.1)+
  scale_color_manual("",values=c("#56B4E9", "#888888", "#332288"),labels=c("Buds", "Culms", "Seeds")) +  
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[1]," (IS,...)")))+
  theme_bw() +
  mytheme +
  scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
  scale_y_continuous(limits=c(0,400), sec.axis=sec_axis(~.*10))
plot1 

plot2 <-
  D.results[,1:5] %>% 
  gather(key = "stage", value = "value", Culms, Seeds, Buds) %>% 
  filter(rot == 2)  %>% 
  ggplot(aes(x=time, y=value, group = stage))+
  geom_line(aes(y = value, color=stage), size = 1.1)+
  scale_color_manual("",values=c("#56B4E9", "#888888", "#332288"),labels=c("Buds", "Culms", "Seeds")) +  
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[2]," (IS,CS, ...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
  scale_y_continuous(limits=c(0,400), sec.axis=sec_axis(~.*10))

plot3 <-
  D.results[,1:5] %>% 
  gather(key = "stage", value = "value", Culms, Seeds, Buds) %>% 
  filter(rot == 3)  %>% 
  ggplot(aes(x=time, y=value, group = stage))+
  geom_line(aes(y = value, color=stage), size = 1.1)+
  scale_color_manual("",values=c("#56B4E9", "#888888", "#332288"),labels=c("Buds", "Culms", "Seeds")) +  
    labs(x = "", y = "") +
  ggtitle(expression(paste(omega[3]," (IS,B, ...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
  scale_y_continuous(limits=c(0,400), sec.axis=sec_axis(~.*10))

plot4 <-
  D.results[,1:5] %>% 
  gather(key = "stage", value = "value", Culms, Seeds, Buds) %>% 
  filter(rot == 4)  %>% 
  ggplot(aes(x=time, y=value, group = stage))+
  geom_line(aes(y = value, color=stage), size = 1.1)+
  scale_color_manual("",values=c("#56B4E9", "#888888", "#332288"),labels=c("Buds", "Culms", "Seeds")) +  
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[4]," (IS,B,CS,B...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,400), sec.axis=sec_axis(~.*10))

plot5 <-
  D.results[,1:5] %>% 
  gather(key = "stage", value = "value", Culms, Seeds, Buds) %>% 
  filter(rot == 5)  %>% 
  ggplot(aes(x=time, y=value, group = stage))+
  geom_line(aes(y = value, color=stage), size = 1.1)+
  scale_color_manual("",values=c("#56B4E9", "#888888", "#332288"),labels=c("Buds", "Culms", "Seeds")) +  
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[5]," (IS,-,W, ...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x), 
                     labels = trans_format("log10", math_format(10^.x)),
                     limits=c(10^-10,1000),
                     sec.axis=sec_axis(trans=~.*10,
                                       breaks = trans_breaks("log10", function(x) 10^x), 
                                       labels = trans_format("log10", math_format(10^.x))))+
  annotation_logticks(sides="l") +
  annotation_logticks(sides="r")

plot6 <-
  D.results[,1:5] %>% 
  gather(key = "stage", value = "value", Culms, Seeds, Buds) %>% 
  filter(rot == 6)  %>% 
  ggplot(aes(x=time, y=value, group = stage))+
  geom_line(aes(y = value, color=stage), size = 1.1)+
  scale_color_manual("",values=c("#56B4E9", "#888888", "#332288"),labels=c("Buds", "Culms", "Seeds")) +  
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[6]," (IS,-,W,CS,-,W,...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x), 
                     labels = trans_format("log10", math_format(10^.x)),
                     limits=c(10^-10,1000),
                     sec.axis=sec_axis(trans=~.*10,
                                       breaks = trans_breaks("log10", function(x) 10^x), 
                                       labels = trans_format("log10", math_format(10^.x))))+
  annotation_logticks(sides="l")+
  annotation_logticks(sides="r")

#Assembling subplots to produce Figure 3

plot.all=
  grid.arrange(grobs=list(plot1,plot2,plot3,plot4,plot5,plot6),
               ncol=2,vp=viewport(width=0.9,height=0.9),
               left=textGrob(expression(paste("Johnsongrass buds and culms per ",m^2)), 
                             gp=gpar(fontsize=25), rot=90),
               right=textGrob(expression(paste("Johnsongrass seeds per ",m^2)), 
                              gp=gpar(fontsize=25), rot=90),
               bottom=textGrob("Time (years)",gp=gpar(fontsize=25))) 

#Saving Figure 3
ggsave("Det_Culms2y.pdf",plot.all,width=10,height = 11, 
       dpi = 300, device = cairo_pdf)



