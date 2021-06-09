# Stochastic Johnson Grass Model 
#Programmed by K Harrison Holmes in 2018
#University of Nebraska - Lincoln
#Modified by Brigitte Tenhumberg 2021

### Changes relative to deterministic model ###
# The model appears twice in this script. The first time it is looped to record the stochastic trials, the second time to run a deterministic simulation.
# The parameter random.tests determines the number of stochastic trials.
# At the end of the stochastic trials, population numbers and allele frequencies are sorted in ascending order for each time step in order to deterrmine percentiles for the plot.
# This script is for all stochastic plots. Scroll to bottom of script to see instructions on how to edit plot code by commenting.
# Variance is defined within the two stochastic functions in the loop below. Standard variance is set to 0.0025. This appears four times in beta functions.


source("parameterfunc.r")
parameter.vector <- parameterfunc() #Calls parameters from external file

#estimating the parameters of the beta distribution given a mean mu and variance var
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

############
#Parameters#
############

### Let h -> -inf if resistance is recessive, h -> inf if resistance is dominant, or h in reals if resistance is partially dominant. ###

h = parameter.vector[1]			#dominance
h.0 = (0/4)^(exp(-h))
h.1 = (1/4)^(exp(-h))
h.2 = (2/4)^(exp(-h))
h.3 = (3/4)^(exp(-h))
h.4 = (4/4)^(exp(-h))

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

alpha.3A = parameter.vector[16]    #parameterrs for self thinning
alpha.4A = alpha.3A

delta.3A = parameter.vector[17]    #parameterrs for self thinning
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

### For genetic model make set.1=1, for same model with allele frequency plot make set.1=2, otherwise make set.1=0 for only demography model ###

set.1 = 2

### Set max time ###
T = 25
random.tests = 10000
AG.rot1 = matrix(0,T+1,6*random.tests)
p.rot1 = matrix(0,T+1,6*random.tests)

#crop rotations
#Note that 1 = Inzen, 2 = Sorghum, 3 = Soybean, 4 = Wheat, and 5 = Fallow

rotation = list()
rotation[[1]] = c(1)
rotation[[2]] = c(1,2)
rotation[[3]] = c(1,3)
rotation[[4]] = c(1,3,2,3)
rotation[[5]] = c(1,5,4)
rotation[[6]] = c(1,5,4,2,5,4)
#crop density at different rotations
d.c = c(parameter.vector[31],parameter.vector[32],parameter.vector[33],parameter.vector[34],parameter.vector[35])	#Density values for Inzen, Sorghum, Soybean, Wheat and Fallow respectively.
#culm and bud mortality at different rotations
sigma = c(parameter.vector[37],parameter.vector[38],parameter.vector[39],parameter.vector[40],parameter.vector[41])	#Survival of herbicide by culm
tau = c(parameter.vector[42],parameter.vector[43],parameter.vector[44],parameter.vector[45],parameter.vector[46])	#Survival of herbicide by bud
mu = parameter.vector[47]		#limited fecundidty for resistant weed



for (rotation.loop in 1:6) {     #rotation loop
  ### Choose crop rotation below ###
  
  set.2 = rotation.loop               #Choose one of the crop rotations below
  
  rot.control = rep(rotation[[set.2]],times=(ceiling(T/length(rotation[[set.2]]))+1))
 
  #which crop rotation includes gene flow
  p.crop = c(1)
  for (i in 1:(T+1)) {
    if (rot.control[i] == 1) {
      p.crop[i] = 1
    }
    if (rot.control[i] > 1) {
      p.crop[i] = 0
    }
  }
  
  if (set.1 == 0) {
    gamma = 0
  }
  if (set.1 > 0) {
    gamma = parameter.vector[36]			#gene flow
  }
  
  d.c.rot = d.c[rot.control[(1:(T))]]
  sigma.rot = sigma[rot.control[(1:(T+1))]]
  tau.rot = tau[rot.control[(1:(T+1))]]
  sigma.shape1 = estBetaParams(sigma.rot,0.0025)$alpha
  sigma.shape2 = estBetaParams(sigma.rot,0.0025)$beta
  tau.shape1 = estBetaParams(tau.rot,0.0025)$alpha
  tau.shape2 = estBetaParams(tau.rot,0.0025)$beta
  
  #estimating the parameters of the beta distribution given a mean mu and variance var
  estBetaParams(0.01, 0.0025)
for (r.t in 1:random.tests) {    #replications  
#creating population vector
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

#starting condition
S.0[1] = parameter.vector[48]
S.1[1] = parameter.vector[52]
S.2[1] = parameter.vector[53]
S.3[1] = parameter.vector[54]
S.4[1] = parameter.vector[55]
B.0[1] = parameter.vector[49]
P.S.0[1] = parameter.vector[50]
P.B.0[1] = parameter.vector[51]
A.G[1] = P.S.0[1] + P.B.0[1]
AG.rot1[,(rotation.loop-1)*random.tests+r.t] = A.G[1]
p.rot1[,(rotation.loop-1)*random.tests+r.t] = p.Rcount

for (t in 1:T) {
	sigma.rot[t] = rbeta(1,sigma.shape1[t],sigma.shape2[t])
	tau.rot[t] = rbeta(1,tau.shape1[t],tau.shape2[2])
	#adjusting mortality based on number of resistance genes if Inzen is planted
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
	if (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4<=0) {p.tild = 0} else 
	{p.tild = (1-gamma)*(((1/4)*P.fec.1+(1/2)*P.fec.2+(3/4)*P.fec.3+P.fec.4)/(P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4))+gamma*p.crop[t]}
		
#	if (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4<=0) stop("P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4<=0")
	
	culm.spr = (s.1J*s.1N*g.1temp*S.0[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*((1-p.tild)^4)*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.1[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*p.tild*((1-p.tild)^3))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.2[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(6*(p.tild^2)*((1-p.tild)^2))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.3[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*(p.tild^3)*(1-p.tild))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.4[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(p.tild^4)*s.1N*g.1temp) + (s.2J*(tau.t[1])*s.2N2*g.2temp*B.0[t] + s.3J*(sigma.t[1])*f.SB*s.2N1*g.2temp*P.S.0[t] + s.4J*(sigma.t[1])*f.BB*s.2N1*g.2temp*P.B.0[t]) + (s.2J*(tau.t[2])*s.2N2*g.2temp*B.1[t] + s.3J*(sigma.t[2])*f.SB*s.2N1*g.2temp*P.S.1[t] + s.4J*(sigma.t[2])*f.BB*s.2N1*g.2temp*P.B.1[t]) + (s.2J*(tau.t[3])*s.2N2*g.2temp*B.2[t] + s.3J*(sigma.t[3])*f.SB*s.2N1*g.2temp*P.S.2[t] + s.4J*(sigma.t[3])*f.BB*s.2N1*g.2temp*P.B.2[t]) + (s.2J*(tau.t[4])*s.2N2*g.2temp*B.3[t] + s.3J*(sigma.t[4])*f.SB*s.2N1*g.2temp*P.S.3[t] + s.4J*(sigma.t[4])*f.BB*s.2N1*g.2temp*P.B.3[t]) + (s.2J*(tau.t[5])*s.2N2*g.2temp*B.4[t] + s.3J*(sigma.t[5])*f.SB*s.2N1*g.2temp*P.S.4[t] + s.4J*(sigma.t[5])*f.BB*s.2N1*g.2temp*P.B.4[t])
	
	culm.spr.seeds = (s.1J*s.1N*g.1temp*S.0[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*((1-p.tild)^4)*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.1[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*p.tild*((1-p.tild)^3))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.2[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(6*(p.tild^2)*((1-p.tild)^2))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.3[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(4*(p.tild^3)*(1-p.tild))*s.1N*g.1temp) + (s.1J*s.1N*g.1temp*S.4[t] + (P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4)*(p.tild^4)*s.1N*g.1temp)
	
	culm.spr.buds = (s.2J*(tau.t[1])*s.2N2*g.2temp*B.0[t] + s.3J*(sigma.t[1])*f.SB*s.2N1*g.2temp*P.S.0[t] + s.4J*(sigma.t[1])*f.BB*s.2N1*g.2temp*P.B.0[t]) + (s.2J*(tau.t[2])*s.2N2*g.2temp*B.1[t] + s.3J*(sigma.t[2])*f.SB*s.2N1*g.2temp*P.S.1[t] + s.4J*(sigma.t[2])*f.BB*s.2N1*g.2temp*P.B.1[t]) + (s.2J*(tau.t[3])*s.2N2*g.2temp*B.2[t] + s.3J*(sigma.t[3])*f.SB*s.2N1*g.2temp*P.S.2[t] + s.4J*(sigma.t[3])*f.BB*s.2N1*g.2temp*P.B.2[t]) + (s.2J*(tau.t[4])*s.2N2*g.2temp*B.3[t] + s.3J*(sigma.t[4])*f.SB*s.2N1*g.2temp*P.S.3[t] + s.4J*(sigma.t[4])*f.BB*s.2N1*g.2temp*P.B.3[t]) + (s.2J*(tau.t[5])*s.2N2*g.2temp*B.4[t] + s.3J*(sigma.t[5])*f.SB*s.2N1*g.2temp*P.S.4[t] + s.4J*(sigma.t[5])*f.BB*s.2N1*g.2temp*P.B.4[t])
	
	s.3A = ifelse(1+alpha.3A*(culm.spr.seeds + delta.3A*culm.spr.buds)==0,1, 1/(1+alpha.3A*(culm.spr.seeds + delta.3A*culm.spr.buds))) 	
	s.4A = ifelse(1+alpha.4A*(culm.spr.seeds + delta.4A*culm.spr.buds)==0,1, 1/(1+alpha.4A*(culm.spr.seeds + delta.4A*culm.spr.buds))) 	
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
}

AG.rot1[,(rotation.loop-1)*random.tests+r.t] = A.G			#These matrices contain all times steps for a trial in a particular column. Every random.tests we switch to a new rotation strategy.
p.rot1[,(rotation.loop-1)*random.tests+r.t] = p.Rcount

}}

AG.rot6=AG.rot1[,((6-1)*random.tests+1):(6*random.tests)]
AG.rot5=AG.rot1[,((5-1)*random.tests+1):(5*random.tests)]
AG.rot4=AG.rot1[,((4-1)*random.tests+1):(4*random.tests)]
AG.rot3=AG.rot1[,((3-1)*random.tests+1):(3*random.tests)]
AG.rot2=AG.rot1[,((2-1)*random.tests+1):(2*random.tests)]
AG.rot1=AG.rot1[,1:random.tests]
p.rot6=p.rot1[,((6-1)*random.tests+1):(6*random.tests)]
p.rot5=p.rot1[,((5-1)*random.tests+1):(5*random.tests)]
p.rot4=p.rot1[,((4-1)*random.tests+1):(4*random.tests)]
p.rot3=p.rot1[,((3-1)*random.tests+1):(3*random.tests)]
p.rot2=p.rot1[,((2-1)*random.tests+1):(2*random.tests)]
p.rot1=p.rot1[,1:random.tests]




#Running deterministic version
# model parameters stay teh same as in teh stochastic model

#create result matrix
AG.matrix = matrix(0,T+1,6)
Seeds.matrix = matrix(0,T+1,6)
Buds.matrix = matrix(0,T+1,6)
Rfreq.matrix = matrix(0,T+1,6)
p.Rcount.matrix = matrix(0,T+1,6)

for (plot.loop in 1:6) {

### Choose crop rotation below ###

set.2 = plot.loop               	#Choose one of the crop rotations below
rot.control = rep(rotation[[set.2]],times=(ceiling(T/length(rotation[[set.2]]))+1))
d.c.rot = d.c[rot.control[(1:(T))]]

p.crop = c(1)
for (i in 1:(T+1)) {
	if (rot.control[i] == 1) {
		p.crop[i] = 1
	}
	if (rot.control[i] > 1) {
		p.crop[i] = 0
	}
}

if (set.1 == 0) {
	gamma = 0
}
if (set.1 > 0) {
	gamma = parameter.vector[36]			#gene flow
}

sigma.rot = sigma[rot.control[(1:(T+1))]]
tau.rot = tau[rot.control[(1:(T+1))]]

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

S.0[1] = parameter.vector[48]
S.1[1] = parameter.vector[52]
S.2[1] = parameter.vector[53]
S.3[1] = parameter.vector[54]
S.4[1] = parameter.vector[55]
All.Seeds[1] = S.0[1]
B.0[1] = parameter.vector[49]
All.Buds[1] = B.0[1]
P.S.0[1] = parameter.vector[50]
P.B.0[1] = parameter.vector[51]
A.G[1] = P.S.0[1] + P.B.0[1]

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
	p.tild = (1-gamma)*(((1/4)*P.fec.1+(1/2)*P.fec.2+(3/4)*P.fec.3+P.fec.4)/(P.fec.0+P.fec.1+P.fec.2+P.fec.3+P.fec.4))+gamma*p.crop[t]
	
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

p.Rcount.matrix[,plot.loop] = p.Rcount

}

library(dplyr)
library(tidyverse)
U=T+1
results=data.frame(rot=c(rep(1,U),rep(2,U),rep(3,U),rep(4,U),
                         rep(5,U),rep(6,U)),time=rep(1:U,6),
                   median=rep(0,6*U), LCI=rep(0,6*U),
                   UCI=rep(0,6*U), 
                   det=rep(0,6*U))
dum = results$rot==1 
for (t in 1:U){
  results$median[t] = median(p.rot1[t,], na.rm=TRUE)
  results$LCI[t] = quantile(p.rot1[t,], probs = 0.025)
  results$UCI[t] = quantile(p.rot1[t,], probs = 0.975)
}
results[dum,]$det[1:U] = p.Rcount.matrix[,1]

dum = results$rot==2 
for (t in 1:U){
  results[dum,]$median[t] = median(p.rot2[t,], na.rm=TRUE)
  results[dum,]$LCI[t] = quantile(p.rot2[t,], probs = 0.025)
  results[dum,]$UCI[t] = quantile(p.rot2[t,], probs = 0.975)
}
results[dum,]$det[1:U] = p.Rcount.matrix[,2]

dum = results$rot==3 
for (t in 1:U){
  results[dum,]$median[t] = median(p.rot3[t,], na.rm=TRUE)
  results[dum,]$LCI[t] = quantile(p.rot3[t,], probs = 0.025)
  results[dum,]$UCI[t] = quantile(p.rot3[t,], probs = 0.975)
}
results[dum,]$det[1:U] = p.Rcount.matrix[,3]

dum = results$rot==4 
for (t in 1:U){
  results[dum,]$median[t] = median(p.rot4[t,], na.rm=TRUE)
  results[dum,]$LCI[t] = quantile(p.rot4[t,], probs = 0.025)
  results[dum,]$UCI[t] = quantile(p.rot4[t,], probs = 0.975)
}
results[dum,]$det[1:U] = p.Rcount.matrix[,4]

dum = results$rot==5 
for (t in 1:U){
  results[dum,]$median[t] = median(p.rot5[t,], na.rm=TRUE)
  results[dum,]$LCI[t] = quantile(p.rot5[t,], probs = 0.025, na.rm=TRUE)
  results[dum,]$UCI[t] = quantile(p.rot5[t,], probs = 0.975, na.rm=TRUE)
}
results[dum,]$det[1:U] = p.Rcount.matrix[,5]

dum = results$rot==6 
for (t in 1:U){
  results[dum,]$median[t] = median(p.rot6[t,], na.rm=TRUE)
  results[dum,]$LCI[t] = quantile(p.rot6[t,], probs = 0.025,na.rm=TRUE)
  results[dum,]$UCI[t] = quantile(p.rot6[t,], probs = 0.975,na.rm=TRUE)
}
results[dum,]$det[1:U] = p.Rcount.matrix[,6]

library(tidyverse)
library(grid)
#library (dplyr)
#library(cowplot)
library(gridExtra)
mytheme=theme(plot.title = element_text(size = 18, face = "bold"),
              axis.text.y = element_text(size=18),
              axis.text.x = element_text(size=18),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              legend.position=c(0.7,0.2),
              legend.title =  element_blank(),
              legend.text = element_text(size=12))

mythemeNL=theme(plot.title = element_text(size = 18, face = "bold"),
                axis.text.y = element_text(size=18),
                axis.text.x = element_text(size=18),
                panel.border = element_rect(color = "black", fill = NA, size = 1),
                legend.position="none",
                legend.title =  element_blank())


plot1 <- results %>% 
  gather(key = "statistic", value = "value", median, det) %>% 
  filter(rot == 1) %>%
  ggplot(aes(x=time, y=value, group = statistic))+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(statistic)), alpha = 0.2,show.legend = FALSE) +
  scale_color_manual(values=c("#D55E00","#44AA99"),labels = c("No Variation", "Variation"))+
  geom_line(aes(y = value, color=statistic),size = 1.3)+
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[1]," (IS,...)")))+
  theme_bw() +
  mytheme +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))
#scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0))
plot1

#scale_y_log10()
plot2 <- results %>% 
  gather(key = "statistic", value = "value", median, det) %>% 
  filter(rot == 2) %>%
  ggplot(aes(x=time, y=value, group = statistic))+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(statistic)), alpha = 0.2,show.legend = FALSE) +
  scale_color_manual(values=c("#D55E00","#44AA99"))+
  geom_line(aes(y = value, color=statistic),size = 1.3)+
  #scale_color_grey()+
  labs(x = "", y = "") +
  ggtitle(expression(paste(omega[2]," (IS,CS,...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))
plot2

plot3 <- results %>% 
  gather(key = "statistic", value = "value", median, det) %>% 
  filter(rot == 3) %>%
  ggplot(aes(x=time, y=value, group = statistic))+
  #geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "grey70")+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(statistic)), alpha = 0.2,show.legend = FALSE) +
  scale_color_manual(values=c("#D55E00","#44AA99"))+
  geom_line(aes(y = value, color=statistic),size = 1.3)+
  #scale_color_grey()+
  labs(x = "", y = "") + 
  ggtitle(expression(paste(omega[3]," (IS,B,...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))

plot4 <- results %>% 
  gather(key = "statistic", value = "value", median, det) %>% 
  filter(rot == 4) %>%
  ggplot(aes(x=time, y=value, group = statistic))+
  #geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "grey70")+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(statistic)), alpha = 0.2,show.legend = FALSE) +
  scale_color_manual(values=c("#D55E00","#44AA99"))+
  geom_line(aes(y = value, color=statistic),size = 1.3)+
  #scale_color_grey()+
  labs(x = "", y = "")+ 
  ggtitle (expression(paste(omega[4],"  (IS,B,CS,B...)")))+
  theme_bw() +
  mythemeNL+
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))

plot5 <- results %>% 
  gather(key = "statistic", value = "value", median, det) %>% 
  filter(rot == 5) %>%
  ggplot(aes(x=time, y=value, group = statistic))+
  #geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "grey70")+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(statistic)), alpha = 0.2,show.legend = FALSE) +
  scale_color_manual(values=c("#D55E00","#44AA99"))+
  geom_line(aes(y = value, color=statistic),size = 1.3)+
  #scale_color_grey()+
  labs(x = "", y = "") + 
  ggtitle(expression(paste(omega[5],"  (IS,-,W,...)")))+
  theme_bw() +
  mythemeNL+
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))

plot6 <- results %>% 
  gather(key = "statistic", value = "value", median, det) %>% 
  filter(rot == 6) %>%
  ggplot(aes(x=time, y=value, group = statistic))+
  #geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "grey70")+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(statistic)), alpha = 0.2,show.legend = FALSE) +
  scale_color_manual(values=c("#D55E00","#44AA99"))+
  geom_line(aes(y = value, color=statistic),size = 1.3)+
 # scale_color_grey()+
  labs(x = "", y = "")+
  ggtitle(expression(paste(omega[6]," (IS,-,W,CS,-,W,...)")))+
  theme_bw() +
  mythemeNL +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))

plot.all=
  grid.arrange(grobs=list(plot1,plot2,plot3,plot4,plot5,plot6),
               ncol=2,vp=viewport(width=0.9,height=0.9),
               left=textGrob("Resistance Allele Frequency", 
                           gp=gpar(fontsize=25), rot=90),
               bottom=textGrob("Time (years)",gp=gpar(fontsize=25))) 

ggsave("Stochastic_ptest0.0006.pdf",plot.all,width=8.5,height = 11,
       dpi = 300, device = cairo_pdf)

ggsave("Stochastic_ptest.pdf",plot.all,width=8.5,height = 11,
       dpi = 300, device = cairo_pdf)
