parameterfunc <- function(x) {
	
	h <- 7					#dominance
	
	predation <- 0.5		#Multiply this against s.1J for s.1P
	
	s.1J <- 0.7 				#seed survival June - Oct
	s.1N <- 0.25 			#seed survival Nov - March
	s.1A <- 0.97 			#seed survival April - May
	
	s.2J <- 0.5				#bud survival June - Oct
	s.2N1 <- 0.5				#Bud survival Oct - March
	s.2A <- 0.98				#Bud survival April - May
	s.2N2 <- 0.1				#Decreased bud survival Oct - March
	
	s.3J = 0.9  				#Seedling survival June - flowering
		
	s.4J = 0.95				#Sprout survival June - flowering
	
	g.1temp = 0.35			#Seed germination April - May
	g.2temp = 0.73			#Bud sprouting April - May
	
	alpha.3A = 0.008 	#1/maximum density per area	
	
	delta.3A = 2      #seedling competitiveness
	
	delta.4A = 0.5     #sprout competitiveness	

	#makes it easier to change parameters 
  #old parameters
	#k.w = 1            
#	k.c = 0.8392    
#	f = 0.125     
	#new parameters
	k.w = 0.00549   
  k.c = 0.029
	f = 1
	
	#Seedling to seed transition
	m.SS = 200000 *f				#max production
	k.w.SS = k.w					#weed competitiveness
	k.c.SS = k.c				#crop competitiveness
	
	#Sprout to seed transition
	m.BS = 280000*f		 			#max production
	k.w.BS = k.w					#weed competitiveness
	k.c.BS = k.c				#crop competitiveness
	
	#Seedling to bud transition
	m.SB = 4000	*f	 			#max production
	k.w.SB = k.w					#weed competitiveness
	k.c.SB = k.c				#crop competitiveness
	
	#Sprout to bud transition
	m.BB = 4800	*f	 			#max production
	k.w.BB = k.w					#weed competitiveness
	k.c.BB = k.c				#crop competitiveness

	d.c=c(10,10,33,222,100)	#Density values for Inzen, Sorghum, Soybean, Wheat and Fallow respectively.
	
	gamma = 0.0001 				#gene flow
	sigma =  c(0.01,0.4,0.01,0.005,0.005)			#Survival of herbicide by culm (same order as above)
	tau =  c(0.04,0.4,0.05,0.02,0.01)					#Survival of herbicide by bud (same order as above)
	mu = 1						#limited fecundidty for resistant weed
	
	#S.1.init=S.2.init=S.3.init=S.4.init=0.000001
	#==> stochastic freq is sometimes higher than det.
	#S.1.init=S.2.init=S.3.init=S.4.init=0.00001
	#==> stochastic freq is always higher than det., but not always by much
	S.0.init = 1                
	S.1.init = 0.001     
	S.2.init = 0.001
	S.3.init = 0.001
	S.4.init = 0.001
	B.0.init = 1
	P.S.0.init = 1
	P.B.0.init = 1
	
#res.gene.freq	= (S.1.init/4+S.2.init/2+S.3.init*3/4+S.4.init)/(S.1.init*3/4+S.2.init/2+S.3.init/4+S.0.init+B.0.init+P.S.0.init+P.B.0.init)
	
	#Place holders for parameter vector, parameters not currently used.	
##############################################################################	
place.holder1 = 0 		#place-holder previously used for a parameter	##
place.holder2 = 0 		#place-holder previously used for a parameter	##
##############################################################################
	
	
	return(c(h,predation,s.1J,s.1N,s.1A,s.2J,s.2N1,s.2A,s.2N2,s.3J,place.holder1,s.4J,place.holder2,g.1temp,g.2temp,alpha.3A,delta.3A,delta.4A,m.SS,k.w.SS,k.c.SS,m.BS,k.w.BS,k.c.BS,m.SB,k.w.SB,k.c.SB,m.BB,k.w.BB,k.c.BB,d.c[1],d.c[2],d.c[3],d.c[4],d.c[5],gamma,sigma[1],sigma[2],sigma[3],sigma[4],sigma[5],tau[1],tau[2],tau[3],tau[4],tau[5],mu,S.0.init,B.0.init,P.S.0.init,P.B.0.init,S.1.init, S.2.init,S.3.init,S.4.init))
}