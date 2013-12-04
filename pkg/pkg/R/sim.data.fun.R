###################################################################################
#Title: sim.data.fun (Simulating data function)
#Author: Lindsey Dietz
#Last Updated: 11/5/13
#Description:
#Function to simulate the data for the experiment given:
#num_subs-Index of i, scalar number of subjects
#obs_per_sub-Index of j, vector of number of observations per subject
#true.mu- True value of mu provided by user
#true.sigma- True value of sigma provided by user
###################################################################################
#Simulated Data Function
sim.data.fun<-function(num_subs =NULL, obs_per_sub =NULL,true.mu=NULL,true.sigma=NULL,family="binomial", set.seed=NULL){
	
	#Confirming the user has entered all required arguments
	stopifnot(!is.null(num_subs))
	stopifnot(!is.null(obs_per_sub))
	stopifnot(length(obs_per_sub)== num_subs)
	stopifnot(!is.null(true.mu))
	stopifnot(!is.null(true.sigma))
	
	#Acknowledging if the user has set a seed
	if(!is.null(set.seed)){
		set.seed(set.seed) 
	}

  if(family=="binomial"){	
	logit.p.y<-vector(length= num_subs)
	#logit[P(y_ij=1|eta)]=true.mu+true.sigma*eta_i
	#eta_i~N(0,1)
	logit.p.y<-true.mu +true.sigma*rnorm(num_subs,0,1)

	#Matrix for random sample of y_ij's
	y<-matrix(ncol=max(obs_per_sub),nrow= num_subs)
	for(i in 1:num_subs){
		for(j in 1: obs_per_sub[i]){
			y[i,j]<-rbinom(1,1, inv.logit(logit.p.y[i]))
		}
	}
	
	#Vector y_ij's summed over j
	y.i<-vector(length= num_subs)
	n.i<-vector(length= num_subs)
	for(i in 1:num_subs){
		y.i[i]<-sum(y[i,],na.rm=T)
		n.i[i]<-max(obs_per_sub)-sum(is.na(y[i,]))
	}
  }
  
  if(family=="poisson"){	
	log.mu<-vector(length= num_subs)
	#logit[P(y_ij=1|eta)]=true.mu+true.sigma*eta_i
	#eta_i~N(0,1)
	log.mu<-true.mu +true.sigma*rnorm(num_subs,0,1)

	#Matrix for random sample of y_ij's
	y<-matrix(ncol=max(obs_per_sub),nrow= num_subs)
	for(i in 1:num_subs){
		for(j in 1: obs_per_sub[i]){
			y[i,j]<-rpois(1,exp(log.mu[i]))
		}
	}
	
	#Vector y_ij's summed over j
	y.i<-vector(length= num_subs)
	n.i<-vector(length= num_subs)
	for(i in 1:num_subs){
		y.i[i]<-sum(y[i,],na.rm=T)
		n.i[i]<-max(obs_per_sub)-sum(is.na(y[i,]))
	}
  }
return(list(y=y,y.i=y.i,n.i=n.i))	
}