###################################################################################
#Title: msm (Method of Simulated Moments function)
#Author: Lindsey Dietz
#Last Updated: 11/12/13
#Description:
#Function to produce the estimates given :
#family- exponential family of the top layer in the hierarchical model
#nsim- Number of simulations of the MSM average
#num_MC_sims-Number of values to produce one MSM
#num_subs-Index of i (number of subjects)
#obs_per_sub-Index of j (number of observations per subject)
#y.i-Sums over j of the y_ij, produced by sim.data.fun or provided by user
#start- Vector of starting values for (mu,sigma)
#true.mu- True value of mu provided by user
#true.sigma- True value of sigma provided by user
#method- one of three methods which include ”multiroot”,”optim”,”nleqslv”
#message-logical to display messages
###################################################################################
msm<-function(family="binomial", nsim=1, num_MC_sims =10000, num_subs =NULL, obs_per_sub =NULL,y.i=NULL,start=c(0,1),message=TRUE,true.mu=NULL, true.sigma=NULL,method="nleqslv"){

	#Confirming the user has entered all required arguments
	stopifnot(!is.null(num_subs))
	stopifnot(!is.null(obs_per_sub))
	stopifnot(!is.null(y.i))
	stopifnot(length(y.i)== num_subs)

	#Vectors to store estimates for mu, sigma, sigma^2 and mses of each
	mu<-vector(length=nsim)
	mu.mse<-vector(length=nsim)
	sigma2<-vector(length=nsim)
	sigma<-vector(length=nsim)
	sigma.mse<-vector(length=nsim)
	sigma2.mse<-vector(length=nsim)

	#Running simulations
#Currently only binomial and poisson are supported, this will be extended to other exponential families in future
 
 if (family=="binomial" || family=="poisson")
   {
   if(family=="binomial"){
 		for(i in 1: nsim){
		pars<-solver.sim(num_MC_sims= num_MC_sims, num_subs= num_subs, obs_per_sub= obs_per_sub,y.i= y.i,start=start, true.mu, true.sigma,method=method,family="binomial")

		mu[i]<-pars$par.1.mu
		mu.mse[i]<-pars$mu.mse
		sigma[i]<-pars$par.1.sigma
		sigma.mse[i]<-pars$sigma.mse
		sigma2[i]<-pars$par.1.sigma2
		sigma2.mse[i]<-pars$sigma2.mse
		remove(pars)
	}
  }
  
  if(family=="poisson"){
	for(i in 1: nsim){
		pars<-solver.sim(num_MC_sims= num_MC_sims, num_subs= num_subs, obs_per_sub= obs_per_sub,y.i= y.i,start=start, true.mu, true.sigma,method=method,family="poisson")

		mu[i]<-pars$par.1.mu
		mu.mse[i]<-pars$mu.mse
		sigma[i]<-pars$par.1.sigma
		sigma.mse[i]<-pars$sigma.mse
		sigma2[i]<-pars$par.1.sigma2
		sigma2.mse[i]<-pars$sigma2.mse
		remove(pars)
	}
  }
  
  #Finding the means and ses of the estimates of mu,sigma, and sigma^2
	
	#If no starting values provided for both mu and sigma, se is estimated by sample sd of estimates
	if(is.null(true.mu) && is.null(true.sigma)){
		 if(nsim>1){
		 if(message){print("nsim>1; Simulation SE= sd(parameter est.)/sqrt(nsim)")}	
		 mu.mse.est<-sd(mu.mse)/sqrt(length(mu.mse))
		 sigma2.mse.est<-sd(sigma2.mse)/sqrt(length(mu.mse))
		 sigma.mse.est<-sd(sigma.mse)/sqrt(length(mu.mse))
		 }
		 else {
		 if(message){print("nsim=1; No standard errors")}	
		 mu.mse.est<-NA
		 sigma2.mse.est<-NA
		 sigma.mse.est<-NA
		 }
	}
	#If there are starting values provided for both mu and sigma, se is calculated using true value	
	else{
		 if(nsim>1){
		 if(message){print("nsim>1 with true values; Simulation SE= sqrt(MSE)/sqrt(nsim)")}	
		 mu.mse.est<-sd(mu.mse)/sqrt(length(mu.mse))
		 sigma2.mse.est<-sd(sigma2.mse)/sqrt(length(sigma2.mse))
		 sigma.mse.est<-sd(sigma.mse)/sqrt(length(sigma.mse))	
		 }
		 else {
		 if(message){print("nsim=1 with true values; Simulation SE=sqrt(MSE)")}	
		 mu.mse.est<-sqrt(mu.mse)
		 sigma2.mse.est<-sqrt(sigma2.mse)
		 sigma.mse.est<-sqrt(sigma.mse)
		 }
	}
	return(list(mu=mean(mu),mu.se=mu.mse.est,sigma=mean(sigma),sigma.se=sigma.mse.est,sigma2=mean(sigma2),sigma2.se=sigma2.mse.est))
	
   }
  
  else{
		print("ERROR: Family not supported")
	}
}