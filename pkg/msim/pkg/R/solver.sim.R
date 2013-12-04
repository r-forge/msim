###################################################################################
#Title: msm (Method of Simulated Moments function)
#Author: Lindsey Dietz
#Last Updated: 11/12/13
#Description:
#Function to solve the system of equations given:
#family- exponential family of the top layer in the hierarchical model
#num_MC_sims-Number of MC simulated values to produce one MSM
#num_subs-Index of i, scalar number of subjects
#obs_per_sub-Index of j, vector of number of observations per subject
#y.i-Sums over j of the y_ij, produced by sim.data.fun or provided by user
#start- Vector of starting values for (mu,sigma)
#true.mu -True value of mu provided by user
#true.sigma -True value of sigma provided by user
#method- one of three methods which include ”multiroot”,”optim”,”nleqslv”
###################################################################################
solver.sim<-function(num_MC_sims =10000, num_subs=NULL,obs_per_sub=NULL,y.i=NULL,start=c(0,1), true.mu=NULL, true.sigma=NULL,method="nleqslv",family="binomial"){
	
	#Confirming the user has entered all required arguments
	stopifnot(!is.null(num_subs))
	stopifnot(!is.null(obs_per_sub))
	stopifnot(length(obs_per_sub)== num_subs)
	stopifnot(!is.null(y.i))

library(nleqslv)
library(rootSolve)	

	#Generate num_MC_sims eta's following N(0,1)
	eta<-rnorm(num_MC_sims,0,1)
		
	#Intial values of mu and sigma provided by user
	start <- start
	
	
  if (family=="binomial"){

    if (method=="nleqslv"){
		NRsolver<-function(x){
			mu<-x[1]
			sigma<-x[2]
			
			eq<-numeric(2)
			eq[1]<-(
					(1/num_MC_sims) * sum(exp(mu + sigma*eta)/(exp(mu + sigma*eta)+1))
					-(1/num_subs) * sum(y.i/obs_per_sub)
					)
			eq[2]<-(
					(1/num_MC_sims) * sum((exp(mu + sigma*eta)/(exp(mu + sigma*eta)+1))^2)
					-(1/num_subs) * sum((y.i^2-y.i)/(obs_per_sub*(obs_per_sub-1)))
					)
			eq
		}
	#Solve the system of equations and store parameters
	pars<-nleqslv(start, NRsolver,method="Newton")
	par.1.mu<-pars$x[1]
	par.1.sigma<-pars$x[2]
	par.1.sigma2<-par.1.sigma^2
	}
	
    if (method=="multiroot"){
		NRsolver<-function(x){
			mu<-x[1]
			sigma<-x[2]
			
			eq<-numeric(2)
			eq[1]<-(
					(1/num_MC_sims) * sum(exp(mu + sigma*eta)/(exp(mu + sigma*eta)+1))
					-(1/num_subs) * sum(y.i/obs_per_sub)
					)
			eq[2]<-(
					(1/num_MC_sims) * sum((exp(mu + sigma*eta)/(exp(mu + sigma*eta)+1))^2)
					-(1/num_subs) * sum((y.i^2-y.i)/(obs_per_sub*(obs_per_sub-1)))
					)
			eq
		}
	#Solve the system of equations and store parameters
	pars<-multiroot(f=NRsolver, start= start)
	par.1.mu<-pars$root[1]
	par.1.sigma<-pars$root[2]
	par.1.sigma2<-par.1.sigma^2
	}
		
	else{
		euclid.norm.sq<-function(par){
		mu<-par[1]
		sigma<-par[2]
		(((1/num_MC_sims)*sum(exp(mu+sigma*eta)/(exp(mu+sigma*eta)+1))-(1/num_subs)*sum(y.i/obs_per_sub)))^2+(((1/num_MC_sims)*sum((exp(mu+sigma*eta)/(exp(mu+sigma*eta)+1))^2)-(1/num_subs)*sum((y.i ^2-y.i)/(obs_per_sub*(obs_per_sub-1)))))^2
		}
	#Solve the equations and store parameters
	pars<-optim(start, fn=euclid.norm.sq)
	par.1.mu<-pars$par[1]
	par.1.sigma<-pars$par[2]
	par.1.sigma2<-par.1.sigma^2
	}
}
  if (family=="poisson"){

    if (method=="nleqslv"){
		NRsolver<-function(x){
			mu<-x[1]
			sigma<-x[2]
			
			eq<-numeric(2)
			eq[1]<-(
					(1/num_MC_sims) * sum(exp(mu + sigma*eta))
					-(1/num_subs) * sum(y.i/obs_per_sub)
					)
			eq[2]<-(
					(1/num_MC_sims) * sum((exp(mu + sigma*eta))^2)
					-(1/num_subs) * sum((y.i^2-y.i)/(obs_per_sub^2))
					)
			eq
		}
	#Solve the system of equations and store parameters
	pars<-nleqslv(start, NRsolver,method="Newton")
	par.1.mu<-pars$x[1]
	par.1.sigma<-pars$x[2]
	par.1.sigma2<-par.1.sigma^2
	}
	
    if (method=="multiroot"){
				NRsolver<-function(x){
			mu<-x[1]
			sigma<-x[2]
			
			eq<-numeric(2)
			eq[1]<-(
					(1/num_MC_sims) * sum(exp(mu + sigma*eta))
					-(1/num_subs) * sum(y.i/obs_per_sub)
					)
			eq[2]<-(
					(1/num_MC_sims) * sum((exp(mu + sigma*eta))^2)
					-(1/num_subs) * sum((y.i^2-y.i)/(obs_per_sub^2))
					)
			eq
		}
	#Solve the system of equations and store parameters
	pars<-multiroot(f=NRsolver, start= start)
	par.1.mu<-pars$root[1]
	par.1.sigma<-pars$root[2]
	par.1.sigma2<-par.1.sigma^2
	}
		
	else{
		euclid.norm.sq<-function(par){
		mu<-par[1]
		sigma<-par[2]
		(
		(((1/num_MC_sims) * sum(exp(mu + sigma*eta))
					-(1/num_subs) * sum(y.i/obs_per_sub)))^2
		+(((1/num_MC_sims) * sum((exp(mu + sigma*eta))^2)
					-(1/num_subs) * sum((y.i^2-y.i)/(obs_per_sub^2))))^2
		)
		}
	#Solve the equations and store parameters
	pars<-optim(start, fn=euclid.norm.sq)
	par.1.mu<-pars$par[1]
	par.1.sigma<-pars$par[2]
	par.1.sigma2<-par.1.sigma^2
	}
 }

	#If no true parameters were entered, the mse will not be calculated here
	if(is.null(true.mu) && is.null(true.sigma)){
		mu.mse<-par.1.mu
		sigma2.mse<-par.1.sigma2
		sigma.mse<-par.1.sigma
	}
	
	#If true parameters were entered, the mse will be calculated here
	else{
		mu.mse<-(par.1.mu-true.mu)^2
		sigma2.mse<-(par.1.sigma2-true.sigma^2)^2
		sigma.mse<-(par.1.sigma-true.sigma)^2
		}
		
return(list(par.1.mu=par.1.mu,mu.mse=mu.mse,par.1.sigma=par.1.sigma,par.1.sigma2=par.1.sigma2,sigma2.mse=sigma2.mse,sigma.mse=sigma.mse))	
}
  