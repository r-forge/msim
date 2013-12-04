###################################################################################
#Title: boot.msm (Bootstrapping Bias Correction)
#Author: Lindsey Dietz
#Last Updated: 11/5/13
#Description:
#Function to simulate the data for the experiment given:
#msm.est- list generated from msm() function
#boot.sim- number of bootstrap samples
#nsim-number of times to compute the estimates for one bootstrap sample
#num_MC_sims-number of MC samples for estimating expected value of moments
#num_subs-Index of j, scalar number of subjects
#obs_per_sub-Index of i, vector of number of observations per subject
#start- starting values of mu, sigma
#boot.message- logical to display iteration counts
###################################################################################
boot.msm<-function(msm.est=NULL, boot.message=TRUE ,boot.sim=10, family="binomial",nsim=1, num_MC_sims = 10000, num_subs =NULL, obs_per_sub =NULL,start=c(0,1),method="nleqslv"){
	
	#Confirming the user has entered all required arguments
	stopifnot(!is.null(num_subs))
	stopifnot(!is.null(obs_per_sub))
	stopifnot(length(obs_per_sub)==num_subs)
	stopifnot(!is.null(msm.est))
	
	#Vectors to contain the estimates for each bootstrap sample
	vect1<-vector(length= boot.sim)
	vect2<-vector(length= boot.sim)
	vect3<-vector(length= boot.sim)

	for(i in 1:boot.sim){
	if(boot.message){print(paste("Bootstrap",i))}
	#Simulate data from the distribution with the parameter estimates from msm object	
ys<-sim.data.fun(num_subs= num_subs,obs_per_sub= obs_per_sub,true.mu= msm.est$mu,true.sigma= msm.est$sigma)

	#Run msm for the bootstrap sample
	ms2<-msm(family= family, nsim= nsim,num_MC_sims= num_MC_sims, num_subs= num_subs,obs_per_sub= obs_per_sub,y.i= ys$y.i ,start= start,true.mu= msm.est$mu,true.sigma= msm.est$sigma, message=FALSE,method= method)
	
	#Correct the bias in the original estimate by using the bootstrap estimate
	vect1[i]<- 2* msm.est$mu-ms2$mu
	vect2[i]<- 2* msm.est$sigma-ms2$sigma
	vect3[i]<- vect2[i]^2
	}
	if(boot.message){print("Bootstrap SE is used")}

	#Calculate estimates of parameters and standard errors
	boot.mu<-mean(vect1)
	boot.mu.se<-sd(vect1)/sqrt(boot.sim)
	boot.sigma <-mean(vect2)
	boot.sigma.se<-sd(vect2)/sqrt(boot.sim)
	boot.sigma2<-mean(vect3)
	boot.sigma2.se<-sd(vect3)/sqrt(boot.sim)
	return(list(boot.mu= boot.mu, boot.mu.se= boot.mu.se, boot.sigma= boot.sigma, boot.sigma.se= boot.sigma.se, boot.sigma2= boot.sigma2,boot.sigma2.se=boot.sigma2.se))
}
