##############################################
#Title: Inverse-Logit function         		 #
#Author: Lindsey Dietz          			 #
#Last Modified: 4/8/13					     #
#Description: Computes f(x)=exp(x)/(1+exp(x))#
##############################################
inv.logit<-function(x){
	exp(x)/(1+exp(x))
}