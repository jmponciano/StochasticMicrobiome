###  Parametric Bootstrap Likelihood Ratio test, an example using
###  the BBS time series from the North American Breeding Bird Survey,
###  record # 02014 3328 08636, 1966-95 (Table 1 in Dennis et al 2006)

#----------------------------------------------------------------------
#  Load all the OUSS model functions and needed packages
#----------------------------------------------------------------------
library("MASS")
source("ROUSSE-1.0.R")


#----------------------------------------------------------------------
#        USER INPUT SECTION
#----------------------------------------------------------------------
#  User supplies time series data here into the vector "Observed.t.  
#  User can substitute R statements to read population abundance data 
#  from a file into the vector "Observed.t". Do not change the object name "Observed.t"
#  Times of observation are entered into the vector "Time.t". Do not change the object name "Time.t"

Observed.t <- c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6);
Time.t    <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29);


#-------- Log-transform the observations to carry all the calculations in this program--------------------#
log.obs    <- log(Observed.t); 
#---------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------
#        PARAMETER ESTIMATION AND PARAMETRIC BOOTSTRAP LIKELIHOOD RATIO TEST
#----------------------------------------------------------------------
#  The output contains 6 objects:
#  Object 1:  "egssml": this is a list that contains the ml estimates of the EGSS model, the maximized log likelihod
#             under that model and the AIC value under that model  
#
#  Object 2:  "roussml": this is a list that contains the ml estimates of the OUSS model, the maximized log likelihood
#             under that model and the AIC value under that model 
#
#  Object 3:  "Lam.obs":  this is the observed -2*ln[L(Ho/L(H1))]
#
#  Object 4:  "Lam.vec":  this is a vector of the bootstrapped -2*ln[L(Ho/L(H1))] values
#
#  Object 5:  "pvalue":  this is the proportion of the bootstrap simulations that the bootstrapped -2*ln[L(Ho/L(H1))] 
#                        is more extreme than the observed -2*ln[L(Ho/L(H1))] 
#
#  Object 6:  "Decision.rule": this is a character object that prints out the decision of the test (Fail to reject or
#              reject the null hypothesis of density independence)
#
# Before doing the calculations, the user has to specify ONLY the following 3 options:

# 1. How many bootstrap replicates do you want to use? (Default is a thousand)
NBoot  <- 1000;

# 2. At what level alpha do you want to test the null hypothesis of density independence vs. density dependence?
my.alpha <- 0.05;

# 3. Do you want to plot the sampling distribution of the log likelihood ratio -2*ln[L(Ho/L(H1))]

plot.it <- FALSE

pblrt.trial <- PBLRT.ddpcont(B=NBoot, Yobs=log.obs, Tvec=Time.t, alpha=my.alpha, plot.PBLR=plot.it)

# Printing results:

#  Object 1:  "egssml": this is a list that contains the ml estimates of the EGSS model, the maximized log likelihod
#             under that model and the AIC value under that model  
pblrt.trial$egssml

#  Object 2:  "roussml": this is a list that contains the ml estimates of the OUSS model, the maximized log likelihood
#             under that model and the AIC value under that model 
pblrt.trial$roussml

#  Object 3:  "Lam.obs":  this is the observed -2*ln[L(Ho)/L(H1)]
pblrt.trial$Lam.obs

#  Object 5:  "pvalue":  this is the proportion of the bootstrap simulations that the bootstrapped -2*ln[L(Ho)/L(H1)] 
#                        is more extreme than the observed -2*ln[L(Ho)/L(H1)] 
pblrt.trial$pvalue

#  Object 6:  "Decision.rule": this is a character object that prints out the decision of the test (Fail to reject or
#              reject the null hypothesis of density independence)
pblrt.trial$Decision.rule

#  To plot the sampling distribution of the likelihood ratio test, I recommend to do a trial histogram first, and after
#  a visual inspection, eliminate from the plotting range the numerical extremes, which should be very few.
#  here is an example with this particular data set

#  This is the histogram of the pboot LRT with the observed LRT in red
hist(pblrt.trial$Lam.vec);abline(v=pblrt.trial$Lam.obs,lwd=2,col="red")

# After plotting, note that we have some numerical extremes around -600 and around +600.  
# Let's re-do the plot zooming into the bulk of the bootstrap distribution:

hist(pblrt.trial$Lam.vec[abs(pblrt.trial$Lam.vec)<100], xlab="Bootstrap values of -2*ln[L(Ho)/L(H1)]", main="");
abline(v=pblrt.trial$Lam.obs,lwd=2,col="red")


