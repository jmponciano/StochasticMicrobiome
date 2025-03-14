#----------------------------------------------------------------------
#  Load all the OUSS model functions and needed packages
#----------------------------------------------------------------------
library("MASS")
source("ROUSSE-1.0.R")

#----------------------------------------------------------------------
#        USER INPUT SECTION
#----------------------------------------------------------------------
#  User supplies time series data here into the vecor "Observed.t.  
#  User can substitute R statements to read population abundance data 
#  from a file into the vector "Observed.t". Do not change the object name "Observed.t"
#  Times of observation are entered into the vector "Time.t". Do not change the object name "Time.t"

#  Example data are bobcat (Lynx rufus) in Florida and in Idaho, data sets 211 and 212 from the Global
#  Population Dynamics Database. Pick one of these two data sets and associated sampling years to run the example

#  Linx rufus, from Florida.  GPPD data set 211
Time.t <- c(1946,1947,1948,1949,1950,1954,1955,1956,1957,1958,1959,1960,1961,1963,1964,1965,1966,1967,1968,1975,1976,1977,1978,1979,1980,1981);
Observed.t <- c(672,1028,538,566,300,400,400,400,400,300,250,450,450,13,23,23,2,400,20,389,537,983,1698,1132,1702,1031);

# Linx rufus, from Idaho, GPPD data set 212
Observed.t=c(346,675,802,1478,1173,756,861,972,854,1161,1318,901,901,1173, 608,811,903,584,1179,1020,1129,966);  #  No zeros! 
Time.t=c(1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981); 

# Population growth in Guatemala city
Observed.t=c(1459368,1504045,1550224,1596787,1646563,1696391,1747542,1800076,1853691,1908085,1962953,2246440,2538227,2580256,2641473,2702257,2762328,2821400,2879664,2937307,	2994047,3049601,3103685);  #  No zeros! 

Observed.t<-c(393343,404632,416273,428015,440453,453854,467321,481323,495794,510733,526249,610395,704731,550772.3197,565109.111,579749.5316,594577.7391,609477.8912,624527.2162,639802.9422,	655189.227,670570.2282,685830.1035)
Time.t=c(1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1995, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010); 


plot(Time.t, Observed.t, pch=16)

#-------- Log-transform the observations to carry all the calculations in this program--------------------#
log.obs    <- log(Observed.t); 
#---------------------------------------------------------------------------------------------------------#



#----------------------------------------------------------------------
#        PARAMETER ESTIMATION, PARAMETRIC BOOTSTRAP AND PREDICTIONS
#----------------------------------------------------------------------
# Before doing the calculations, the user has to specify ONLY the following 4 options:

# 1. Do you want to compute the ML estimates or the REML estimates?

method <- "REML" # alternatively, set method <- "REML"

# 2. Do you want to plot the predictions?
pred.plot <- "TRUE" # Set it to "FALSE" if you do not want to plot the predictions

# 3. Do you want to plot the parametric bootstrap distribution of the estimates?
pboot.plot <- "TRUE" # Set it to "FALSE" if you do not want to plot the bootstrap distribution of the estimates

# 4. How many bootstrap replicates?
NBoot <- 1000; 


#-------------------------------------------------------------------------------------------------#
#  5. RUN THE FOLLOWING LINE OF CODE TO COMPUTE THE ESTIMATES, PREDICTIONS,
#     AND RUN A PARAMETRIC BOOTSTRAP. THE USER DOES NOT NEED TO MODIFY THIS LINE OF CODE.
#     THE OUTPUT OF THE FUNCTION 'ROUSS.CALCS' IS A LIST AND THE USER CAN RETRIEVE EACH OF THE   
#     LIST ELEMENTS PRINTED AND SAVED IN THE OBJECT "all.results". 
#     THE 95% PARAMETRIC BOOTSTRAP FOR BOTH, THE PARAMETERS AND THE PREDICTIONS ARE COMPUTED 
#     BY THE FUNCTION "ROUSS.CALCS" 	
#      
#-------------------------------------------------------------------------------------------------#

all.results <- ROUSS.CALCS(Yobs=log.obs,Tvec=Time.t, pmethod=method, nboot=NBoot, plot.pred=pred.plot, plot.bootdists = pboot.plot);


