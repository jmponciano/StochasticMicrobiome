#  Ornstein Uhlenbeck State Space Model version date 122010.R:
#
#  Program for calculating maximum likelihood (ML) or restricted maximum
#  likelihood (REML) estimates of unknown parameters for the Ornstein
#  Uhlenbeck State Space (OUSS) model of stochastic population growth.  
#  The model is
#
#  dX(t) = theta*[mu - X(t)]*dt  +  dW(t)
#               with dW(t) ~ normal(0,dt*betasq),
#  Y(t) = X(t) + F(t)
#               with F(t) ~ normal(0,tausq).
#
#  Here X(t) is log-population abundance, Y(t) is observed or estimated value
#  of X(t), theta, mu, betasq, tausq are parameters.  The parameter betasq
#  is the variance of the process noise, and tausq is the variance of the
#  observation error.
#
#  The model takes population abundance N(t) = exp(X(t)) to be governed by a
#  stochastic, continuous time density dependent model, with the observed
#  abundances O(t) = N(t)*exp(F(t)) arising from lognormal sampling error.
#
#  User provides time series of observed population abundances o(0), o(1),
#  ..., o(q), which are log-transformed by the program into y(0), y(1), ...,
#  y(q), assumed to be a time series realization of Y(t).  Likelihood
#  function of y(0), y(1), ..., y(q) is that of a multivariate normal
#  distribution.  The observation times t_0, t_1, t_2, ..., t_q can have
#  unequal intervals.
#
#  Program computes initial parameter values for iterations.  The program
#  should be re-run for several sets of initial values, as the likelihood
#  function for the model frequently has multiple local maxima.
#
#  This program written by Brian Dennis (Dept Fish and Wildlife Resources,
#  Univ Idaho, Moscow, Idaho, 83844-1136 USA  brian@uidaho.edu).
#
#  Citations:
#    Dennis, B. 2010.  A density dependent state space model for population 
#    abundance data with unequal time intervals.  Unpublished manuscript for:
#    Ecology.
#
#----------------------------------------------------------------------
#        USER INPUT SECTION
#----------------------------------------------------------------------
#  User supplies time series data here.  User can substitute R statements to
#    read population abundance data from a file into the vector "Observed.t".
#    Times of observation are entered into the vector "Time.t".

Observed.t=c(346,675,802,1478,1173,756,861,972,854,1161,1318,901,901,1173,
   608,811,903,584,1179,1020,1129,966);  #  No zeros!  (With zeros, you must
                                         #    use another model).
Time.t=c(1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1970,1971,1972,
   1973,1974,1975,1976,1977,1978,1979,1980,1981);  #  Initial time can 
                                                   #    be nonzero.

#  Example data are bobcat (Lynx rufus) in Idaho, data set 212 from the Global
#  Population Dynamics Database.
#----------------------------------------------------------------------
#        PROGRAM INITIALIZATION SECTION
#----------------------------------------------------------------------  
library(MASS);  #  loads miscellaneous functions (ginv, etc.)
T.t=Time.t-Time.t[1]; #  For calculations, time starts at zero.
Y.t=log(Observed.t);  #  Log-transform the observations.
q=length(Y.t)-1;      #  Number of time series transitions, q.
qp1=q+1;              #  q+1 gets used a lot, too.
S.t=T.t[2:qp1]-T.t[1:q];  #  Time intervals.

#----------------------------------------------------------------------
#        SECTION FOR DEFINING ML & REML LOG-LIKELIHOODS
#----------------------------------------------------------------------

#  ML objective function "negloglike.OU.ml" is negative of log-likelihood;
#  the Nelder-Mead optimization routine in R, "optim", is a minimization
#  routine.  The ML objective function uses a multivariate normal log-
#  likelihood (Eq. 16 in Dennis 2010).  The three function arguments are:
#  phi, vector of parameters (transformed to the real line), 
#  yt, vector of time series observations (log scale),
#  tt, vector of observation times.

negloglike.OU.ml=function(phi,yt,tt)  
{
   mu=phi[1];
   theta=exp(phi[2]);         #  Constrains theta > 0.
   betasq=exp(phi[3]);      #  Constrains betasq > 0. 
   tausq=exp(phi[4]);       #  Constrains tausq > 0.
   thing=betasq/(2*theta);    #  Recurring quantity.
   q=length(yt)-1;
   ss=tt[2:qp1]-tt[1:q];      #  Time intervals.
   qp1=q+1;
   vx=matrix(1,qp1,qp1);  # Preallocate matrix for autocorrelations.
   #  Following loop calculates the model autocorrelations and puts 
   #    them in vx.
   for (ti in 1:q)
   {
      vx[(ti+1):qp1,ti]=exp(-theta*cumsum(ss[ti:q]));
      vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
   }
   Sigma.mat=vx*thing;   #  Variance-covariance matrix for X(t).
   #print(Sigma.mat)
   Itausq=matrix(0,qp1,qp1);
   diag(Itausq)=rep(tausq,qp1);
   V=Sigma.mat+Itausq;   #  Variance-covariance matrix for Y(t).
   mm=rep(mu,qp1);
   ofn=(qp1/2)*log(2*pi)+0.5*log(det(V))+
      0.5*(yt-mm)%*%ginv(V)%*%(yt-mm);
   return(ofn);
}

#  REML objective function "negloglike.OU.reml" is negative of log-likelihood
#  for first differences of the log-scale observations.  The REML objective
#  function uses equations 17-19 of Dennis (2010).  The three
#  function arguments are: 
#  phi, vector of parameters (transformed to the real line), 
#  yt, vector of time series observations (log scale),
#  tt, vector of observation times.
#  The function performs the differencing.
negloglike.OU.reml=function(phi,yt,tt)
{
   theta=exp(phi[1]);          #  Constrains th > 0.
   betasq=exp(phi[2]);      #  Constrains betasq > 0. 
   tausq=exp(phi[3]);       #  Constrains tausq > 0.
   thing=betasq/(2*theta);       #  Recurring quantity.
   ss=tt[2:qp1]-tt[1:q];         #  Time intervals.
   q=length(yt)-1;
   qp1=q+1;
   vx=matrix(1,qp1,qp1);  # Preallocate matrix for autocorrelations.
   #  Following loop calculates the autocorrelations and puts 
   #    them in vx.
   for (ti in 1:q)
   {
      vx[(ti+1):qp1,ti]=exp(-theta*cumsum(ss[ti:q]));
      vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
   }
   Sigma.mat=vx*thing;
   Itausq=matrix(0,qp1,qp1);
   diag(Itausq)=rep(tausq,qp1);
   V=Sigma.mat+Itausq;
   Dmat=cbind(-diag(1,q),matrix(0,q,1))+
      cbind(matrix(0,q,1),diag(1,q));    #  Differencing matrix.
   Phi.mat=Dmat%*%V%*%t(Dmat);           #  REML var-cov mattrix.
   wt=yt[2:qp1]-yt[1:q];
   ofn=(q/2)*log(2*pi)+0.5*log(det(Phi.mat))+
      0.5*wt%*%ginv(Phi.mat)%*%wt;
   return(ofn);
}

#----------------------------------------------------------------------
#        SECTION FOR CALCULATING INITIAL VALUES.
#----------------------------------------------------------------------
Ybar=mean(Y.t);
Yvar=sum((Y.t-Ybar)*(Y.t-Ybar))/q;

mu1=Ybar;
th1=-mean(log(abs((Y.t[2:qp1]-mu1)/
   (Y.t[1:q]-mu1)))/S.t);            # Kludge an initial value for theta
                                     #  based on mean of Y(t+s) given Y(t).
bsq1=2*th1*Yvar/(1+2*th1);         #  Moment estimate using stationary
tsq1=bsq1;                         #   variance, with betasq=tausq.

#----------------------------------------------------------------------
#        SECTION FOR CALCULATING ML & REML PARAMETER ESTIMATES
#----------------------------------------------------------------------

# The ML estimates.
OUSSml=optim(par=c(mu1,log(th1),log(bsq1),log(tsq1)),
   negloglike.OU.ml,NULL,method="Nelder-Mead",yt=Y.t,tt=T.t);
params.ml=c(OUSSml$par[1],exp(OUSSml$par[2]),exp(OUSSml$par[3]),
   exp(OUSSml$par[4]));
lnlike.ml=-OUSSml$value[1];
AIC.ouss=-2*lnlike.ml+2*length(params.ml);

mu.ml=params.ml[1];           # These are the ML estimates.
theta.ml=params.ml[2];        #          --
betasq.ml=params.ml[3];       #          --
tausq.ml=params.ml[4];        #          --

# The REML estimates.
OUSSreml=optim(par=c(log(th1),log(bsq1),log(tsq1)),
   negloglike.OU.reml,NULL,method="Nelder-Mead",yt=Y.t,tt=T.t);
params.reml=c(exp(OUSSreml$par[1]),exp(OUSSreml$par[2]),
   exp(OUSSreml$par[3]));
theta.reml=params.reml[1];    #  These are the REML estimates.
betasq.reml=params.reml[2];   #           --
tausq.reml=params.reml[3];    #           --

#  Calculate REML estimate of mu using Eq. 20 of Dennis (2010).
thing=betasq.reml/(2*theta.reml);
vx=matrix(1,qp1,qp1);
for (ti in 1:q)
{
   vx[(ti+1):qp1,ti]=exp(-theta.reml*cumsum(S.t[ti:q]));
   vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
}
Sigma.mat=vx*thing;
Itausq=matrix(0,qp1,qp1);
diag(Itausq)=rep(tausq.reml,qp1);
V.reml=Sigma.mat+Itausq;
j=matrix(1,qp1,1);
Vinv=ginv(V.reml);
mu.reml=(t(j)%*%Vinv%*%Y.t)/(t(j)%*%Vinv%*%j);  #  REML estimate of mu.

Var_mu.reml=1/(t(j)%*%Vinv%*%j);        #  Variance of mu
mu_hi.reml=mu.reml+1.96*sqrt(Var_mu.reml); #  95% CI for mu
mu_lo.reml=mu.reml-1.96*sqrt(Var_mu.reml); #       --

#  Calculate predicted population sizes for the OUSS model
#  (X(tj) given all the observations except for Y(tj))
#  with multivariate normal distribution, for plotting.
#
#  Choose ML or REML estimates here (by commenting out the unwanted).
#  mu=mu.ml;  theta=theta.ml;  betasq=betasq.ml;  tausq=tausq.ml;
mu=mu.reml;  theta=theta.reml;  betasq=betasq.reml;  tausq=tausq.reml;

thing=betasq/(2*theta);
vx=matrix(1,qp1,qp1);
for (ti in 1:q)
{
   vx[(ti+1):qp1,ti]=exp(-theta*cumsum(S.t[ti:q]));
   vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
}
Sigma.mat=vx*thing;
Itausq=matrix(0,qp1,qp1);
diag(Itausq)=rep(tausq,qp1);
V=Sigma.mat+Itausq;

Predict.t=rep(0,qp1);
Muvec=rep(mu,q);
for (tj in 1:qp1)
{
Y.omitj=Y.t[-tj];    #  Omit observation at time tj.
V.omitj=V[-tj,-tj];  #  Omit row tj and col tj from var-cov matrix.
V12=V[tj,-tj];       #  Submatrix:  row tj without col tj.
Predict.t[tj]=mu+V12%*%ginv(V.omitj)%*%(Y.omitj-Muvec);  #  Usual expression
                                                         #  for conditional
                                                         #  MV normal mean.
}
Predict.t=exp(Predict.t);

#  Plot the data & model-fitted values
plot(Time.t,Observed.t,xlab="time",ylab="population abundance",
   type="o",pch=1,cex=1.5);      #  Population data are circles.
par(lty="dashed");               #  Predicted abundances are dashed line.
points(Time.t,Predict.t, type="l", lwd=1);

#  Print the parameter estimates
parms.reml=c(mu.reml,theta.reml, betasq.reml,tausq.reml); #  Collect for 
                                                          #    printing.
parms.ml=c(mu.ml,theta.ml, betasq.ml,tausq.ml);           #      --
names=c("mu","theta","betasq","tausq");                   #      --
types=c("OUSS-ML","OUSS-REML");                           #      --
matrix(cbind(parms.ml,parms.reml),
   nrow=2,ncol=4,byrow=TRUE,dimnames=list(types,names));       #  Print stuff

matrix(cbind(mu_lo.reml,mu_hi.reml),nrow=1,ncol=2,byrow=TRUE,
   dimnames=list("95% CI for MU",c("LO","HI")));               #  Print stuff

matrix(cbind(lnlike.ml,AIC.ouss),nrow=1,ncol=2,byrow=TRUE,
   dimnames=list("OUSS ML RESULTS",c("LN-LIKELIHOOD","AIC"))); #  Print stuff
