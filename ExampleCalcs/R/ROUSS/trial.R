# Trial

Time.t <- c(1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1970,1971,1972,
   1973,1974,1975,1976,1977,1978,1979,1980,1981);
tt  <- Time.t-Time.t[1];
q   <- lenght(tt) -1; 
qp1 <- q+1;
theta <- 1.5

# This is from your function:
vx=matrix(1,qp1,qp1);  # Preallocate matrix for autocorrelations.
#  Following loop calculates the model autocorrelations and puts 
#    them in vx.
for (ti in 1:q)
   {
      #print(cumsum(ss[ti:q]))
      vx[(ti+1):qp1,ti]= exp(-theta*cumsum(ss[ti:q]));
      vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
   }


#Here's a faster way of doing the same (useful to speed up optimization)

t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE)
t.rows     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=TRUE)
abs.diffs  <- abs(t.rows-t.cols)
theothervx <- exp(-theta*abs.diffs)

vx-theothervx # 0 indeed.