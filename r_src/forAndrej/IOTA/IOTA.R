#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################



beta_calc_weight <- function(ts,BE,method)
{
 #return(sum((ts[BE[,2]+1] > ts[BE[,1]] & ts[BE[,2]] < ts[BE[,1]]) | (ts[BE[,2]+1] < ts[BE[,1]] & ts[BE[,2]] > ts[BE[,1]])))
 if(method == "slope") w <- abs(ts[BE[,2]+1]-ts[BE[,2]])
 if(method == "sqrt") w <- (ts[BE[,2]+1]-ts[BE[,2]])^2
 if(method == "am") w <- (ts[BE[,2]+1]+ts[BE[,2]])/2
 if(method == "gm") w <- sqrt(ts[BE[,2]+1]*ts[BE[,2]])
 if(method == "hm") w <- 2/((1/ts[BE[,2]+1])+(1/ts[BE[,2]]))
 if(sum(c("slope","sqrt","am","gm","hm") == method) == 0) print(" impemented weighting functions are only: 'slope', 'sqrt', 'am', 'gm' and 'hm' ")

 beta <- sum(w*((ts[BE[,2]+1] > ts[BE[,1]] & ts[BE[,2]] < ts[BE[,1]]) | (ts[BE[,2]+1] < ts[BE[,1]] & ts[BE[,2]] > ts[BE[,1]])))
 return(beta)
}
beta_calc <- function(ts,BE)
{
 beta <- sum((ts[BE[,2]+1] > ts[BE[,1]] & ts[BE[,2]] < ts[BE[,1]]) | (ts[BE[,2]+1] < ts[BE[,1]] & ts[BE[,2]] > ts[BE[,1]]))
 return(beta)
}
ordering <- function(ord,timeseries,BE,delta)
{
 return(1-(apply(timeseries[,ord],1,beta_calc,BE)/delta))
}
crossing2 <- function(timeseries,BE,delta)
{
 ord <- t(apply(timeseries,1,order,decreasing="FALSE"))
 return(t(apply(ord,1,ordering,timeseries,BE,delta)))
}

#beta_calc_weight <- function(ts,BE)
#{
# return(sum(((ts[BE[,2]+1]-ts[BE[,2]])^2)*((ts[BE[,2]+1] > ts[BE[,1]] & ts[BE[,2]] < ts[BE[,1]]) | (ts[BE[,2]+1] < ts[BE[,1]] & ts[BE[,2]] > ts[BE[,1]])))
#)
#}
ordering_weight <- function(ord,timeseries,BE,delta,method)
{
 return(1-(apply(timeseries[,ord],1,beta_calc_weight,BE,method)/delta))
}
crossing2_weight <- function(timeseries,BE,delta,method)
{
 ord <- t(apply(timeseries,1,order,decreasing="FALSE"))
 return(t(apply(ord,1,ordering_weight,timeseries,BE,delta,method)))
}

IOTA <- function(TimeSeries,method="both")
{

 #################################################
 # each row of TimeSeries should contain values at different time points for one variable 
 #################################################

 n <- length(TimeSeries[1,])
 delta <- (n-1)*(n-2)/2 

 COMP <- array(0,c(((n-1)*(n-2)/2),2))
 c1 <- 0
 for(t in 1:(n-2))
 {
  for(s in (t+1):(n-1))
  {
   c1 <- c1+1
   COMP[c1,1] <- t
   COMP[c1,2] <- s
  }
 }

 #################################################
 if(method == "both")
 {
  mu_org <- crossing2(TimeSeries,COMP,delta)
  mu_org_weight <- crossing2_weight(TimeSeries,COMP,delta,method="sqrt") 
  return(list(basic.measure = mu_org,weighted.measure = mu_org_weight))
 }
 if(method == "orginal")
 {
  mu_org <- crossing2(TimeSeries,COMP,delta)
  return(basic.measure = mu_org)
 } 
 if(method != "both" && method != "original")
 {
  mu_org_weight <- crossing2_weight(TimeSeries,COMP,delta,method) 
  return(weighted.measure = mu_org_weight)
 } 
}



#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################

crossing2_reverse <- function(timeseries,BE,delta)
{
 ord <- t(apply(timeseries,1,order,decreasing="TRUE"))
 return(t(apply(ord,1,ordering,timeseries,BE,delta)))
}
crossing2_weight_reverse <- function(timeseries,BE,delta,method)
{
 ord <- t(apply(timeseries,1,order,decreasing="TRUE"))
 return(t(apply(ord,1,ordering_weight,timeseries,BE,delta,method)))
}

IOTA_reverse <- function(TimeSeries,method="both")
{

 #################################################
 # each row of TimeSeries should contain values at different time points for one variable 
 #################################################

 n <- length(TimeSeries[1,])
 delta <- (n-1)*(n-2)/2 

 COMP <- array(0,c(((n-1)*(n-2)/2),2))
 c1 <- 0
 for(t in 1:(n-2))
 {
  for(s in (t+1):(n-1))
  {
   c1 <- c1+1
   COMP[c1,1] <- t
   COMP[c1,2] <- s
  }
 }

 #################################################
 if(method == "both")
 {
  mu_org <- crossing2_reverse(TimeSeries,COMP,delta)
  mu_org_weight <- crossing2_weight_reverse(TimeSeries,COMP,delta) 
  return(list(basic.measure = mu_org,weighted.measure = mu_org_weight))
 }
 if(method == "orginal")
 {
  mu_org <- crossing2_reverse(TimeSeries,COMP,delta)
  return(basic.measure = mu_org)
 } 
 if(method != "both" && method != "original")
 {
  mu_org_weight <- crossing2_weight_reverse(TimeSeries,COMP,delta,method) 
  return(weighted.measure = mu_org_weight)
 } 
}

#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################


sign_calc <- function(timeseries)
{
 n <- length(timeseries)
 return(sum(timeseries[2:n]-timeseries[1:(n-1)]))
}
sign_order <- function(ord,timeseries)
{
 return(apply(timeseries[,ord],1,sign_calc))
}
sign_allpairs <- function(timeseries)
{
 ord <- t(apply(timeseries,1,order,decreasing="FALSE"))
 return(t(apply(ord,1,sign_order,timeseries)))
}

IOTAsigned <- function(TimeSeries,method="sqrt")
{
 SAP <- sign(sign_allpairs(TimeSeries))
 IOTA_reg <- SAP*IOTA(TimeSeries,method)
}





####################################
####################################
####################################
######### to run the program #######
####################################
####################################
####################################

#########################################################
############## load or define time series ###############
#########################################################

#x7 <- x6 <- x5 <- x4 <- x3 <- x2 <- x1 <- runif(10,0,1)
#TS0 <- rbind(x1,x2,x3,x4,x5,x6,x7)

#TS0 <- matrix(runif(70,0,1),7,10)

#########################################################
#########################################################
#########################################################

#TimeSeries <- (TS0-apply(TS0,1,min,na.rm=TRUE))/apply((TS0-apply(TS0,1,min,na.rm=TRUE)),1,max,na.rm=TRUE)

# to calculate IOTA with squarded slope weight run
# I <- IOTA(TimeSeries,method="weighted")
# where TimeSeries is an array of size m*n, with m number of variables and n number of time points

