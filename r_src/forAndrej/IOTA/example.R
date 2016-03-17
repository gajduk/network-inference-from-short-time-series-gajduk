
# load subroutines R
source('IOTA.R')

# load or simulate time series as an array of size m*n, with m number of variables and n number of time points
TS0 <- matrix(runif(70,0,1),7,10)

# normalize time series
# to estimate the weighted IOTA the time series must have values between zero and one
# depending on the time series different normalization must be used
TimeSeries <- (TS0-apply(TS0,1,min,na.rm=TRUE))/apply((TS0-apply(TS0,1,min,na.rm=TRUE)),1,max,na.rm=TRUE)

# calculates pairwise IOTA as described in Hempel et al., PRL (2011)
# possible options for method (weighting functions) 
# 'both' (default): uniform and squared sloped 
# 'slope': slope 
# 'sqrt': squared slope 
# 'am': arithmetric mean
# 'gm': geometric mean
# 'hm': harmonic mean 
I <- IOTA(TimeSeries,method='sqrt')
# calculates pairwise IOTA based on reversed ordering as described in Hempel et al., EPJB (2013)
# options are the same as for IOTA
Ir <- IOTA_reverse(TimeSeries,method='sqrt')
# calculates signed version of pairwise IOTA to indicate up-/downregulation as described in Hempel et al., EPJB (2013)
# option 'both' does not work in this case
Is <- IOTAsigned(TimeSeries,method='sqrt')



###########################################################################
###########################################################################

# to run C subroutines files must be compiled to get a dynamic library using  "R CMD SHLIB *.c" command rather than calling C-compiler directly

source('iota_subroutines.R')

# number of realizations for significance test
rmax <- 1000
# significance level
alpha <- 0.99
# weighting: uniform (1) or squared slope (2)
w <- 2

# calculated pairwise and partial IOTA and performs a simple permutation-based significance test, only most likely connections are selected while the remaining matix entries are set to zero
I <- IOTA(TimeSeries,rmax,alpha,w)
