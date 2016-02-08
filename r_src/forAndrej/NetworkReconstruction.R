library("Rcpp")
library("apcluster")
library("stringr")
library("igraph")
library("minet")
library("lmtest")
library("stats")
library("MSBVAR")
library("vars")
library("entropy")
library("dtw")

source("sub1.R")
source("mutual_count.R")

source("prob3.R")
source("gc1.R")
source("/home/koseska/Files/EliZamir/FA-Cell1/Data13062013/cell8/symbolvector.R")
source("/home/koseska/Files/EliZamir/FA-Cell1/Data13062013/cell8/plotting.R")


####### read data ###############
initial <- as.matrix(read.table("/home/koseska/Files/EliZamir/FA-Cell1/Data13062013/cell9/FAtracesRAW.txt"))

##### set Na to 0 ##########
initial[is.na(initial)] <- 0

####### number of FAs ##########
fa <- ncol(initial)/4

##### extract FA intensities #####
data <- vector(mode="list", length=fa)
data[[1]] <- initial[,1:4]
k <- 1
for (i in 2:fa){
  data[[i]] <- initial[,(i*2+k):(i*2+k+3)]
  k=k+2
}

####### which are the stationary FAs ###############

stable <- array(1,c(fa))
for (i in 1:fa){
  counter <- 0
  for(j in 1:47){
    for(k in 1:4){
      if (data[[i]][j,k]==0) counter=counter+1
    }
  }
  if (counter>1) stable[i]=0
}

##### normalize as (I-min)/(max-min) within FA ###########

minval <- vector(mode="list",length=fa)
maxval <- vector(mode="list",length=fa)
res <- vector(mode="list",length=fa)

for (i in 1:fa){
  minval[[i]] <- array(0,c(4))
  maxval[[i]] <- array(0,c(4))
  res[[i]] <- array(0,c(nrow(initial),4))
  for(j in 1:4){
    minval[[i]][j] <- min(data[[i]][,j])
    maxval[[i]][j] <- max(data[[i]][,j])
    res[[i]][,j] <- (data[[i]][,j]-minval[[i]][j])/(maxval[[i]][j]-minval[[i]][j])
  }
}

######## Prepare the data set for clustering ###########
dataset <- array(0,c(fa,length(data[[1]])))
for(i in 1:fa){
  dataset[i,] <- c(res[[i]])
}

####### Affinity propagation clustering ###########
s1 <- negDistMat(dataset,r=2)
apres1a <- apcluster(s1,q=0.6)

#####TO DO:combine stationary with ap results and take sets of FA for reconstruction automatically ###############

cl1 <- c(2,4,5,6,8,10,13,16,17,21,22,25,32,33,38)
cl2 <- c(12,15,20,23,26,27,40,41,46)
cl3 <- c(1,3,7,9,11,14,18,19,24,28,29,30,34,35,36,37)
cl <- vector(mode="list",length=3)
cl[[1]] <- cl1###stationary
cl[[2]] <- cl2###appearing
cl[[3]] <- cl3##disapearing
clist <- vector(mode="list",length=3)

###clist contains the 3 lists, and each one is a list of all FAs in the group #########
for(i in 1:3){
  clist[[i]] <- vector(mode="list",length(cl[[i]]))
    for(j in 1:length(cl[[i]])){
      clist[[i]][[j]] <- res[[cl[[i]][j]]] 
    }
}

#### Calculate all BMC Bioinformatics measures #############
#results <- vector(mode="list",length=3)
i<-1
j<-1
k<-1
for(counter1 in 1:3){
  #results[[i]] <- vector(mode="list",length(cl[[i]]))
  for(counter2 in 1:length(cl[[counter1]])){
   EXRATE <- t(clist[[counter1]][[counter2]])
   N_EXP <- length(EXRATE[1,]) 
   N_VAR <- length(EXRATE[,1])
   a1 <-  N_VAR +1
    
   methode <- "MM"
   gc_simp <- array(0,c((a1-1),(a1-1)))
   mi_simp <- array(0,c((a1-1),(a1-1)))
   cor_simp <- array(0,c((a1-1),(a1-1)))
   gc_part <- array(0,c((a1-1),(a1-1),(a1-1)))
   mi_part <- array(0,c((a1-1),(a1-1),(a1-1)))
   cor_part <- array(0,c((a1-1),(a1-1),(a1-1)))
   gc_cond <- array(0,c((a1-1),(a1-1),(a1-1)))
   mi_cond <- array(0,c((a1-1),(a1-1),(a1-1)))
   cor_cond <- array(0,c((a1-1),(a1-1),(a1-1)))
   cor_spearman <- array(0,c((a1-1),(a1-1)))
   cor_kendal <- array(0,c((a1-1),(a1-1)))
   for(i in 1:(a1-1)){
    if(i < (a1-1)){
      for(j in (i+1):(a1-1)){
	for(k in 1:(a1-1)){
	  print(cat("i: ",i," | j: ",j," | k: ",k,"\n"))     
	  T1 <- EXRATE[i,]
	  T2 <- EXRATE[j,]
	  T3 <- EXRATE[k,]
	  if(k == 1){
	    ## granger causality
	    print("simple")
	    l1 <- VARselect(T1,lag.max = 6)$selection[[1]]
	    l2 <- VARselect(T2,lag.max = 6)$selection[[1]]
	    if(is.finite(l1) == 0)  l1 <- NA
	    if(is.finite(l2) == 0)  l2 <- NA
	    LAG <- floor(mean(c(l1,l2),na.rm = TRUE))
	    if(is.na(LAG)) LAG <- 1
	    gc_simp[i,j] <- granger(cbind(T2,T1), L=LAG)
	    gc_simp[j,i] <- granger(cbind(T1,T2), L=LAG)

	    ## mutual information
	    print("simple")
	    mi_simp[i,j] <- mutualinformation(T1,T2,methode)
	    mi_simp[j,i] <- mi_simp[i,j]
 
	    ## pearson correlation
	    print("simple")
	    cor_simp[i,j] <- correl(T1,T2)
	    cor_simp[j,i] <- cor_simp[i,j]

	    ## rank correlation
	    print("spearman")
	    cor_spearman[i,j] <- cor(T1,T2,method="spearman")

	    print("kendal")
	    cor_kendal[i,j] <- cor(T1,T2,method="kendal")
	   }
	   if(k != i & k != j){ 
	    ## granger causality
	    l1 <- VARselect(T1,lag.max = 6)$selection[[1]]
	    l2 <- VARselect(T2,lag.max = 6)$selection[[1]]
	    l3 <- VARselect(T3,lag.max = 6)$selection[[1]]
	    if(is.finite(l1) == 0)  l1 <- NA
	    if(is.finite(l2) == 0)  l2 <- NA
	    if(is.finite(l3) == 0)  l3 <- NA
	    LAG <- floor(mean(c(l1,l2,l3),na.rm = TRUE))
	    if(is.na(LAG)) LAG <- 1
	    print("partial")
	    gc_part[i,j,k] <- granger_part(cbind(T2,T1,T3), L=LAG)
	    gc_part[j,i,k] <- granger_part(cbind(T1,T2,T3), L=LAG)
	    print("contitional")
	    gc_cond[i,j,k] <- granger_cond(cbind(T2,T1,T3), L=LAG)
	    gc_cond[j,i,k] <- granger_cond(cbind(T1,T2,T3), L=LAG)
     
	    ## mutual information
	    print("partial")

	    mi_part[i,j,k] <- mutualinformation_part(T1,T2,T3,methode)
	    mi_part[j,i,k] <- mi_part[i,j,k]

	    print("contitional")

	    mi_cond[i,j,k] <- mutualinformation_cond(T1,T2,T3,methode)
	    mi_cond[j,i,k] <- mi_cond[i,j,k]

	    ## pearson correlation
	    print("partial")
 
	    cor_part[i,j,k] <- correl_part(T1,T2,T3)
	    cor_part[j,i,k] <- cor_part[i,j,k]

	    print("contitional")
	    cor_cond[i,j,k] <- correl_cond(T1,T2,T3)
	    cor_cond[j,i,k] <- cor_cond[i,j,k]
	   }
	  }
	 }
	}
      }
      ## MCIR - coarse grained measures###
      DATA <- t(EXRATE)
      L <- length(DATA[,i])
      Ixy <- array(0,c((a1-1),(a1-1)))
      tau_max <- (L-1)
      for(i in 1:(a1-1)){
	for(j in 1:(a1-1)){ 
	  for(tau in -tau_max: tau_max){
	    if(tau < 0){
	      X <- DATA[(-tau+1):L,i]
	      Y <- DATA[1:(L+tau),j]
	      I <- mutualinformation(X,Y,"MM")
	     }
	     if(tau > 0){
	      X <- DATA[1:(L-tau),i]
	      Y <- DATA[(tau+1):L,j]
	      I <- mutualinformation(X,Y,"MM")
	     }
	     if(tau == 0){
	      I <- 0
	     }
	     Ixy[i,j] <- Ixy[i,j] + I
	    }
	    Ixy[i,j] <- Ixy[i,j]/(2*tau_max)
	  }
      }
   
      ### conditional CIR - coarse grained ###
      i_xy <- 0
      DATA <- t(EXRATE)
      L <- length(DATA[,1])
      i0_xy <- array(0,c((a1-1),(a1-1)))
      i_x <- array(0,c((a1-1),(a1-1)))
      tau_max <- (L-1)
      for(i in 1:(a1-1)){
	for(j in 1:(a1-1)){ 
	  for(tau in 1:tau_max){
	    if(tau > 0){
	      X <- DATA[1:(L-tau),i]
	      X_tau <- DATA[(tau+1):L,j]
	      Y <- DATA[1:(L-tau),j]
	      I0 <- mutualinformation_cond(X,X_tau,Y,"MM")
	      I <- mutualinformation(X,X_tau,"MM")
	    }
	    i0_xy[i,j] <- i0_xy[i,j] + I0
	    i_x[i,j] <- i_x[i,j] + I
	  }
	  i0_xy[i,j] <- i0_xy[i,j]/(tau_max)
	  i_x[i,j] <- i_x[i,j]/(tau_max)
	}
      }
      i_xy <- i0_xy-i_x
   
      ###### Distances #######################
      distances <- function(x,y){
      #Eucidean
      d_euc <- sqrt(sum((x-y)^2))
      #Manhattan
      d_man <- sum(abs(x-y))
      #Lm
      m <- length(x)
      d_lm <- (sum(abs(x - y)^m))^(1/m)
      return(list(euclidean = d_euc, manhattan = d_man, Lm = d_lm))
      }

      D_euc <- array(0,c((a1-1),(a1-1)))
      D_man <- array(0,c((a1-1),(a1-1)))
      D_lm <- array(0,c((a1-1),(a1-1)))
      for(i in 1:(a1-1)){
	for(j in 1:(a1-1)){
	  print(cat(i," | ",j,"\n"))
	  D <- distances(EXRATE[i,],EXRATE[j,])
	  D_euc[i,j] <- D$euclidean
	  D_man[i,j] <- D$manhattan
	  D_lm[i,j] <- D$Lm
	}
      }
      #D_dtw1 <- dtwDist(EXRATE,method="DTW",dist.method="Euclidean",step.pattern=symmetric1)
      #D_dtw2 <- dtwDist(EXRATE,method="DTW",dist.method="Euclidean",step.pattern=symmetric2)
      #D_dtw3 <- dtwDist(EXRATE,method="DTW",dist.method="Euclidean",step.pattern=asymmetric)
   
      results.name <- sprintf("Cell9Results%s-fa%s.txt",counter1,cl[[counter1]][[counter2]])
      sink(results.name)
      print("MI-SIMP")
      print(mi_simp)	
      print("MI-PART")
      print(mi_part)
      print("MI-COND")
      print(mi_cond)
      print("COR-SIMP")
      print(cor_simp)
      print("COR-PART")
      print(cor_part)
      print("COR-COND")
      print(cor_cond)
      print("COR-Spearman")
      print(cor_spearman)
      print("COR-Kendal")
      print(cor_kendal)
      print("Granger-simple")
      print(gc_simp)
      print("Granger-partial")
      print(gc_part)
      print("Granger-cond")
      print(gc_cond)
      print("Ixy")
      print(Ixy)
      print("ixy")
      print(i_xy)
      print("Euclidean")
      print(D_euc)
      print("Manhatan")
      print(D_man)
      print("lm")
      print(D_lm)
      #print("DTW1")
      #print(D_dtw1)
      #print("DTW2")
      #print(D_dtw2)
      #print("DTW3")
      #print(D_dtw3)
      sink()
   }
}








































