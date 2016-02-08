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
initial <- as.matrix(read.table("/home/koseska/Files/EliZamir/FA-Cell1/Data13062013/cell8/FAtracesRAW.txt"))

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

cl1 <- c(1,2,5,11,13,18,19,24,28,30,31)
cl2 <- c(33:39)
cl3 <- c(3,7,8,9,10,14,16,17,21,22,23,25,26,27,29,32)
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
#### Calculate symbolic sequences measures #############
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
    
    nn <- 3
    SVEC <- symbolvector(EXRATE,N_EXP,nn)
    A <- SVEC$A
    A2 <- SVEC$A2
    l_pattern <- SVEC$l_pattern
    ### order pattern similarity ###
    P1 <- array(0,c((a1-1),(a1-1)))
    P2 <- array(0,c((a1-1),(a1-1)))
    P3 <- array(0,c((a1-1),(a1-1)))
    for(i in 1:(a1-1)){
      if((i+1) <= (a1-1)){
	for(j in (i+1):(a1-1)){
	  p1 <- 0
	  p2 <- 0
	  for(pl in 1:l_pattern){
	    p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
	    p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i])) 
	  }
	  P1[i,j] <- p1
	  P1[j,i] <- p1
	  P2[i,j] <- p2
	  P2[j,i] <- p2
	  P3[i,j] <- max(c(P1[i,j],P2[i,j]))
	  P3[j,i] <- max(c(P1[j,i],P2[j,i]))
	}
      }
     }
    P3_norm <- P3 

    for(i in 1:(a1-1)){
      P3_norm[,i] <- P3[,i]/sum(P3[,i])
    }
    P3_norm <- (P3_norm/max(P3_norm))
    ### order pattern mi ###
    MI_A <- mutinformation(discretize(A,disc="equalwidth"),method="mm")
    ### order pattern mi normalized ###
    MI_A_norm <- MI_A
    
    for(i in 1:(a1-1)){
      MI_A_norm[,i] <- MI_A[,i]/sum(MI_A[,i])
    }
    MI_A_norm <- (MI_A_norm/max(MI_A_norm))
    ### order pattern similarity+mi ###
    C_A <- ((P3+MI_A)/2)
    ### order pattern similarity+mi normalized ###
    C_A_norm <- C_A 
    for(i in 1:(a1-1)){
      C_A_norm[,i] <- C_A[,i]/sum(C_A[,i])
    }
    C_A_norm <- (C_A_norm/max(C_A_norm))
    ### order pattern similarity+mi normalized reduced ###
    C_A_norm2 <- C_A_norm/sqrt(2) 
    results.name <- sprintf("ResultsSymbols%s-fa%s.txt",counter1,cl[[counter1]][[counter2]])
    sink(results.name)
    print("SS-Similarity")
    print(P3)
    print("SS-Similarity-Normalized")
    print(P3_norm)
    print("MI")
    print(MI_A)
    print("MI-norm")
    print(MI_A_norm)
    print("SS-Similarity-MI")
    print(C_A)
    print("SS-Similarity-MI-norm")
    print(C_A_norm)
    print("SS-Similarity-MI-norm-reduced")
    print(C_A_norm2)
    sink()
        }
}








































