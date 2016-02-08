entropy_information <- function(X)
{
 L1 <- length(X)
 TPvectorX <- table(X)/length(X)                        # p(x)
 SUMvector <- TPvectorX * log10(TPvectorX)
 # sum_x (p_x*log(p_x))
 return(-sum(SUMvector))
}


joint_entropy_information <- function(X)
{
 n <- length(X[1,])  # col = variable, row = outcome
 L1 <- length(X[,1])
 TPvectorX <- array(0,L1)
 for(i in 1:L1)
 {
  cc <- array("i",((2*n)-1))
  for(m in 1:n)
  {
   cc[((2*m)-1)] <- X[i,m]
  }
  TPvectorX[i] <- paste(cc,collapse="")
 }
 TPvectorX <- table(TPvectorX)/length(TPvectorX)     # p(x1,..,xn)

 SUMvector <- TPvectorX * log10(TPvectorX)
 # sum_x (p_x*log(p_x))
 return(-sum(SUMvector))
}


mutual_information <- function(X,Y)
{
 L1 <- length(X)
 TPvectorXY <- array(0,L1)
 for(i in 1:L1)
 {
  TPvectorXY[i] <- paste(c(X[i],"i",Y[i]),collapse="")
 }
 TPvectorXY <- table(TPvectorXY)/length(TPvectorXY)     # p(x,y)
 TPvectorX <- table(X)/length(X)                        # p(x)
 TPvectorY <- table(Y)/length(Y)                        # p(y)

 SUMvector <- array(0,length(TPvectorXY))

 for(n in 1:length(TPvectorXY))                         # all possible combinations
 {
  m <- length(strsplit(names(TPvectorXY)[1],"i")[[1]])
 
  SUMvector[n] <- TPvectorXY[n] * log10(TPvectorXY[n]/( TPvectorX[unlist(strsplit(names(TPvectorXY)[n],"i"))[1]] * TPvectorY[unlist(strsplit(names(TPvectorXY)[n],"i"))[m]]))
  # sum_y sum_x (p_xy*log(p_xy/(p_x*p_y)))
 }
 return(sum(SUMvector))
}


cond_mutual_information <- function(X,Y,Z)
{
 L1 <- length(X)
 
 TPvectorXYZ <- array(0,L1)
 for(i in 1:L1)
 {
  TPvectorXYZ[i] <- paste(c(X[i],"i",Y[i],"i",Z[i]),collapse="")
 }
 TPvectorXYZ <- table(TPvectorXYZ)/length(TPvectorXYZ)  # p(x,y,z)
 #TPvectorX <- table(X)/length(X)                       # p(x)
 #TPvectorY <- table(Y)/length(Y)                       # p(y)
 TPvectorZ <- table(Z)/length(Z)                        # p(z)
 
 #TPvectorXY <- array(0,L1)
 #for(i in 1:L1)
 #{
 # TPvectorXY[i] <- paste(c(X[i],"i",Y[i]),collapse="")
 #}
 #TPvectorXY <- table(TPvectorXY)/length(TPvectorXY)    # p(x,y)

 TPvectorXZ <- array(0,L1)
 for(i in 1:L1)
 {
  TPvectorXZ[i] <- paste(c(X[i],"i",Z[i]),collapse="")
 }
 TPvectorXZ <- table(TPvectorXZ)/length(TPvectorXZ)     # p(x,z)

 TPvectorYZ <- array(0,L1)
 for(i in 1:L1)
 {
  TPvectorYZ[i] <- paste(c(Y[i],"i",Z[i]),collapse="")
 }
 TPvectorYZ <- table(TPvectorYZ)/length(TPvectorYZ)     # p(y,z)


 SUMvector <- array(0,length(TPvectorXYZ))
 for(n in 1:length(TPvectorXYZ))                        # all possible combinations
 {
  m <- strsplit(names(TPvectorXYZ)[n],"i")
 
  SUMvector[n] <- TPvectorXYZ[n] * log10(((TPvectorXYZ[n]) * (TPvectorZ[(unlist(m))[3]])) / ((TPvectorXZ[paste((unlist(m))[1],"i",(unlist(m))[3],sep="",collapse="")]) * (TPvectorYZ[paste((unlist(m))[2],"i",(unlist(m))[3],sep="",collapse="")])))
  # sum_z sum_y sum_x p(x,y,z)*log((p(x,y,z) * p(z)) / (p(x,y) * p(y,z)))
 }

 return(sum(SUMvector))
}