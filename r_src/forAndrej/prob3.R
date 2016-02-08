library("infotheo")

prob_x <- function(X)
{
 #dig <- ceiling(log10(length(X)))
 #X <- round(X,dig)
 dyn.load("prob_x.so")
 gl <- array(0,length(X))
 px1 <- array(0,length(X))
 px2 <- array(0,length(X)) 
 ii <- 0
 re <- .C("prob_x",as.double(X), as.integer(length(X)), as.double(px1), as.double(px2), as.integer(ii), as.integer(gl))
 px <- array(0,c(re[[5]],2))
 px[,1] <- re[[3]][1:re[[5]]] 
 px[,2] <- re[[4]][1:re[[5]]] 
 return(px)
}
prob_xy <- function(X,Y)
{
 #dig <- ceiling(log10(length(X)))
 #X <- round(X,dig)
 #Y <- round(Y,dig)
 dyn.load("prob_xy.so")
 gl <- array(0,length(X))
 px1 <- array(0,length(X))
 px2 <- array(0,length(X)) 
 px3 <- array(0,length(X)) 
 ii <- 0
 re <- .C("prob_xy",as.double(X), as.double(Y), as.integer(length(X)), as.double(px1), as.double(px2), as.double(px3), as.integer(ii), as.integer(gl))
 px <- array(0,c(re[[7]],3))
 px[,1] <- re[[4]][1:re[[7]]] 
 px[,2] <- re[[5]][1:re[[7]]] 
 px[,3] <- re[[6]][1:re[[7]]] 
 return(px)
}
prob_xyz <- function(X,Y,Z)
{
 #dig <- ceiling(log10(length(X)))
 #X <- round(X,dig)
 #Y <- round(Y,dig)
 #Z <- round(Z,dig)
 dyn.load("prob_xyz.so")
 gl <- array(0,length(X))
 px1 <- array(0,length(X))
 px2 <- array(0,length(X))
 px3 <- array(0,length(X))
 px4 <- array(0,length(X)) 
 ii <- 0
 re <- .C("prob_xyz",as.double(X), as.double(Y), as.double(Z), as.integer(length(X)), as.double(px1), as.double(px2), as.double(px3), as.double(px4), as.integer(ii), as.integer(gl))
 px <- array(0,c(re[[9]],4))
 px[,1] <- re[[5]][1:re[[9]]] 
 px[,2] <- re[[6]][1:re[[9]]] 
 px[,3] <- re[[7]][1:re[[9]]] 
 px[,4] <- re[[8]][1:re[[9]]]
 return(px)
}

expect_value <- function(X)
{
 px <- prob_x(X)
 M1 <- 0
 for(i in 1:length(px[,1]))
 {
  M1 <- M1 + px[i,1]*px[i,2]
 }
 return(M1)
}

#expect_value_cond <- function(X,Y)
#{
# px_test <- prob_x(X)
# py_test <- prob_x(Y)
# pxy_test <- prob_xy(X,Y)
# M12 <- 0
# for(n in 1:length(pxy_test[,1]))                         # all possible combinations
# {
#  pxy_s <- 0
#  py_s <- 0
#  a[1:length(pxy_test[,1])] <- pxy_test[n,1]
#  b <- pxy_test[,2]
#  for(j in 1:length(pxy_test[,1]))
#  {
#   #for(i in 1:length(pxy_test[,1]))
#   #{
#   # if(pxy_test[i,1] == a[j] & pxy_test[i,2] == b[j])  pxy_s <- pxy_s + pxy_test[i,3] 
#   #}
#   #for(i in 1:length(py_test[,1]))
#   #{
#   # if(py_test[i,1] == b[j])  py_s <- py_s + py_test[i,2]
#   #}
#   pxy_s <- pxy_s + sum(pxy_test[which(pxy_test[,1] == a[j] & pxy_test[,2] == b[j]),3])
#   py_s <- py_s + sum(py_test[which(py_test[,1] == b[j]),2])
#  }
#  M12 <- M12 + pxy_test[n,1] * (pxy_s/py_s)
# }
# return(M12)
#} 
 
expect_value_cond <- function(X,Y)
{
 px_test <- prob_x(X)
 py_test <- prob_x(Y)
 pxy_test <- prob_xy(X,Y)
 px1 <- px_test[,1]
 px2 <- px_test[,2]
 py1 <- py_test[,1]
 py2 <- py_test[,2]
 pxy1 <- pxy_test[,1]
 pxy2 <- pxy_test[,2]
 pxy3 <- pxy_test[,3]
 M12 <- 0
 dyn.load("exp_c.so")
 ec <- .C("exp_c",as.double(py1), as.double(py2), as.double(pxy1), as.double(pxy2), as.double(pxy3), as.integer(length(py1)), as.integer(length(pxy1)), as.double(M12))
 return(ec[[8]])
}


residual <- function(X,Y) #residum of x from linear regression of x with y
{
 Rx <- X - expect_value(X) - (((expect_value((X-expect_value(X))*(Y-expect_value(Y))))/(expect_value((Y-expect_value(Y))*(Y-expect_value(Y))))) * (Y - expect_value(Y)))
 return(Rx)
}

################################
################################
################################

covar <- function(X,Y)
{
 C12 <- expect_value((X-expect_value(X))*(Y-expect_value(Y)))
 return(C12)
}

correl <- function(X,Y)
{
 C12_n <- covar(X,Y)/sqrt(covar(X,X)*covar(Y,Y))
 return(C12_n)
}

covar_cond <- function(X,Y,Z)
{
 C12_3 <- expect_value_cond(((X-expect_value_cond(X,Z))*(Y-expect_value_cond(Y,Z))),Z)
 return(C12_3)
}

correl_cond <- function(X,Y,Z)
{
 C12_3_c <- covar_cond(X,Y,Z)/sqrt(covar_cond(X,X,Z)*covar_cond(Y,Y,Z))
 return(C12_3_c)
}

correl_part <- function(X,Y,Z)
{
 correl(residual(X,Z),residual(Y,Z))
}

#######################################
#######################################
#######################################




mutualinformation <- function(X,Y,methode)
{
 Xd <- unlist(discretize(X,disc="equalwidth"))
 Yd <- unlist(discretize(Y,disc="equalwidth"))
 XYd <- array(0,c(length(X),2))
 XYd[,1] <- Xd
 XYd[,2] <- Yd

 I <- entropy(Xd,method=methode) + entropy(Yd,method=methode) - entropy(XYd,method=methode)
 return(I)
}

mutualinformation_cond <- function(X,Y,Z,methode)
{
 Xd <- unlist(discretize(X,disc="equalwidth"))
 Yd <- unlist(discretize(Y,disc="equalwidth"))
 Zd <- unlist(discretize(Z,disc="equalwidth"))
 XZd <- array(0,c(length(X),2))
 XZd[,1] <- Xd
 XZd[,2] <- Zd
 YZd <- array(0,c(length(X),2))
 YZd[,1] <- Yd
 YZd[,2] <- Zd
 XYZd <- array(0,c(length(X),3))
 XYZd[,1] <- Xd
 XYZd[,2] <- Yd
 XYZd[,3] <- Zd
 Ic <- -entropy(XYZd,method=methode) + entropy(XZd,method=methode) + entropy(YZd,method=methode) - entropy(Zd,method=methode)
 return(Ic)
}

mutualinformation_part <- function(X,Y,Z,methode)
{
 mutualinformation(residual(X,Z),residual(Y,Z),methode)
}
