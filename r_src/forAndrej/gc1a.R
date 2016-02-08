library(vars)
granger <- function(d, L){
#d is a bivariate time-series:  regress d[,k] on L lags of d[,1] and d[,2].
  names.d <- dimnames(d)[[2]]
  D <- d
  for(i in 1:L){
    D <- ts.intersect(D, lag(d,  - i))
#don't you find -1 bizarre here and it annoying to need the loop
  }
  dimnames(D)[[2]] <- paste(rep(names.d, L + 1), "_", rep(0:L, times = rep(2, L + 1)), sep = "")
  y  <- D[, 1] # ==========>>>>>>>> first ts
  n  <- length(y)
  x1 <- D[,  - (1:2)] # =========>>>>>>> first and second ts shifted
  x0 <- x1[, ((1:L) * 2) - 1] # =========>>>>>>> first ts shifted
  z1 <- lm(y ~ x1)
  z0 <- lm(y ~ x0)
#S1 <- sum(z1$resid^2) # =========>>>>>>> prediction error taking both ts
#S0 <- sum(z0$resid^2) # =========>>>>>>> prediction error taking only first ts
#ftest <- ((S0 - S1)/L)/(S1/(n - 2 * L - 1)) # ==========>>>>>>>> p0 = L-1, p1 = 2L-1
#list(ftest = ftest, p.val = 1 - pf(ftest, L, n - 2 * L - 1), R2 = summary(z1)$r.squared)
  F <- log(var(z0$resid)/var(z1$resid))  # 2 gc 1 if F > 0
} 

gc <- function(x,y){
# x -> y
  A1 <- ts(x)
  A2 <- ts(y)
  LAG <- floor(mean(c(VARselect(A1)$selection[[1]],VARselect(A2)$selection[[1]]),na.rm = TRUE))
  F <- granger(cbind(A2,A1), L=LAG)
}

granger_cond <- function(d, L){
#d is a trivariate time-series:  regress d[,k] on L lags of d[,1] and d[,2]and d[,3].
  names.d <- dimnames(d)[[2]]
  D <- d
  for(i in 1:L){
    D <- ts.intersect(D, lag(d,  - i))
#don't you find -1 bizarre here and it annoying to need the loop
  }
  dimnames(D)[[2]] <- paste(rep(names.d, L + 1), "_", rep(0:L, times = rep(3, L + 1)), sep = "")
  y  <- D[, 1] # ==========>>>>>>>> first ts
  n  <- length(y)
  x1 <- D[,  - (1:3)] # =========>>>>>>> first and second and third ts shifted
  x0 <- x1[, c(((1:L) * 3) - 2,((1:L) * 3))] # =========>>>>>>> first and thrid ts shifted
  z1 <- lm(y ~ x1)
  z0 <- lm(y ~ x0)
#S1 <- sum(z1$resid^2) # =========>>>>>>> prediction error taking all three ts
#S0 <- sum(z0$resid^2) # =========>>>>>>> prediction error taking only first and third ts
#ftest <- ((S0 - S1)/L)/(S1/(n - 2 * L - 1)) # ==========>>>>>>>> p0 = L-1, p1 = 2L-1
#list(ftest = ftest, p.val = 1 - pf(ftest, L, n - 2 * L - 1), R2 = summary(z1)$r.squared)
  F <- log(var(z0$resid)/var(z1$resid)) # 2 gc 1 if F > 0
} 

gc_c <- function(x,y,z){
# x -> y | z
  A1 <- ts(x)
  A2 <- ts(y)
  A3 <- ts(z)
  LAG <- floor(mean(c(VARselect(A1)$selection[[1]],VARselect(A2)$selection[[1]],VARselect(A3)$selection[[1]]),na.rm = TRUE))
  F <- granger_cond(cbind(A2,A1,A3), L=LAG)
}

granger_part <- function(d, L) 
{
#d is a trivariate time-series:  regress d[,k] on L lags of d[,1] and d[,2]and d[,3].
names.d <- dimnames(d)[[2]]
D <- d
for(i in 1:L)
{
D <- ts.intersect(D, lag(d,  - i))
#don't you find -1 bizarre here and it annoying to need the loop
}
dimnames(D)[[2]] <- paste(rep(names.d, L + 1), "_", rep(0:L, times = rep(3, L + 1)), sep = "")
y  <- D[, 1] # ==========>>>>>>>> first ts
n  <- length(y)
x1 <- D[,  - (1:3)] # =========>>>>>>> first and second and third ts shifted
x0 <- x1[, c(((1:L) * 3) - 2,((1:L) * 3))] # =========>>>>>>> first and thrid ts shifted
z1 <- lm(y ~ x1)
z0 <- lm(y ~ x0)
y_2  <- D[, 3] # ==========>>>>>>>> third ts
n  <- length(y)
x1_2 <- D[,  - (1:3)] # =========>>>>>>> first and second and third ts shifted
x0_2 <- x1_2[, c(((1:L) * 3) - 2,((1:L) * 3))] # =========>>>>>>> first and thrid ts shifted
z1_2 <- lm(y_2 ~ x1_2)
z0_2 <- lm(y_2 ~ x0_2)

F <- log((var(z0$resid) - (cov(z0$resid,z0_2$resid)*((var(z0_2$resid))^(-1))*cov(z0_2$resid,z0$resid)))/(var(z1$resid) - (cov(z1$resid,z1_2$resid)*((var(z1_2$resid))^(-1))*cov(z1_2$resid,z1$resid))))  # 2 gc 1 if F > 0
}

gc_p <- function(x,y,z) # x -> y | z
{
 A1 <- ts(x)
 A2 <- ts(y)
 A3 <- ts(z)
 LAG <- floor(mean(c(VARselect(A1)$selection[[1]],VARselect(A2)$selection[[1]],VARselect(A3)$selection[[1]]),na.rm = TRUE))
 F <- granger_part(cbind(A2,A1,A3), L=LAG)
}