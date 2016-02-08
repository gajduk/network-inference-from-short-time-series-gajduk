##########################################################################
##########################################################################
##########################################################################
################### begin subroutines ####################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
crops_all_tp <- function(EXRATE,n,m,rmax,alpha,w,Ta)
{
 print("calculate IOTA")
 ###########################################################
 ##### use C program to calculate crossing points ##########
 ###########################################################
 # R CMD SHLIB crossing.c
 dyn.load("crossing.so")
 ###########################################################
 M <- as.numeric(t(EXRATE))  
 MU <- array(0,(m*m))
 MUsig <- array(0,(m*m))
 MUsigd <- array(0,(m*m))
 ###########################################################
 re <- .C("crossing", as.double(M), as.integer(n), as.integer(m), as.double(MU), as.double(MUsig), as.integer(rmax), as.double(MUsigd), as.integer(w))
 ########################################################### 
 MU_m <- array(0,c(m,m))
 MUsig_m <- array(0,c(m,m))
 MUsigd_m <- array(0,c(m,m))

 MU_m <- matrix(re[[4]],m,m)
 MUsig_m <- matrix(re[[5]],m,m)
 MUsigd_m <- matrix(re[[7]],m,m)

 ###########################################################
 MU_m_2 <- MU_m
 MU_m_2[which( ((MU_m-t(MU_m)) < 0 &  MUsigd_m  <= alpha) | (MUsig_m > alpha) )] <- 0
 ###########################################################
 ###########################################################
 ##### use C program to reveal autoregulation ##############
 ###########################################################
 # R CMD SHLIB autoreg.c
 #dyn.load("autoreg.so")
 ###########################################################
 MU_m_3 <- MU_m_2
 MU_m_3[(Ta-1)*length(MU_m_3[,1])+Ta] <- 1
 ###########################################################

 return(MU_m_3)
}
crops_part_tp <- function(EXRATE,n,m,rmax,alpha,w,Ta,e)
{
 print("calculate crops")
 ###########################################################
 ##### use C program to calculate crossing points ##########
 ###########################################################
 # R CMD SHLIB crossing.c
 dyn.load("crossing_narm.so")
 ###########################################################
 M <- as.numeric(t(EXRATE))  
 MU <- array(0,(m*m))
 MUsig <- array(0,(m*m))
 MUsigd <- array(0,(m*m))
 ###########################################################
 re <- .C("crossing_narm", as.double(M), as.integer(n), as.integer(m), as.double(MU), as.double(MUsig), as.integer(rmax), as.double(MUsigd), as.integer(w), as.integer(e))
 ########################################################### 
 MU_m <- array(0,c(m,m))
 MUsig_m <- array(0,c(m,m))
 MUsigd_m <- array(0,c(m,m))

 MU_m <- matrix(re[[4]],m,m)
 MUsig_m <- matrix(re[[5]],m,m)
 MUsigd_m <- matrix(re[[7]],m,m)

 ###########################################################
 MU_m_2 <- MU_m
 MU_m_2[which( ((MU_m-t(MU_m)) < 0 &  MUsigd_m  <= alpha) | (MUsig_m > alpha) )] <- 0
 ###########################################################
 ###########################################################
 ##### use C program to reveal autoregulation ##############
 ###########################################################
 # R CMD SHLIB autoreg.c
 #dyn.load("autoreg.so")
 ###########################################################
 MU_m_3 <- MU_m_2
 MU_m_3[(Ta-1)*length(MU_m_3[,1])+Ta] <- 1
 ###########################################################

 return(MU_m_3)
}
#####################################################################
findpairs <- function(ts,TimeSeries)
{
 return(which(is.finite(rowSums(ts[which(is.finite(ts))]+TimeSeries[,which(is.finite(ts))]))))
}
#####################################################################
IOTA <- function(TimeSeries,rmax,alpha,w)
{
 ## TimeSeries[number_of_genes,number_of_timepoints]  ... time series of gene expression data
 ## rmax ... number of random permutations for significance test
 ## alpha ... significance level
 ## w ... choice of weighting for inner composition alignment (1 ... slope, 2 ... squared slope)

 ## define squared matrix (dimension: number_of_genes) containing chance of regulatory linkage infered from similarity of time series

print("pairswise, full length")

   CROPS <- array(NA,c(length(TimeSeries[,1]),length(TimeSeries[,1])))

 ## find all genes with complete time series

   finite <- array(0,length(TimeSeries[,1]))
   finite[which(!is.finite(rowSums(TimeSeries)))] <- 1
   s <- which(finite==0)

 ## choose only genes with complete time series 

   ## time series of full length
   TS0 <- TimeSeries[s,]
   ## number of time points
   n <- length(TS0[1,])
   ## number of genes
   m <- length(TS0[,1])
   ## set diagonal to one to allow for autoregulation
   Ta <- seq(1,m)
   
 ## calculate inner composition alignment among genes with complete time series
 
   CROPS_tmp <- crops_all_tp(TimeSeries[s,],n,m,rmax,alpha,2,Ta) 

   ## insert chance of regulatory linkage infered from similarity of time series for gene pairs with complete time series both
   
   CROPS[s,s] <- CROPS_tmp

print("pairswise, incomplete ts")

 ## find all genes with incomplete time series

   e <- which(!finite==0)

 ## find gene pairs with incomplete time series but matching time points

   list1 <- apply(TimeSeries[e,],1,findpairs,TimeSeries)

 ## loop over all genes with incomplete time series

   if(length(e)>0)
   {
    for(pp in 1:length(e))
    {

     ## gene with incomplete time series and all genes with matching time series

     s1 <- as.numeric(unlist(list1[pp]))

     if(is.finite(s1))
     {
 
      ## choose matching time points

      if(length(s1)>1)
      {
       tp1 <- which(is.finite(colSums(TimeSeries[s1,])))
      }
      if(length(s1)==1)
      {
       tp1 <- which(is.finite(TimeSeries[s1,]))
      }

      ## calculate inner composition alignment among genes based on series of matching time points

      CROPS_tmp <- crops_part_tp(TimeSeries[s1,tp1],length(tp1),length(s1),rmax,alpha,2,Ta,which(s1==e[pp]))  

      ## insert chance of regulatory linkage infered from similarity of time series for gene pairs with incomplete time series 

      CROPS[e[pp],s1] <- CROPS_tmp[which(s1==e[pp]),] 
      CROPS[s1,e[pp]] <- CROPS_tmp[,which(s1==e[pp])] 

     }
    }
   }
 ###########################################################
 print("calculate IOTA partial")
 dyn.load("superfluous.so")
 ###########################################################
 ###########################################################
 ##### use C program to identify superfluous ###############
 ###########################################################
 ###########################################################

 ## maximal number to forward to C program

   max_triple <-  500000

 ## define squared matrix (dimension: number_of_genes) containing chance of regulatory linkage infered from similarity of time series and superfluous removed

   CROPSsfrm <- CROPS

 ## choose all gene pairs with complete time series both and select the related submatrix from that infered by inner composition alignment

   N0 <- CROPS[s,s]
   
 ## find potentially linked gene triplets

   N0[which(CROPS[s,s]>0)] <- 1
   ab <- (N0==1)
   i2 <- ceiling(which(ab)/length(ab[,1]))
   i1 <- which(ab)-((ceiling(which(ab)/length(ab[,1]))-1)*length(ab[,1]))
   abc <- (N0[,i1]==1 & N0[,i2]==1)
   ii1 <- ceiling(which(abc)/length(abc[,1]))
   r1 <- which(abc)-((ceiling(which(abc)/length(abc[,1]))-1)*length(abc[,1]))

print("triplets, full length")

 ## choose time series of all genes with complete time series

   M <- as.numeric(t(TimeSeries[s,])) 

 ###########################################################
 ## find superfluous links
 ###########################################################

 ## define temporary squared matrix (dimension: number_of_genes_with_complete_timeseries) containing chance of regulatory linkage infered from similarity of time series and superfluous removed


 CROPSsfrm_tmp <- find_superfluous(CROPS[s,s],max_triple,r1,M,s,rmax,w,i1,i2,ii1,n,alpha)

 
 #  CROPSsfrm_tmp <- CROPS[s,s]
 
 ## define segments of the list of triplets not longer than max_triple

 #  ll <- ceiling(length(r1)/max_triple)
 #  ii_s <- 1+max_triple*(seq(0,ll-1))
 #  ii_e <- max_triple+max_triple*(seq(0,ll-1))
 #  ii_e[which(ii_e > length(r1))] <- length(r1) 
 
 ## find superfluous for all the segments of the list of triplets

 #  for(index in 1:ll)
 #  {
 #   pa <- ii_s[index]:ii_e[index] ## segment
 #   MU_tri <- array(0,length(i1[ii1[pa]]))
 #   MUsig_tri <- array(0,length(i1[ii1[pa]]))
 
    ## calculate inner composition alignment for triplets, n time points and m genes involved  

 #     re <- .C("superfluous", as.double(M), as.integer(n), as.integer(m), as.double(MU_tri), as.integer(i1[ii1[pa]]), as.integer(i2[ii1[pa]]), as.integer(r1[pa]), as.integer(length(i1[ii1[pa]])), as.integer(rmax), as.double(MUsig_tri), as.integer(w))

    ## find significant values

 #     ii <- which(re[[10]] > alpha)

    ## remove superfluous in temporary matrix

 #     CROPSsfrm_tmp[(i2[ii1[pa]][ii]-1)*m+i1[ii1[pa]][ii]] <- 0
 #  }

  ## remove superfluous 

    CROPSsfrm[s,s][which(CROPSsfrm_tmp==0)] <- 0

  ###########################################################
  ###########################################################
  ###########################################################

  ## find all potentially linked gene triplets 

   N <- CROPS
   N[which(CROPS>0)] <- 1
   ab <- (N==1)
   i2 <- ceiling(which(ab)/length(ab[,1]))
   i1 <- which(ab)-((ceiling(which(ab)/length(ab[,1]))-1)*length(ab[,1]))

   abc <- (N[,i1]==1 & N[,i2]==1)
   ii1 <- ceiling(which(abc)/length(abc[,1]))
   r1 <- which(abc)-((ceiling(which(abc)/length(abc[,1]))-1)*length(abc[,1]))
 
print("triplets, incomplete ts")

  ## identify potentially linked gene triplets where at least one gene has an incomplete time series
  
    j1 <- which(!is.finite(rowSums(TimeSeries[i1[ii1],])) | !is.finite(rowSums(TimeSeries[i2[ii1],])) | !is.finite(rowSums(TimeSeries[r1,])))

    i1n <- i1[ii1][j1]
    i2n <- i2[ii1][j1]
    i3n <- r1[j1]

    ts <- is.finite(TimeSeries[i1n,]+TimeSeries[i2n,]+TimeSeries[i3n,])

  if(length(ts)>0)  
  {
  ## define matrix of possible NA position in the time series, two time points NA at most 

    nnanz <- length(TimeSeries[1,])*(length(TimeSeries[1,])-1)/2
    napos <- array(1,c((length(TimeSeries[1,])+nnanz),length(TimeSeries[1,])))
    diag(napos[seq(1,length(TimeSeries[1,])),seq(1,length(TimeSeries[1,]))]) <- 0
    count <- length(TimeSeries[1,])
    for(count1 in 1:(length(TimeSeries[1,])-1))
    {
     for(count2 in (count1+1):length(TimeSeries[1,]))
     {
      count <- count + 1
      napos[count,c(count1,count2)] <- 0
     }
    }

  ## loop over possible NA position

    for(count in 1:length(napos[,1]))
    {
     
     ## find all triplets from the list with incomplete time series, which have NA at the defined position

       j2 <- which(colSums(napos[count,]!=t(ts)) == 0)

     if(length(j2)>0)
     {
      ## infer genes from the remaining list of triplets
 
        TSf <- array(0,length(TimeSeries[,1]))
        TSf[i1n[j2]] <- 1
        TSf[i2n[j2]] <- 1
        TSf[i3n[j2]] <- 1
        se <- which(TSf==1)

      ## choose gene pairs with incomplete time series, but matching NA positions and select the related submatrix from that infered by inner composition alignment
      
        N1 <- CROPS[se,se]
      
      ## find potentially linked gene triplets
      
        N1[which(CROPS[se,se]>0)] <- 1
        ab <- (N1==1)
        i2 <- ceiling(which(ab)/length(ab[,1]))
        i1 <- which(ab)-((ceiling(which(ab)/length(ab[,1]))-1)*length(ab[,1]))
        abc <- (N1[,i1]==1 & N1[,i2]==1)
        ii1 <- ceiling(which(abc)/length(abc[,1]))
        r1 <- which(abc)-((ceiling(which(abc)/length(abc[,1]))-1)*length(abc[,1]))
      
      ## choose time series of all genes with incomplete time series, but matching NA positions and remove time points containing NA
      
        M <- as.numeric(t(TimeSeries[se,which(napos[count,]==1)])) 
      
      ###########################################################
      ## find superfluous links 
      ###########################################################

      ## define temporary squared matrix (dimension: number_of_genes_with_incomplete_timeseries_matchingNA) containing chance of regulatory linkage infered from similarity of time series and superfluous removed
  
    CROPSsfrm_tmp <- find_superfluous(CROPS[se,se],max_triple,r1,M,se,rmax,w,i1,i2,ii1,length(which(napos[count,]==1)),alpha)
    
 #       CROPSsfrm_tmp <- CROPS[se,se]
      
      ## define segments of the list of triplets not longer than max_triple
      
 #      ll <- ceiling(length(r1)/max_triple)
 #      ii_s <- 1+max_triple*(seq(0,ll-1))
 #      ii_e <- max_triple+max_triple*(seq(0,ll-1))
 #      ii_e[which(ii_e > length(r1))] <- length(r1) 
      
      ## find superfluous for all the segments of the list of triplets

 #     for(index in 1:ll)
 #     {
 #      pa <- ii_s[index]:ii_e[index]
 #      MU_tri <- array(0,length(i1[ii1[pa]]))
 #      MUsig_tri <- array(0,length(i1[ii1[pa]]))

       ## calculate inner composition alignment for triplets, length(which(napos[count,]==1)) time points and length(se) genes involved  
       
 #        re <- .C("superfluous", as.double(M), as.integer(length(which(napos[count,]==1))), as.integer(length(se)), as.double(MU_tri), as.integer(i1[ii1[pa]]), as.integer(i2[ii1[pa]]), as.integer(r1[pa]), as.integer(length(i1[ii1[pa]])), as.integer(rmax), as.double(MUsig_tri), as.integer(w))
       
       ## find significant values
       
 #        ii <- which(re[[10]] > alpha)
       
       ## remove superfluous in temporary matrix
       
 #        CROPSsfrm_tmp[(i2[ii1[pa]][ii]-1)*length(se)+i1[ii1[pa]][ii]] <- 0
 #     }

      ## remove superfluous 
      
      CROPSsfrm[se,se][which(CROPSsfrm_tmp==0)] <- 0
     }
   
    }

  }

 return(CROPSsfrm)

}



find_superfluous <- function(CROPSsfrm_tmp,max_triple,r1,M,se,rmax,w,i1,i2,ii1,n,alpha)	{
	
	## Needs the superfluous C library
	 dyn.load("superfluous.so")
	
	## define segments of the list of triplets not longer than max_triple
	ll <- ceiling(length(r1)/max_triple)
	ii_s <- 1+max_triple*(seq(0,ll-1))
	ii_e <- max_triple+max_triple*(seq(0,ll-1))
	ii_e[which(ii_e > length(r1))] <- length(r1) 
	
	## find superfluous for all the segments of the list of triplets
	for(index in 1:ll)
	{
		pa <- ii_s[index]:ii_e[index]
		MU_tri <- array(0,length(i1[ii1[pa]]))
		MUsig_tri <- array(0,length(i1[ii1[pa]]))
		
		## calculate inner composition alignment for triplets, length(which(napos[count,]==1)) time points and length(se) genes involved	
		re <- .C("superfluous", as.double(M), as.integer(n) , as.integer(length(se)), as.double(MU_tri), as.integer(i1[ii1[pa]]), as.integer(i2[ii1[pa]]), as.integer(r1[pa]), as.integer(length(i1[ii1[pa]])), as.integer(rmax), as.double(MUsig_tri), as.integer(w))
		
		## find significant values
		
		ii <- which(re[[10]] > alpha)
		
		## remove superfluous in temporary matrix
		
		CROPSsfrm_tmp[(i2[ii1[pa]][ii]-1)*length(se)+i1[ii1[pa]][ii]] <- 0
	}
	return(CROPSsfrm_tmp)
}









####################################
####################################
####################################
######### to run the program #######
####################################
####################################
####################################

#rmax <- 1000   # realizations
#alpha <- 0.99   # significance
#w <- 2          # weighting: squared slope

#source("iota_subroutines.R")


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

#I <- IOTA(TimeSeries,rmax,alpha,w)











