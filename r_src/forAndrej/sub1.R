calc_act_tr <- function(TF1, const_p, CON)
{
 act_tr <- 1
 li_tf <- length(which(is.na(const_p[,1,TF1])==0))
 reg_g <- which(CON[TF1,]==1)
 li_g <- 0
 if(length(reg_g)>0)
 {
  for(i_rg in 1:length(reg_g))
  {
   if( length(which(is.na(const_p[,1,reg_g[i_rg]])==0)) > li_g)
   {
    li_g <- length(which(is.na(const_p[,1,reg_g[i_rg]])==0))
   }
  }
  for(l1 in 1:li_tf)
  {
   for(l2 in 1:li_g)
   {
    for(i_rg in 1:length(reg_g))
    {
     if(is.na(const_p[l2,4,reg_g[i_rg]])==0)
     {
      if( (const_p[l1,4,TF1] == 1) & (const_p[l2,4,reg_g[i_rg]] == 1) )
      {
       if( (const_p[l2,1,reg_g[i_rg]] >= const_p[l1,1,TF1]) & (const_p[l2,1,reg_g[i_rg]] <= const_p[l1,2,TF1]) )
       {
        if(act_tr > const_p[l1,3,TF1])
        {
         act_tr <- const_p[l1,3,TF1]
        }
       }
      }
     }
    }
   }
  }
 }
 return(act_tr)
}



#########################################################################
#########################################################################
#########################################################################


reg_direction <- function(i1,i2,xc1_p,xc1)
{
 #VS <- xc1[which(xc1_p[,i1,i2]<= 0.05),i1,i2]
 VS <- xc1[which( (xc1_p[,i1,i2]<= 0.05) & (abs(xc1[,i1,i2]) > 0.5) ),i1,i2]
 VVV <- NA
 if(length(VS > 0))
 {
  V <- ((N_shift+1)-which(abs(abs(VS)) == max(abs(VS))))
  if(length(V) > 0)
  {
   VV <- V[which(min(abs(V)) == abs(V))]
   if(length(VV) == 1)
   {
    if(VV > 0) 
    {
     # i1 may activate i2
     VVV <- 1
    }
    if(VV < 0)
    {
     # i2 may activate i1
     VVV <- 2 
    }
    if(VV == 0)
    {
     # i2 may activate i1 or otherway around
     VVV <- 0    
    }
   }
   else
   {
    if(VV == 0)
    {
     # i2 may activate i1 or otherway around
     VVV <- 0
    }
   }
  }
 }
 return(VVV)
}






