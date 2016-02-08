#include <R.h>
void prob_xy(double *X, double *Y, int *L, double *px1, double *px2, double *px3, int *ii, int *gl)
{  
 int i,j;
 for(i=0; i < *L; i++)
 { 
  if(gl[i] == 0)
  {  
   px1[*ii] = X[i];
   px2[*ii] = Y[i];   
   gl[i] = 1;
   for(j=0; j < *L; j++)
   {
    if(X[i] == X[j] & Y[i] == Y[j]) 
    {
     px3[*ii] = px3[*ii]+1;
     gl[j] = 1;
    }
   }
   px3[*ii] = px3[*ii] / *L;
   *ii += 1;
  }
 }
}
 

