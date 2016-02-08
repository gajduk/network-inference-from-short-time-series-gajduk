#include <R.h>
void prob_xyz(double *X, double *Y, double *Z, int *L, double *px1, double *px2,
double *px3, double *px4, int *ii, int *gl)
{  
 int i,j;
 for(i=0; i < *L; i++)
 { 
  if(gl[i] == 0)
  {  
   px1[*ii] = X[i];
   px2[*ii] = Y[i]; 
   px3[*ii] = Z[i];    
   gl[i] = 1;
   for(j=0; j < *L; j++)
   {
    if(X[i] == X[j] & Y[i] == Y[j] & Z[i] == Z[j]) 
    {
     px4[*ii] = px4[*ii]+1;
     gl[j] = 1;
    }
   }
   px4[*ii] = px4[*ii] / *L;
   *ii += 1;
  }
 }
}
 

