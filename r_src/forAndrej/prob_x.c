#include <R.h>
void prob_x(double *X, int *L, double *px1, double *px2, int *ii, int *gl)
{  
 int i,j;
 for(i=0; i < *L; i++)
 { 
  if(gl[i] == 0)
  {  
   px1[*ii] = X[i];
   gl[i] = 1;
   for(j=0; j < *L; j++)
   {
    if(X[i] == X[j]) 
    {
     px2[*ii] = px2[*ii]+1;
     gl[j] = 1;
    }
   }
   px2[*ii] = px2[*ii] / *L;
   *ii += 1;
  }
 }
}
 

