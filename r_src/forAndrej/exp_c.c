#include <R.h>
void exp_c(double *py1, double *py2, double *pxy1, double *pxy2, double *pxy3, int *Lpy, int *Lpxy, double *M12)
{ 
 int n, i, j;
 double pxy_s, py_s;
 for(n = 0; n < *Lpxy; n++)                        
 {
  pxy_s = 0.;
  py_s = 0.;
  for(j = 0; j < *Lpxy; j++)
  {
   for(i = 0; i < *Lpxy; i++)
   {
    if(pxy1[i] == pxy1[n] & pxy2[i] == pxy2[j])  pxy_s += pxy3[i];
   }
   for(i = 0; i < *Lpy; i++)
   {
    if(py1[i] == pxy2[j])  py_s += py2[i];
   }
  }
  *M12 +=  (pxy1[n] * (pxy_s/py_s));
 }
}



