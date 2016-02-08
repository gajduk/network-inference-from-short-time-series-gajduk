#include <R.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "rngs.h"
#include "rvgs.h"


void quickSort( double a[], int l, int r, int o[]);
int partition( double a[], int l, int r, int o[]) ;
double cross( double g[], int n, double delta, int w); 
double Random(void);
void PlantSeeds(long x);
void PutSeed(long x);
void GetSeed(long *x);
void SelectStream(int index);
void TestRandom(void);
double Uniform(double a, double b);
void sample( int s[], int n);
void crops_p(double MM[], int nn, int mm, double mu[], double delta, int o_init[], int gg_i[], int gg_j[], int gg_k[], int w);

 
void superfluous(double *M, int *n, int *m, double *MU, int *g_i, int *g_j, int *g_k, int *l_tri, int *rmax, double *MUsig, int *weight)
{ 
 int i, j, o_init[*n], c_sig[*l_tri], r, s[*n];
 double delta, M_r[(*n * *m)], MU_r[*l_tri];
 
 /*
 #############################################################
 initial sequence (permutation)
 #############################################################
 */ 
 for(i = 0; i < *n; i++)
 {
  o_init[i] = (i+1);
 }
 /*
 #############################################################
 normalization constant
 #############################################################
 */  
 delta = (double)(*n-1) * (double)(*n-2) / 2.; 
 /*
 #############################################################
 calculate crossing points partial 
 #############################################################
 */  
 crops_p(M, *n, *l_tri, MU, delta, o_init, g_i, g_j, g_k, *weight);
 /*
 #############################################################
 initialize sig counter
 #############################################################
 */ 
 for(i = 0; i <  *l_tri; i++)
 {
  c_sig[i] = 0;   
 }
 /*
 #############################################################
 calculate significance of crossing points measure
 #############################################################
 */ 
 for(r = 0; r < *rmax; r++)
 {
  printf(" %d \n ", r);
  for(j = 0; j < *m; j++)
  { 
   for(i = 0; i < *n; i++)
   {
    s[i] = (o_init[i]-1);
   }  
   /*
   #############################################################
   calculate random permutation
   #############################################################
   */
   sample( s, *n);
   /*
   #############################################################
   randomize time series 
   #############################################################
   */   
   for(i = 0; i < *n; i++)
   {
    M_r[((j * *n) + i)] = M[((j * *n) + s[i])];
   }     
  }
  /*
  #############################################################
  calculate crossing points for random series
  #############################################################
  */  
  crops_p(M_r, *n, *l_tri, MU_r, delta, o_init, g_i, g_j, g_k, *weight);    
  for(i = 0; i < *l_tri; i++)
  { 
   /*
   #############################################################
   significance of superfluity
   #############################################################
   */  
   c_sig[i] += (fabs(MU_r[i]) >= fabs(MU[i])) ? 1 : 0;
   /* c_sig[i] += (MU_r[i] <= MU[r]) ? 1 : 0;  */ 
  }
 }
 for(i = 0; i < *l_tri; i++)
 {
  MUsig[i] = (double)(c_sig[i]+1)/(double)(*rmax+1);
 }
}



/*
#############################################################################################
#############################################################################################
######################################### subroutines  ######################################
#############################################################################################
#############################################################################################
*/

/*
#############################################################################################
########################### routines for measure evaluation #################################
#############################################################################################
*/
void crops_p(double MM[], int nn, int mm, double mu[], double delta, int o_init[], int gg_i[], int gg_j[], int gg_k[], int w)
{ 
 int i, o1[nn], o2[nn], l;
 double a[nn], g1[nn], g2[nn];
 
 /*
 #############################################################
 loop for triplets k -> i, k -> j, i -> j
 #############################################################
 */
 for(l = 0; l < mm; l++)
 { 
  /*
  #############################################################
  sort k
  #############################################################
  */ 
  for(i = 0; i < nn; i++)
  {
   a[i] = MM[((gg_k[l]-1) * nn)+i];
   o1[i] = o_init[i];
  }  
  quickSort( a, 0, (nn-1), o1);
  /*
  #############################################################
  sort i
  #############################################################
  */   
  for(i = 0; i < nn; i++)
  {
   a[i] = MM[((gg_i[l]-1) * nn)+i];
   o2[i] = o_init[i];
  }
  quickSort( a, 0, (nn-1), o2);  
  /*
  #############################################################
  reorder j according to i[k] and according to i
  #############################################################
  */  
  for(i = 0; i < nn; i++)
  {
   g1[i] = MM[(((gg_j[l]-1) * nn) + o2[o1[i]-1] - 1)];
   g2[i] = MM[(((gg_j[l]-1) * nn) + o2[i] - 1)];
  } 
  /*
  #############################################################
  calculate difference between crops measure including i and k and i only 
  #############################################################
  */     
  mu[l] = cross( g1, nn, delta, w) - cross( g2, nn, delta, w);
 }  
}
/*
#############################################################################################
#############################################################################################
*/
#define MAX(a, b) ((a) > (b) ? (a) : (b))
double cross( double g[], int n, double delta, int w) 
{
  int tau,s;
  double omega = 0., mu;
     
  /* no weight */ 
  if(w == 0)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? 1. : 0.;  
    }
   } 
  }
  /* slope */  
  if(w == 1)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? (fabs(g[s+1]-g[s])) : 0.;  
    }
   } 
  } 
  /* squared slope */   
  if(w == 2)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? ( (g[s+1]-g[s])*(g[s+1]-g[s]) ) : 0.;  
    }
   } 
  }
  /* arithmetic mean */   
  if(w == 3)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? ( (g[s+1]+g[s])/2. ) : 0.;  
    }
   } 
  }
  /* geometric mean */  
  if(w == 4)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? ( sqrt(g[s+1] * g[s]) ) : 0.;  
    }
   } 
  } 
  /* harmonic mean */
  if(w == 5)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? ( (2./((1./g[s+1])+(1./g[s]))) ) : 0.;  
    }
   } 
  }   
  /* max */
  if(w == 6)
  {
   for(tau = 0; tau < (n-2); tau++)
   {
    for(s = (tau+1); s < (n-1); s++)
    {
     omega += ( (g[s+1] > g[tau] & g[s] < g[tau]) | (g[s+1] < g[tau] & g[s] > g[tau]) ) ? ( (MAX((fabs(g[s+1] - g[tau])), (fabs(g[s] - g[tau])))) ) : 0.;  
    }
   } 
  }      

  mu = 1. - (omega/delta);
  return(mu); 
}
/*
#############################################################################################
########################### routines for sorting process ####################################
#############################################################################################
*/
/*
#############################################################################################
#############################################################################################
*/
void quickSort( double a[], int l, int r, int o[])
{
   int j;

   if( l < r ) 
   {
       // divide and conquer
       j = partition( a, l, r, o);
       quickSort( a, l, j-1, o);
       quickSort( a, j+1, r, o);
   }
	
}
/*
#############################################################################################
#############################################################################################
*/
int partition( double a[], int l, int r, int o[]) 
{
   int i, j, tt; 
   double pivot, t;
   pivot = a[l];
   i = l; j = r+1;
		
   while( 1)
   {
   	do ++i; while( a[i] <= pivot && i <= r );
   	do --j; while( a[j] > pivot );
   	if( i >= j ) break;
   	t = a[i]; a[i] = a[j]; a[j] = t;
  	t = o[i]; o[i] = o[j]; o[j] = t;	
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   tt = o[l]; o[l] = o[j]; o[j] = tt;  
   return j;
}
/*
#############################################################################################
########################### routines for randomization process ##############################
#############################################################################################
*/
/*
#############################################################################################
#############################################################################################
*/
void sample( int s[], int n)
{
 int i;
 int nn = n;
 int t, tt; 
 
 for(i = 0; i < n; i++)
 {
  nn--;
  t = (int)round(Uniform(0., (double)nn));
  tt = s[t]; s[t] = s[nn]; s[nn] = tt;
 }
}
/*
#############################################################################################
#############################################################################################
*/
#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
      
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */
/*
#############################################################################################
#############################################################################################
*/
   double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0) 
    seed[stream] = t;
  else 
    seed[stream] = t + MODULUS;
  return ((double) seed[stream] / MODULUS);
}
/*
#############################################################################################
#############################################################################################
*/
   void PlantSeeds(long x)
/* ---------------------------------------------------------------------
 * Use this function to set the state of all the random number generator 
 * streams by "planting" a sequence of states (seeds), one per stream, 
 * with all states dictated by the state of the default stream. 
 * The sequence of planted states is separated one from the next by 
 * 8,367,782 calls to Random().
 * ---------------------------------------------------------------------
 */
{
  const long Q = MODULUS / A256;
  const long R = MODULUS % A256;
        int  j;
        int  s;

  initialized = 1;
  s = stream;                            /* remember the current stream */
  SelectStream(0);                       /* change to stream 0          */
  PutSeed(x);                            /* set seed[0]                 */
  stream = s;                            /* reset the current stream    */
  for (j = 1; j < STREAMS; j++) {
    x = A256 * (seed[j - 1] % Q) - R * (seed[j - 1] / Q);
    if (x > 0)
      seed[j] = x;
    else
      seed[j] = x + MODULUS;
   }
}
/*
#############################################################################################
#############################################################################################
*/
   void PutSeed(long x)
/* ---------------------------------------------------------------
 * Use this function to set the state of the current random number 
 * generator stream according to the following conventions:
 *    if x > 0 then x is the state (unless too large)
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively
 * ---------------------------------------------------------------
 */
{
  char ok = 0;

  if (x > 0)
    x = x % MODULUS;                       /* correct if x is too large  */
  if (x < 0)                                 
    x = ((unsigned long) time((time_t *) NULL)) % MODULUS;              
  if (x == 0)                                
    while (!ok) {
      printf("\nEnter a positive integer seed (9 digits or less) >> ");
      scanf("%ld", &x);
      ok = (0 < x) && (x < MODULUS);
      if (!ok)
        printf("\nInput out of range ... try again\n");
    }
  seed[stream] = x;
}
/*
#############################################################################################
#############################################################################################
*/
   void GetSeed(long *x)
/* ---------------------------------------------------------------
 * Use this function to get the state of the current random number 
 * generator stream.                                                   
 * ---------------------------------------------------------------
 */
{
  *x = seed[stream];
}
/*
#############################################################################################
#############################################################################################
*/
   void SelectStream(int index)
/* ------------------------------------------------------------------
 * Use this function to set the current random number generator
 * stream -- that stream from which the next random number will come.
 * ------------------------------------------------------------------
 */
{
  stream = ((unsigned int) index) % STREAMS;
  if ((initialized == 0) && (stream != 0))   /* protect against        */
    PlantSeeds(DEFAULT);                     /* un-initialized streams */
}
/*
#############################################################################################
#############################################################################################
*/
   void TestRandom(void)
/* ------------------------------------------------------------------
 * Use this (optional) function to test for a correct implementation.
 * ------------------------------------------------------------------    
 */
{
  long   i;
  long   x;
  double u;
  char   ok = 0;  

  SelectStream(0);                  /* select the default stream */
  PutSeed(1);                       /* and set the state to 1    */
  for(i = 0; i < 10000; i++)
    u = Random();
  GetSeed(&x);                      /* get the new state value   */
  ok = (x == CHECK);                /* and check for correctness */

  SelectStream(1);                  /* select stream 1                 */ 
  PlantSeeds(1);                    /* set the state of all streams    */
  GetSeed(&x);                      /* get the state of stream 1       */
  ok = ok && (x == A256);           /* x should be the jump multiplier */    
  if (ok)
    printf("\n The implementation of rngs.c is correct.\n\n");
  else
    printf("\n\a ERROR -- the implementation of rngs.c is not correct.\n\n");
}
/*
#############################################################################################
#############################################################################################
*/
   double Uniform(double a, double b)
/* =========================================================== 
 * Returns a uniformly distributed real number between a and b. 
 * NOTE: use a < b
 * ===========================================================
 */
{ 
  return (a + (b - a) * Random());
}
