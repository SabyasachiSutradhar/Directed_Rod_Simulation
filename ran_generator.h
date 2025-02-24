/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
Part of this code is taken from Numerical receipies in C, William H. Press, Saul A. Teukolsky, William T. Vetterling,
and Brian P. Flannery, CAMBRIDGE UNIVERSITY PRESS
Copyright @ Sabyasachi Sutradhar
*/
  ///************************************************************************//
void nrerror0(const char error_text[])
{
        printf("Numerical Recipes run-time error...\n");
        printf("%s\n",error_text);
        printf("...now exiting to system...\n");
        exit(1);
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran2(long *idum){
int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
double temp;
 if (*idum <= 0) { //Initialize.
    if (-(*idum) < 1) *idum=1; //Be sure to prevent idum = 0.
else *idum = -(*idum);
idum2=(*idum);
 for (j=NTAB+7;j>=0;j--) { //Load the shuffle table (after 8 warm-ups).
k=(*idum)/IQ1;
*idum=IA1*(*idum-k*IQ1)-k*IR1;
if (*idum < 0) *idum += IM1;
if (j < NTAB) iv[j] = *idum;
}
iy=iv[0];
}
 k=(*idum)/IQ1; //Start here when not initializing.
 *idum=IA1*(*idum-k*IQ1)-k*IR1; //Compute idum=(IA1*idum) % IM1 without
 if (*idum < 0) *idum += IM1; //overflows by Schrage’s method.
k=idum2/IQ2;
 idum2=IA2*(idum2-k*IQ2)-k*IR2;// Compute idum2=(IA2*idum) % IM2 likewise.
if (idum2 < 0) idum2 += IM2;
 j=iy/NDIV; //Will be in the range 0..NTAB-1.
 iy=iv[j]-idum2; //Here idum is shuffled, idum and idum2 are
 iv[j] = *idum; //combined to generate output.
if (iy < 1) iy += IMM1;
 if ((temp=AM*iy) > RNMX) return RNMX; //Because users don’t expect endpoint values.
else return temp;
}


//**************************************************************************//
double gasdev()
//Returns a normally distributed deviate with zero mean and unit variance, using ran2(idum) as the source of uniform deviates.
{
static int iset=0;
static double gset;
double fac,rsq,v1,v2;
if (&idum < 0) iset=0;
if (iset == 0) {
do {
  v1=2.0*ran2(&idum)-1.0;
  v2=2.0*ran2(&idum)-1.0;
  rsq=v1*v1+v2*v2;
} while (rsq >= 1.0 || rsq == 0.0);
fac=sqrt(-2.0*log(rsq)/rsq);
gset=v1*fac;
iset=1;
return v2*fac;
 } else {
    iset=0;
return gset;
}
}



double gaussdev(double mu,double sigma)
//Returns a normally distributed deviate with zero mean and unit variance, using ran2(idum)as the source of uniform deviates.
{
  if(mu==0.0 && sigma==0.0){
    return 0.0;
}else{
return mu+sigma*gasdev();
}
}

double logndev(double mu,double sigma)
//Returns a normally distributed deviate with zero mean and unit variance, using ran2(idum)as the source of uniform deviates.
{
return exp(mu+sigma*gasdev());
}


double expdev(double a){
double dum;
do{
dum=ran2(&idum);
}while (dum == 0.0);
 return (-a*log(dum));
}
