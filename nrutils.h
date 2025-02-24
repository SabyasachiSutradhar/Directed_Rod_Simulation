#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NR_END 1
#define FREE_ARG char*
#define TINY 1.0e-20
static double dmaxarg1,dmaxarg2;

double sqr(double q){
  return (q*q);
}
double cube(double p){
  return(p*p*p);
}
double cube_root(double p){
  return(pow(p,1.0/3.0));
}
int find_min(int i,int j){
  if(i>=j) return(j);
  else return(i);
}
int find_max(int i,int j){
  if(i>=j) return(i);
  else return(j);
}
int find_dmin(double i,double j){
  if(i<j) return(i);
  else return(j);
}
int find_dmax(double i,double j){
  if(i>=j) return(i);
  else return(j);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(char const error_text[])
/* Numerical Recipes standard error handler */ {
printf("Numerical Recipes run-time error...\n");
printf("%s\n",error_text);
printf("...now exiting to system...\n");
exit(1);
}
////////////////////////////////////////
double ***dtensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) /* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1; double ***t;
/* allocate pointers to pointers to rows */
t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**))); if (!t) nrerror("allocation failure 1 in f3tensor()");
t += NR_END;
t -= nrl;
/* allocate pointers to rows and set pointers to them */
t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*))); if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
t[nrl] += NR_END;
t[nrl] -= ncl;
/* allocate rows and set pointers to them */
t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double))); if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
 t[nrl][ncl] += NR_END; t[nrl][ncl] -= ndl;
for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep; for(i=nrl+1;i<=nrh;i++) {
t[i]=t[i-1]+ncol; t[i][ncl]=t[i-1][ncl]+ncol*ndep; for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
}
/* return pointer to array of pointers to rows */
return t; }
//////////////////////////////////////////////////////

int ***itensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) /* allocate a int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1; int ***t;
/* allocate pointers to pointers to rows */
t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**))); if (!t) nrerror("allocation failure 1 in f3tensor()");
t += NR_END;
t -= nrl;
/* allocate pointers to rows and set pointers to them */
t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*))); if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
t[nrl] += NR_END;
t[nrl] -= ncl;
/* allocate rows and set pointers to them */
t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int))); if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");

 t[nrl][ncl] += NR_END; t[nrl][ncl] -= ndl;
for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep; for(i=nrl+1;i<=nrh;i++) {
t[i]=t[i-1]+ncol; t[i][ncl]=t[i-1][ncl]+ncol*ndep; for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
}
/* return pointer to array of pointers to rows */
return t; }

//////////////////////////////////////////////////////////////////////////////////
double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */ {
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1; double **m;
/* allocate pointers to rows */
m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*))); if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double))); if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
return m; }
//////////////////////////////////////////////

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */ {
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1; int **m;
/* allocate pointers to rows */
m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*))); if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int))); if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
return m; }

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */ {
double *v;
v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); if (!v) nrerror("allocation failure in dvector()");
return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */ {
int *v;
v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int))); if (!v) nrerror("allocation failure in ivector()");
return v-nl+NR_END;
}



void free_dvector(double *v, long nl, long nh)
/* free a float vector allocated with vector() */ {
free((FREE_ARG) (v+nl-NR_END)); }
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */ {
free((FREE_ARG) (v+nl-NR_END)); }

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) /* free a double matrix allocated by dmatrix() */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END)); }

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) /* free an int matrix allocated by imatrix() */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END)); }

void free_dtensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */ {
free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END)); free((FREE_ARG) (t[nrl]+ncl-NR_END)); free((FREE_ARG) (t+nrl-NR_END));
}

void free_itensor(int***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */ {
free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END)); free((FREE_ARG) (t[nrl]+ncl-NR_END)); free((FREE_ARG) (t+nrl-NR_END));
}

double distance2D(double xi,double yi,double xf,double yf){
 return(sqrt((xi-xf)*(xi-xf)+(yi-yf)*(yi-yf)));
}
double distance(double xi,double yi,double zi ,double xf,double yf,double zf){
 return(sqrt(sqr(xi-xf)+sqr(yi-yf)+sqr(zi-zf)));
}
double random_sign(){
  if(ran2(&idum)>0.5) return(1.0);
  else return(-1.0);
}
double Volume(double r){
  return 4.0*pi*r*r*r/3.0;
}
double Surface_area(double r){
  return 4.0*pi*r*r;
}

double wrap_angle(double theta){
  if(theta<0.0) {theta=2.0*pi+theta;}
  if(theta>2.0*pi){theta=fmod(theta,2.0*pi);}
  return(theta);
}
int random_int(int a,int b){
  return (a+floor(ran2(&idum)*(b-a+1)));
}

int lineseg_intersection(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4){
int ind=0;
double den=(x4-x3)*(y1-y2)-(x1-x2)*(y4-y3);
if(den!=0.0000){
double ta,tb;
ta=((y3-y4)*(x1-x3)+(x4-x3)*(y1-y3))/den;
tb=((y1-y2)*(x1-x3)+(x2-x1)*(y1-y3))/den;
if(ta>=0 && ta<=1){
if(tb>=0 && tb<=1){
ind=1;
}
}
}
return ind;
}

double **equispace_points(double *x,double *y,int N1,int N2){
  double **p,**q;
  p=dmatrix(0,N1,0,2);
  q=dmatrix(0,N2,0,2);
  for(int i=1;i<=N1;i++){
    p[i][1]=x[i];
    p[i][2]=y[i];
  }
  double *currentpt,totaldist,intv,disttmp,distsum,*ptnow,*pttarget,remainder,*newpt;
  currentpt=dvector(0,2);
  ptnow=dvector(0,2);
  pttarget=dvector(0,2);
  newpt=dvector(0,2);

  double l,m,norm;
  int indfirst,len,kk,ind;
  indfirst=2;

  for(unsigned int j=1;j<=2;j++){
    currentpt[j]=p[1][j];
    q[1][j]=currentpt[j];
  }

  totaldist=0.0;
  for(int i=1;i<=N1-1;i++){
    totaldist+=distance2D(p[i][1],p[i][2],p[i+1][1],p[i+1][2]);
  }
  intv=totaldist/(double)(N2-1);
///////////////////////////////////////////////
  for(int k=1;k<=N2-1;k++){
//printf("%d %d\n",k,indfirst);
    for(int j=1;j<=2;j++){
      ptnow[j]=currentpt[j];
      pttarget[j]= p[indfirst][j];
    }
    remainder = intv;
    kk = 0;
    ind=0;
    distsum = 0.0;
    while(ind==0){
      disttmp=distance2D(ptnow[1],ptnow[2],pttarget[1],pttarget[2]);
      distsum  = distsum+disttmp;
      if(distsum >= intv){
        norm=distance2D(ptnow[1],ptnow[2],pttarget[1],pttarget[2]);
        l=(pttarget[1]-ptnow[1])/norm;
        m=(pttarget[2]-ptnow[2])/norm;
        newpt[1]=remainder*l+ptnow[1];newpt[2]=remainder*m+ptnow[2];
        ind=1;
      }else{
        remainder = remainder - disttmp;
        for(int j=1;j<=2;j++){ptnow[j] = pttarget[j];}
        kk=kk+1;
        if(indfirst+kk>N1){
          for(int j=1;j<=2;j++){newpt[j] = p[N1][j];}
                  ind=1;
        }else{
          for(int j=1;j<=2;j++){pttarget[j] = p[indfirst+kk][j];}
        }
      }
    }

    for(int j=1;j<=2;j++){
      q[k+1][j]=newpt[j];
      currentpt[j] = newpt[j];
    }
  indfirst = indfirst + kk;
  }
return q;
free_dmatrix(p,0,N1,0,2);
free_dmatrix(q,0,N2,0,2);
free_dvector(currentpt,0,2);
free_dvector(ptnow,0,2);
free_dvector(pttarget,0,2);
free_dvector(newpt,0,2);
}

void ludcmp(double **a, int n, int *indx, double *d){
/*Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation effected by the partial
pivoting; d is output as ±1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.*/
int i,imax,j,k;
double big,dum,sum,temp;
double *vv; //vv stores the implicit scaling of each row.
vv=dvector(1,n);
*d=1.0;
for (i=1;i<=n;i++) {
  big=0.0;
  for (j=1;j<=n;j++)
  if ((temp=fabs(a[i][j])) > big) big=temp;
  if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
  vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
  for (i=1;i<j;i++) { sum=a[i][j];
    for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
    a[i][j]=sum;
  }
big=0.0;
for (i=j;i<=n;i++) {
  sum=a[i][j];
  for (k=1;k<j;k++)
  sum -= a[i][k]*a[k][j]; a[i][j]=sum;
if ( (dum=vv[i]*fabs(sum)) >= big) {
  big=dum;
imax=i;
}
}
if (j != imax) {
for (k=1;k<=n;k++) { dum=a[imax][k];
a[imax][k]=a[j][k];
        a[j][k]=dum;
    }
    *d = -(*d);
    vv[imax]=vv[j];
  }
    indx[j]=imax;
if (a[j][j] == 0.0) a[j][j]=TINY;
if (j != n) {
  dum=1.0/(a[j][j]);
for (i=j+1;i<=n;i++) a[i][j] *= dum;
}
}
free_dvector(vv,1,n);

}


void lubksb(double **a, int n, int *indx, double b[]){
/*Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, and indx are not modified by this routine
and can be left in place for successive calls with different right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion.*/
int i,ii=0,ip,j;
double sum;
for (i=1;i<=n;i++) { ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
      }
      for (i=n;i>=1;i--) {
        sum=b[i];
for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
b[i]=sum/a[i][i];
}

}



void Solve_Lin(double **a,int N,double *b){
double d;
int *indx;
indx=ivector(1,N);

ludcmp(a,N,indx,&d); //Decompose the matrix just once.
lubksb(a,N,indx,b);

}

/////////////////////// Matrix inverse
void matrix_inv(double **a,double **ainv,int N){
double d,*col;
int *indx;

col=dvector(1,N);
indx=ivector(1,N);


ludcmp(a,N,indx,&d); //Decompose the matrix just once.

for(int j=1;j<=N;j++) { //Find inverse by columns.
for(int i=1;i<=N;i++) col[i]=0.0;
col[j]=1.0;
lubksb(a,N,indx,col);
for(int i=1;i<=N;i++) ainv[i][j]=col[i];
}

free_dvector(col,1,N);
free_ivector(indx,1,N);
}
////////////////////////////////////////////////////////Determinant
double det(double **a,int N){
double d;
int *indx;
indx=ivector(1,N);

ludcmp(a,N,indx,&d); //This returns d as ±1.
for(int j=1;j<=N;j++) {
d*= a[j][j];
}
free_ivector(indx,1,N);
return d;
}

/////////////////////////////////////////////////// Add matrices c=alpha+a alpha being scalar
void AddScalMat(double **a,double alpha,int n1,int n2){
  for (int i=1;i<=n1;i++){
    for(int j=1;j<=n2;j++){
      a[i][j]=alpha+a[i][j];
    }
  }
}
/////////////////////////////////////////////////// Add matrices c=alpha*a+beta*b, alpha & beta being scalar
void AddMatMat(double **a,double **b,double **c,int n1,int n2,double alpha,double beta){

  for (int i=1;i<=n1;i++){
    for(int j=1;j<=n2;j++){
      c[i][j]=alpha*a[i][j]+beta*b[i][j];
    }
  }
}
/////////////////////////////////////////////////// Multiply matriceby a scalar c=alpha*a
void MulScalVec(double *a,double alpha,int n1){
  for (int i=1;i<=n1;i++){
      a[i]=alpha*a[i];
  }
}

void MulScalMat(double **a,double alpha,int n1,int n2){
  for (int i=1;i<=n1;i++){
    for(int j=1;j<=n2;j++){
      a[i][j]=alpha*a[i][j];
    }
  }
}

void MulMatMatElem(double **a,double **b,double **c,int n1,int n2){

    for (int i=1;i<=n1;i++){
      for(int j=1;j<=n2;j++){
        c[i][j]=a[i][j]*b[i][j];
      }
    }

}

void MulVecVecElem(double *a,double *b,double *c,int n1){
    for (int i=1;i<=n1;i++){
        c[i]=a[i]*b[i];
    }

}

/////////////////////////////////////////////////// Multiply two matrices c=alpha*a
void MulMatMat(double **a,double **b,double **c,int n1,int n2,int n3){
  for (int i=1;i<=n1;i++){
    for(int j=1;j<=n3;j++){
          c[i][j]=0.0;
      for(int k=1;k<=n2;k++){
      c[i][j]+=a[i][k]*b[k][j];
    }
  }
}

}
///////////////////// transpose matrix //////////////////////////////
void TransposeMat(double **a,double **c,int n1,int n2){
  for (int i=1;i<=n1;i++){
    for(int j=1;j<=n2;j++){
      c[j][i]=a[i][j];
    }
  }

}
/////////////////// outer product
void OuterProductMat(double **mat,double **c,int n1,int n2){

for (int i=1;i<=n1;i++){
  for(int j=1;j<=n1;j++){
        c[i][j]=0.0;
    for(int k=1;k<=n2;k++){
    c[i][j]+=mat[i][k]*mat[j][k];
  }
}
}
}

void OuterProductOneRow(double **vec,double **c,int fixedrow,int n1){
for (int i=1;i<=n1;i++){
  for(int j=1;j<=n1;j++){
    c[i][j]=vec[fixedrow][i]*vec[fixedrow][j];
}
}
}
double mod_vector(double **v,int i,int d){
  double dx=0.0;
for (int j=1;j<=d;j++){
dx+=v[i][j]*v[i][j];
}
return(sqrt(dx));
}

double KronerkarDelta(int i,int j){
  if(i==j)return 1.0;
  else return 0.0;
}
