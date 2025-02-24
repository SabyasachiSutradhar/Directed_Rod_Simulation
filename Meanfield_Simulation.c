/*
C code written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates a mean field approximation of branching morphogensis of Drosophila class-IV dendritic arbor.
branches are generates from a random point with rate k inside a box of length L and then allowed to grow with a growth velocity=drift velocity
upon collision the colliding branch is instanteneously deleted.
Copyright @ Sabyasachi Sutradhar
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
static long idum=-123456791;
#include "ran_generator.h"
double pi=3.14159265;
double PI=3.14159265;
#include "nrutils.h"
#define NUMELEM(a) (sizeof(a) / sizeof(*a))
///////////////////////////////// Variables //////////////////////////////////
int Max_Tips=100000;
/////////////////////////////////////////////////////////////////////////////
int N_Tip,N_Sample,NG,NP,NS,CollisionCount,CollisionCountMin;
double *Dendrite_Length,*Xbase,*Ybase,*Xtip,*Ytip,*Dendrite_Direction,Total_Length,Dendrite_Density,Tip_Density;
int *Deletion_index;
double Tip_SS;
double RescueAfterDeath;
double Length_SS,Tip_DensitySS,Dendrite_DensitySS,Avg_LengthSS,NGS,NPS,NSS,CollisionCountSS;
double L_Box,Branching_Rate,Tip_GrowthVelocity;
double T_Max,time,Dt,Area;
/////////////////////// Dynamical parameters //////////////
double Kgp,Kgs,Kpg,Kps,Ksg,Ksp;
double Vg_mean,Vg_sig,Vs_mean,Vs_sig,Vp_sig,VG,VS;
int *state;
double *v;
//////////////////////////////////// initialization Parameters
double Initial_TipDensity;
int N_TipInitial;
////////////////////////////////////////////////////
//////////////////////// File Configuration /////////////////////////////////
FILE *fconf,*fp1,*fp2,*fp3;
char fileconf[1000],file1[1000],file2[1000],file3[1000];
////////////////////////////////////////////////////////////////////////////
////////////////////////////// essential subroutines ///////////////////////
#include "essentials.h"
//////////////////////////////////////////////////
/////////////////// Main function //////////////////
////////////////////////////////////////////////////
int main(int argc, const char *argv[]){

#include "mem_allocation.h"//allocate memory for variables
////////////////////////////// Parameters //////////////////////////////////
L_Box=200.0;// in micron
Area=L_Box*L_Box;
Branching_Rate=0.00954;// /um/micron/
Tip_GrowthVelocity=0.034;////um/min
Initial_TipDensity=0.035580866;
N_TipInitial=ceil(Initial_TipDensity*L_Box*L_Box);
/////////////////////// dynamical parameters ////////////////////////////////
Kgp=0.7839;  Kgs=0.6400;
Kpg=0.3349; Kps=0.3137;
Ksg=0.5978; Ksp=0.9460;

Vg_mean=0.4102; Vg_sig=0.3595;
Vs_mean=0.3543; Vs_sig=0.3749;
Vp_sig=0.3385;

VG=exp(Vg_mean+0.5*Vg_sig*Vg_sig);
VS=exp(Vs_mean+0.5*Vs_sig*Vs_sig);

RescueAfterDeath=(double) 0 *0.1;
/////////////////////////////////////////////////////////////////////////////
T_Max=5000.0;///minutes
Dt=0.025;
int im_mod=(int)(500.0/Dt+TINY);
int ss_counter;
//////////////////////////////////////////////////////////////////
N_Sample=25;
//////////////////////////////////////////////////////////////////
int sample=1;
//for(sample=1;sample<=N_Sample;sample++){

N_Tip=N_TipInitial;
initialize();

sprintf(file1,"TimeData-Rescue%.2f_Sample-%d.dat",RescueAfterDeath,sample);
fp1=fopen(file1,"w");
fprintf(fp1,"# Simulated in a box of size=%f um\n",L_Box);
fprintf(fp1,"#Parameters:Kb=%f;KGP=%f;KGS=%f;KPG=%f;KPS=%f;KSG=%f;KSP=%f;VG=%f,VS=%f,Beta=%f\n",Branching_Rate,Kgp,Kgs,Kpg,Kps,Ksg,Ksp,VG,VS,RescueAfterDeath);
fprintf(fp1,"# All values are in microns and minutes\n");
fprintf(fp1,"#time(min) Tip Total_Length Tip_Density Length_Density Avg_Length CollisionCount CollisionRate NG NP NS\n");
fclose(fp1);

sprintf(file2,"SteadyStateValues-Rescue%.2f_Sample-%d.dat",RescueAfterDeath,sample);
fp2=fopen(file2,"w");
fprintf(fp2,"# Simulated in a box of size=%f um\n",L_Box);
fprintf(fp2,"#Parameters:Kb=%f;KGP=%f;KGS=%f;KPG=%f;KPS=%f;KSG=%f;KSP=%f;VG=%f,VS=%f,Beta=%f\n",Branching_Rate,Kgp,Kgs,Kpg,Kps,Ksg,Ksp,VG,VS,RescueAfterDeath);
fprintf(fp2,"# All values are in microns and minutes\n");
fprintf(fp2,"# Sample Tip Total_Length Tip_Density Length_Density Avg_Length CollisionRate NG NP NS\n");
fclose(fp2);
//////////////// Time evolution /////////////
time=0.0;
int m=0;
int im_counter=1000;
CollisionCountMin=0;
//////////////////////
ss_counter=0;
Tip_SS=0.0;;
Length_SS=0.0;
Tip_DensitySS=0.0;
Dendrite_DensitySS= 0.0;
Avg_LengthSS=0.0;
NGS=0.0;
NPS=0.0;
NSS=0.0;
CollisionCountSS=0.0;
CollisionCountMin=0;
////////////////////////

do{
tip_dynamics();
Remove_DeletedDendrites();
CollisionCountMin+=CollisionCount;

if(m%im_mod==0){
  print_conf(sample,time,RescueAfterDeath);
}

if(m%(int)(1.0/Dt)==0){
  print_data(time);
}
///////////////////////// Calculate steady state densities ////////////////////
if (time>T_Max-200.0){
  ss_counter++;
  Tip_SS+=(double) N_Tip;
  Length_SS+=Total_Length;
  Tip_DensitySS+=(double) N_Tip/Area;
  Dendrite_DensitySS+= Total_Length/Area;
  Avg_LengthSS+=Total_Length/(double)N_Tip;
  CollisionCountSS+=(double)CollisionCount/(Dt*Area);
  NGS+=(double)NG;
  NPS+=(double)NP;
  NSS+=(double)NS;
}
////////////////////////////////////////////////////////////////////////////
time+=Dt;
m++;
}while(time<=T_Max);
//////////////////////// Print the steady state values ////////////////////////
Tip_SS=Tip_SS/(double) ss_counter;
Length_SS=Length_SS/(double) ss_counter;
Tip_DensitySS=Tip_DensitySS/(double) ss_counter;
Dendrite_DensitySS= Dendrite_DensitySS/(double) ss_counter;
Avg_LengthSS=Avg_LengthSS/(double) ss_counter;
CollisionCountSS=CollisionCountSS/(double) ss_counter;
NGS=NGS/(double) ss_counter;
NPS=NPS/(double) ss_counter;
NSS=NSS/(double) ss_counter;

fp2=fopen(file2,"a");
fprintf(fp2,"%d %f %f %f %f %f %f %f %f %f\n",sample,Tip_SS,Length_SS,Tip_DensitySS,Dendrite_DensitySS,Avg_LengthSS,CollisionCountSS,NGS/Area,NPS/Area,NSS/Area);
fclose(fp2);
//////////////////////////////////////////////////////////////////////////////
//}
////////////////// sample loop ends

}

