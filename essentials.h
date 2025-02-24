/*
C code written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates a mean field approximation of branching morphogensis of Drosophila class-IV dendritic arbor.
branches are generates from a random point with rate k inside a box of length L and then allowed to grow with a growth velocity=drift velocity
upon collision the colliding branch is instanteneously deleted.
Copyright @ Sabyasachi Sutradhar
*/
double xb,yb,xt,yt,xbj,ybj,xtj,ytj;
int newtip;

double pbc_direction(double direction,double xin,double yin){
  double dir;
  if(xin<=L_Box/2.0 && xin>=-L_Box/2.0 && yin<=L_Box/2.0 && yin>=-L_Box/2.0 ){
    dir=direction;
  }else{
    dir=pi+direction;
  }
  return dir;
}

double pbc(double xin){
  if(xin>L_Box/2.0){xin=xin-L_Box;}
  if(xin<-L_Box/2.0){xin=L_Box+xin;}
  return(xin);
}

void Create_NewTip(int ind){
  Xbase[ind]=L_Box*ran2(&idum)-L_Box/2.0;
  Ybase[ind]=L_Box*ran2(&idum)-L_Box/2.0;
  Dendrite_Length[ind]=0.0;
  Dendrite_Direction[ind]=2.0*pi*ran2(&idum);
  Xtip[ind]=Xbase[ind];
  Ytip[ind]=Ybase[ind];
  Deletion_index[ind]=0;
  state[ind]=1;
  v[ind]=VG;
}

void initialize(){
for(unsigned int i=1;i<=N_Tip;i++){
Create_NewTip(i);
}
}

void ThreeStateDynamics(int i){
  double r1,rt,pt;
        if(state[i]==1){////Growing phase dynamics
          rt=Kgs+Kgp;//total rate
          r1=Kgs/rt;//
          pt = 1.0-exp(-rt*Dt);
          if(ran2(&idum)<pt){
            if (ran2(&idum)<r1){///G->S transition
              state[i]=-1;
              v[i]=-VS;//1.0*logndev(Vs_mean,Vs_sig);//assign shrinking velocity
            }else{
              state[i]=0;//G->P transition
              v[i]=0.0;//gaussdev(0.0,Vp_sig);//assign paused velocity
            }
          }
        }else if(state[i]==-1){///shrinking state dynamics
          rt=Ksg+Ksp;
          r1=Ksg/rt;
          pt = 1.0-exp(-rt*Dt);
          if(ran2(&idum)<pt){
            if (ran2(&idum)<r1){///S->G transition
              state[i]=1;
              v[i]=VG;//logndev(Vg_mean,Vg_sig);//assign growth velocity
            }else{
              state[i]=0;//S->P
              v[i]=0.0;//gaussdev(0.0,Vp_sig);//assign paused velocity
            }
          }
        }else if(state[i]==0){
          rt=Kpg+Kps;
          r1=Kpg/rt;
          pt = 1.0-exp(-rt*Dt);
          if(ran2(&idum)<pt){
            if (ran2(&idum)<r1){///P->G
              state[i]=1;
              v[i]=VG;//logndev(Vg_mean,Vg_sig);//assign growth velocity
            }else{
              state[i]=-1;//P->S
              v[i]=-VS;//1.0*logndev(Vs_mean,Vs_sig);//assign shrinking velocity
            }
          }
      }
////////////////////////////////////////////////
    Dendrite_Length[i]+=v[i]*Dt;///////Grow tips
///////////////////////  Rescue //////////////////////////
    if(Dendrite_Length[i]<0.0){
      if(ran2(&idum)<RescueAfterDeath){
      Dendrite_Length[i]=0.0;
      state[i]=1;
      v[i]=VG;//assign growth velocity
    }else{
      Deletion_index[i]=1;
    }
  }
}


int InsideBox(double x,double y){
  int ind=0;
  if(x<=L_Box/2.0 && x>=-L_Box/2.0){
    if(y<=L_Box/2.0 && y>=-L_Box/2.0){
      ind=1;
    }
  }
  return(ind);
}

double dirc,dl;
double xtipold,ytipold,w,xm,ym;
void Grow_Tips(int i){
if(Deletion_index[i]==0){
ThreeStateDynamics(i);///////Grow tips
Xtip[i]=Xbase[i]+Dendrite_Length[i]*cos(Dendrite_Direction[i]);
Ytip[i]=Ybase[i]+Dendrite_Length[i]*sin(Dendrite_Direction[i]);
///////////////////////////////////////////  inside box
if(InsideBox(Xtip[i],Ytip[i])==1){
for(unsigned int j=1;j<=N_Tip;j++){
  if(j!=i && Deletion_index[j]==0){
    if(lineseg_intersection(Xbase[i],Ybase[i],Xtip[i],Ytip[i],Xbase[j],Ybase[j],Xtip[j],Ytip[j])==1){
      Deletion_index[i]=1;
      Dendrite_Length[i]=0.0;
      CollisionCount++;
      break;
        }
      }
    }
}else{
/////////////////////////////////////////// periodic boundary condition
xt=pbc(Xtip[i]);yt=pbc(Ytip[i]);
xb=xt+Dendrite_Length[i]*cos(pi+Dendrite_Direction[i]);
yb=yt+Dendrite_Length[i]*sin(pi+Dendrite_Direction[i]);
///////////////////////////////////////////////////////////
for(unsigned int j=1;j<=N_Tip;j++){
  if(j!=i && Deletion_index[j]==0){
    if(lineseg_intersection(xb,yb,xt,yt,Xbase[j],Ybase[j],Xtip[j],Ytip[j])==1 || lineseg_intersection(Xbase[i],Ybase[i],Xtip[i],Ytip[i],Xbase[j],Ybase[j],Xtip[j],Ytip[j])==1){
      Deletion_index[i]=1;
      Dendrite_Length[i]=0.0;
      CollisionCount++;
      break;
          }
        }
      }
    }
  }
}


void tip_dynamics(){
newtip=N_Tip;
CollisionCount=0;
for(unsigned int i=1;i<=N_Tip;i++){
if(ran2(&idum)<1.0-exp(-Dendrite_Length[i]*Dt*Branching_Rate)){
  newtip++;
  Create_NewTip(newtip);
}
Grow_Tips(i);
}
N_Tip=newtip;

}

void Remove_DeletedDendrites(){
newtip=0;
Total_Length=0.0;
NG=0;
NP=0;
NS=0;
  for(unsigned int i=1;i<=N_Tip;i++){
  if(Deletion_index[i]==0){
    newtip++;
    Xbase[newtip]=Xbase[i];
    Ybase[newtip]=Ybase[i];
    Xtip[newtip]=Xtip[i];
    Ytip[newtip]=Ytip[i];
    Dendrite_Length[newtip]=Dendrite_Length[i];
    Dendrite_Direction[newtip]=Dendrite_Direction[i];
    Deletion_index[newtip]=Deletion_index[i];
    state[newtip]=state[i];
    v[newtip]=v[i];
    Total_Length+=Dendrite_Length[i];
    if(state[newtip]==1){NG++;}
    if(state[newtip]==0){NP++;}
    if(state[newtip]==-1){NS++;}
    }
  }

  N_Tip=newtip;
}

void print_data(double time){
fp1=fopen(file1,"a");
fprintf(fp1, "%f %d %f %f %f %f %d %f %d %d %d\n",time,N_Tip,Total_Length,N_Tip/Area,Total_Length/Area,Total_Length/N_Tip,CollisionCountMin,(double)CollisionCountMin/Area,NG,NP,NS);
fclose(fp1);
CollisionCountMin=0;
}


void print_conf(int sample,double time,double RescueAfterDeath){
sprintf(fileconf,"Conf-Rescue-%.2f_Sample-%d-Time%.2f.dat",RescueAfterDeath,sample,time);
fconf=fopen(fileconf,"w");
fprintf(fconf,"# Configuration file for time=%g min\n# All length values are in micron\n",time);
fprintf(fconf,"# Simulated in a box of size=%f um\n",L_Box);
fprintf(fconf,"#Parameters:Kb=%f;KGP=%f;KGS=%f;KPG=%f;KPS=%f;KSG=%f;KSP=%f;VG=%f,VS=%f,Beta=%f\n",Branching_Rate,Kgp,Kgs,Kpg,Kps,Ksg,Ksp,VG,VS,RescueAfterDeath);
fprintf(fconf,"# All values are in microns and minutes\n");
fprintf(fconf,"#x_base y_base dx dy Dendrite_Length\n");
for(unsigned int i=1;i<=N_Tip;i++){
if(Deletion_index[i]==0){
fprintf(fconf,"%f %f %f %f %f\n",Xbase[i],Ybase[i],Xtip[i]-Xbase[i],Ytip[i]-Ybase[i],Dendrite_Length[i]);
}
}
fclose(fconf);
}
