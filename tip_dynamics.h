///////////////////////////////////////////////////////////////////////////
//////////// add length to a growing tip //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void add_tiplength(int i,double l_add){
int brcol,somacol,boundarycol;
double dis,dtheta,dir;
double lold=Tiplength[i];
int n_old=Tipl[i];
double l_end=Tiplength[i]-floor(round(100000*Tiplength[i]/Tip_intv)/100000)*Tip_intv;
double l_totaladd=l_add+l_end;
int n_add=(int) floor(round(100000*l_totaladd/Tip_intv)/100000);
int ne;
if(l_end<=0.000001){ne=n_old;}
else{ne=n_old-1;}
dtheta=gaussdev(0.0,sqrt(2.0*l_totaladd/Tip_Persis));//angle at which the new length is added
dir=wrap_angle(assign_direction(i,ne)+dtheta);
double xt,yt;
xt=Tipx[i][n_old]+l_totaladd*cos(dir);////new tip positions
yt=Tipy[i][n_old]+l_totaladd*sin(dir);
/////////////////////////////////////////// determine if the new position hits the boundary soma or other branch
brcol=branch_collision(i,xt,yt);
somacol=soma_collision(i,xt,yt);

if(somacol==0 && brcol==0){
  if(n_old+n_add<Max_Tippoints){
    int j=0;
    dis=l_totaladd;
    double disadd;
    while(dis>0.0){
      j++;
      dis-=Tip_intv;
      disadd=Tip_intv;if(dis<0.0){disadd=Tip_intv+dis;}
      Tipx[i][ne+j]=Tipx[i][ne+j-1]+disadd*cos(dir);
      Tipy[i][ne+j]=Tipy[i][ne+j-1]+disadd*sin(dir);
      if(dis<0.0){break;}
    }
    Tipl[i]=ne+j;
  }else{
    int n_tot=Max_Tippoints-1;
    n_old=n_tot-n_add;
    for(int j=1;j<=Tipl[i];j++){tipx[j]=Tipx[i][j];tipy[j]=Tipy[i][j];}
    tippoints=equispace_points(tipx,tipy,Tipl[i],n_old);
    for(int j=1;j<=n_old;j++){Tipx[i][j]=tippoints[j][1];Tipy[i][j]=tippoints[j][2];}
    for(int j=1;j<=n_add;j++){
      Tipx[i][n_old+j]=Tipx[i][n_old+j-1]+Tip_intv*cos(dir);
      Tipy[i][n_old+j]=Tipy[i][n_old+j-1]+Tip_intv*sin(dir);
    }
    Tipl[i]=n_tot;
    n_tot=n_tot+1;
      Tipx[i][n_tot]=Tipx[i][n_old]+l_add*cos(dir);
      Tipy[i][n_tot]=Tipy[i][n_old]+l_add*sin(dir);
      Tipl[i]=n_tot;
  }
}
else{
  state[i]=-1;
  v[i]=-1.5;
  collision[i]=1;
  collision_time[i]=0.0;
}
/*
if(brcol==-1){
  state[i]=0;
  v[i]=0;
  collision[i]=0;
}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Collision with boundary /////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double prob,sig;
if(strcmp(Boundary_type,"repulsive")==0 || strcmp(Boundary_type,"static")==0){
//if(Boundary_Cond==1 || Boundary_Cond==-1){
sig=25.0;
prob = 0.5*(1.0+erf((fabs(yt)-SizeLR-1.0*sig)/(sqrt(2.0)*sig)));
if(prob<0.0){prob=0.0;}
  if(ran2(&idum)<1.0-exp(-prob*Dt)){
    if(strcmp(Boundary_type,"static")==0){state[i]=0;v[i]=0;collision[i]=0;  collision_time[i]=0.0;}
    if(strcmp(Boundary_type,"repulsive")==0){state[i]=-1;v[i]=-1.5;collision[i]=0;  collision_time[i]=0.0;}
  }

  sig=15.0;
  prob = 0.5*(1.0+erf((fabs(xt)-SizeAP-1.0*sig)/(sqrt(2.0)*sig)));
  if(prob<0.0){prob=0.0;}
  if(ran2(&idum)<1.0-exp(-prob*Dt)){
    if(strcmp(Boundary_type,"static")==0){state[i]=0;v[i]=0;collision[i]=0;  collision_time[i]=0.0;}
    if(strcmp(Boundary_type,"repulsive")==0){state[i]=-1;v[i]=-1.5;collision[i]=0;  collision_time[i]=0.0;}
  }
}
if(strcmp(Boundary_type,"hybrid")==0){
  sig=15.0;
  prob = 0.5*(1.0+erf((fabs(xt)-SizeAP-1.0*sig)/(sqrt(2.0)*sig)));
  if(prob<0.0){prob=0.0;}
  if(ran2(&idum)<1.0-exp(-prob*Dt)){
  state[i]=0;v[i]=0;collision[i]=0;  collision_time[i]=0.0;
  }
}


}

//////////////////////////////////////////////////////////////////////////////////
//////////////// delete tiplength if the tip is shrinking ////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void delete_length(int i, double l_del){
  double lold=Tiplength[i];
    l_del=fabs(l_del);
      double dis=l_del;
    int n_del=(int) floor(round(100000*l_del/Tip_intv)/100000);
      int j=Tipl[i];
    if(Tipl[i]>=3){///////delete until tip length is > 0.1 micron
      while(dis>0.0){
        dis-=distance2D(Tipx[i][j],Tipy[i][j],Tipx[i][j-1],Tipy[i][j-1]);
        j--;
      }
      Tipl[i]=j+1;
      if(Tipl[i]>=3){///////delete until tip length is > 0.1 micron
      double th=assign_direction(i,Tipl[i]);
      Tipx[i][Tipl[i]]=Tipx[i][j]+fabs(dis)*cos(th);
      Tipy[i][Tipl[i]]=Tipy[i][j]+fabs(dis)*sin(th);
    }else{
      Tipl[i]=0;Tipage[i]=0.0;
    }
    }else{
        Tipl[i]=0;////delete the tip
        Tipage[i]=0.0;
    }
}
///////////////////////////////////////////////////////////////////////////////////
////////////// Simulation of free tips with the three state (G,P,S) model /////////
///////////////////////////////////////////////////////////////////////////////////
void Tip_Dynamics(){
  int ll,mm,nn=0;
  double r,rt,r1,pt;
  for (int i=1; i<=n_tip; i++){
    if (Tiptype[i]==1){///if the tip is free then do the dynamics
      Collision_Time(i);
      if(Tipl[i]>0 && Tipage[i]>0.0){///three state tip dynamics
        if(state[i]==1){////Growing phase dynamics
          rt=Kgs+Kgp;//total rate
          r1=Kgs/rt;//
          pt = 1.0-exp(-rt*Dt);
          if(ran2(&idum)<pt){
            if (ran2(&idum)<r1){///G->S transition
              state[i]=-1;
              v[i]=-1.0*logndev(Vs_mean,Vs_sig);//assign shrinking velocity
            }else{
              state[i]=0;//G->P transition
              v[i]=gaussdev(0.0,Vp_sig);//assign paused velocity
            }
          }
        }else if(state[i]==-1){///shrinking state dynamics
          rt=Ksg+Ksp;
          r1=Ksg/rt;
          pt = 1.0-exp(-rt*Dt);
          if(ran2(&idum)<pt){
            if (ran2(&idum)<r1){///S->G transition
              state[i]=1;
              v[i]=logndev(Vg_mean,Vg_sig);//assign growth velocity
            }else{
              state[i]=0;//S->P
              v[i]=gaussdev(0.0,Vp_sig);//assign paused velocity
            }
          }
        }else if(state[i]==0){
          rt=Kpg+Kps;
          r1=Kpg/rt;
          pt = 1.0-exp(-rt*Dt);
          if(ran2(&idum)<pt){
            if (ran2(&idum)<r1){///P->G
              state[i]=1;
              v[i]=logndev(Vg_mean,Vg_sig);//assign growth velocity
            }else{
              state[i]=-1;//P->S
              v[i]=-1.0*logndev(Vs_mean,Vs_sig);//assign shrinking velocity
            }
          }
        }

        if(i<=N_TipInitial && Tiplength[i]<R_Soma+1.0){//intial tips regrow from soma
          state[i]=1;v[i]=logndev(Vg_mean,Vg_sig);collision[i]=0;
        }
        ///////////////////  end of current state and velocity///////////////
        if(v[i]>0.0){//velocity is postive; add length to the tip
          add_tiplength(i,v[i]*Dt);
        }else if (v[i]<0.0){//velocity is negative; subtract length from the tip
          delete_length(i,fabs(v[i]*Dt));
        }
        Tiplength[i]=calculate_curvelength(i,Tipl[i]);
      }
    }
  }
}
