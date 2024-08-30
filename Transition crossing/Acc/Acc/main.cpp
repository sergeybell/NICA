#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define PI 3.1415926535897932385
#define M_max 10    // maximum m for FFT (see below)
#define OutPerDstr 302//Period for printing out the particle distr. and force from Long impedance
#define AMPL_H 1.0 // Factor to scale a number of active harmonics in the impedance
#define TAU_D  10. // Time constand for high pass filter of dipole damper
#define TAU_Q  2000.  // Time constand for high pass filter of quadrupole damper
#define V_SATUR  0.05 // Voltage suturation for quadrupole damper
//========================================================
/*
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform
*/
short FFT(short int dir,long m,double *x,double *y)
{
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++)
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return(true);
}
//========================================================
int main (void){
FILE *fp, *fp1;
char *fnameCP="/Users/Bell 1/GIT_REPOS/NICA/Transition crossing/Acceleration/Acceleration/AccCycleParam.dat";//file name for particle coordinates: phi, DeltaP_P, ro
char *fname=(char *)"/Users/Bell 1/GIT_REPOS/NICA/Transition crossing/Acceleration/Acceleration/PartCoord.dat";      // file name for particle coordinates: phi, DeltaP_P, ro
char *fnameZ=(char *)"/Users/Bell 1/GIT_REPOS/NICA/Transition crossing/Acceleration/Acceleration/impedance.dat";     // file name for effective impedances
char *fnameOut=(char *)"PartCoordOut.dat";  // file name for final particle coordinates
char *fnameCB=(char *)"/Users/Bell 1/GIT_REPOS/NICA/Transition crossing/Acceleration/Acceleration/AccCycleBunchPrm.dat";//file name for resulting turn-by-turn bunch parameters
char buf[256];
int Np, Nt, Nh, Nhh, n, i, k, j, m, nc=0, nh_max, last_turn;
double capa, capa2, Vd, Vdd, Vd_a, x, y, x_p, s, Vb, dphi;
double *phi, *Dpp, *ro;
double *Ekin, *Vrf, *PhiAcc;
double *g_d, *g_q, *DeltaP, *DeltaPhiD, *DeltaVQ, S_DVQ, dphi_tot; //feedbacks
double *eta_p, *p_accept; //non-linear slip factor, and momentum acceptance
double *eta, *gamma, *beta, *PhiCor;
double *Ntot, *PAv, *PhiAv, *SigmaPhi, *SigmaP, *P_Phi;
double *Ref, *Imf, h;
double Asc,*ReZ, *ImZ, *ReZsc, *ImZsc, ReZt, ImZt;
time_t t1, t2;
    time(&t1);
//
    // Booster parameters
   double alpha=0.03346;
   double q=84;
   double Mp=938.256E6;
   //
    // reading beam energy and accelerating voltage through the cycle.
   fp=fopen(fnameCP,"r");
   Nt=0;     // number of turns
   while(fgets(buf, 255, fp)) Nt++;
   fclose(fp);
   Ekin=new double [Nt];
   Vrf=new double [Nt];
   PhiAcc=new double [Nt];
   PhiCor=new double [Nt];
   g_d=new double [Nt];
   g_q=new double [Nt];
   DeltaP=new double [Nt];
   eta_p=new double [Nt];
   p_accept=new double [Nt];
   eta=new double [Nt];
   fp=fopen(fnameCP,"r");
   for(n=0; n<Nt; n++){
       if(fgets(buf, 255, fp)==NULL)break;
      sscanf(buf,"%le %le %le %le %le %le %le %le %le %le",
                  &Ekin[n],&Vrf[n],&PhiAcc[n],&PhiCor[n],&g_d[n],&g_q[n],&DeltaP[n],&eta_p[n],&p_accept[n], &eta[n]);
   }
   fclose(fp);
    // computing beam parameters for acceleration
   gamma=new double [Nt];
   beta=new double [Nt];
   for(n=0; n<Nt; n++){
       gamma[n]=1.+Ekin[n]/Mp;
      beta[n]=sqrt( 1.-1./(gamma[n]*gamma[n]) );}
   // Arrays for tracking reporting
   Ntot=new double [Nt];
   PhiAv=new double [Nt];
   PAv=new double [Nt];
   SigmaPhi=new double [Nt];
   SigmaP=new double [Nt];
   P_Phi=new double [Nt];
   // Arrays and variasbles for corrections coming from dampers
   S_DVQ=0;
   DeltaPhiD=new double [Nt];
   DeltaVQ=new double [Nt];
   //
    // Reading particle coordinates
   fp=fopen(fname,"r");
   Np=0;    // number of particles
   while(fgets(buf, 255, fp)) Np++;
   fclose(fp);
   phi=new double [Np];
   Dpp=new double [Np];
   ro=new double [Np];   // particle wieght
   fp=fopen(fname,"r");
   for(n=0; n<Np; n++){
       if(fgets(buf, 255, fp)==NULL)break;
      sscanf(buf,"%le %le %le", &phi[n],&Dpp[n], &ro[n]);
   }
   fclose(fp);
   //
    // Reading effective impedances
   fp=fopen(fnameZ,"r");
   Nh=0; // number of harmonics
   while(fgets(buf, 255, fp)) Nh++;
   fclose(fp);
   m=1; k=2;
   for(m=1; m<M_max; m++){
       if((Nh-1)==k)break;
      k=2*k;}
   if(m==M_max){
       printf("Incorect array length for the impedance array(L=%d)\n", Nh);
      exit(1);}
   m=m+1;
   Nhh=(Nh-1)*2;
   h=2.*PI/Nhh;
    ReZ = new double[Nh];
    ImZ = new double[Nh];
   ReZsc=new double[Nh];
   ImZsc=new double[Nh];
   fp=fopen(fnameZ,"r");
   for(n=0; n<Nh; n++){
       if(fgets(buf, 255, fp)==NULL)break;
      sscanf(buf,"%le %le %le %le", &ReZ[n],&ImZ[n], &ReZsc[n], &ImZsc[n]);}
   fclose(fp);
   Ref=new double [Nhh]; //Array for longitudinal distribution function (Re)
   Imf=new double [Nhh]; //Array for longitudinal distribution function (Im)
   //
    // Turn-by-turn tracking
   last_turn=1;
    fp1=fopen("distr.dat", "w");
    printf("Tracking starts\n\tNumber of particles=%d, Number of turns=%d\n",Np, Nt);
   for(n=0; n<Nt-2; n++){
      // computation of avarage beam parameters
      Ntot[n]=0;
      PhiAv[n]=0;
         PAv[n]=0;
       for(i=0; i<Np; i++){
         Ntot[n]=Ntot[n]+ro[i];
         PhiAv[n]=PhiAv[n]+phi[i]*ro[i];
         PAv[n]=PAv[n]+Dpp[i]*ro[i];}
      if(Ntot[n]<2.){
          last_turn=0;
         break;}
      PhiAv[n]=PhiAv[n]/Ntot[n];
      PAv[n]=PAv[n]/Ntot[n];
      SigmaPhi[n]=0;
      P_Phi[n]=0;
      SigmaP[n]=0;
       for(i=0; i<Np; i++){
          x=phi[i]-PhiAv[n];
         SigmaPhi[n]=SigmaPhi[n]+x*x*ro[i];
          y=Dpp[i]-PAv[n];
         SigmaP[n]=SigmaP[n]+y*y*ro[i];
         P_Phi[n]=P_Phi[n]+x*y*ro[i];}
      SigmaPhi[n]=sqrt(SigmaPhi[n]/Ntot[n]);
      SigmaP[n]=sqrt(SigmaP[n]/Ntot[n]);
      P_Phi[n]=P_Phi[n]/Ntot[n];
      // Tracking step for particle longitudinal coordinate (phi)
      capa=beta[n+1]/beta[n];
      x=2.*PI*q*eta[n];
      x_p=2.*PI*q*eta_p[n];
      dphi=capa*PhiAcc[n]+PhiCor[n]-PhiAcc[n+1]-PhiCor[n+1];
           for(i=0; i<Np; i++)if(ro[i]>0.){
          phi[i]=capa*phi[i]-Dpp[i]*(x+Dpp[i]*x_p)+dphi;
         if(eta[n]<0.){
             if(phi[i]>PI)phi[i]=phi[i]-2*PI;
             if(phi[i]<-PI)phi[i]=phi[i]+2*PI;}
         else{
             if(phi[i]>2.*PI)phi[i]=phi[i]-2*PI;
             if(phi[i]<0.)phi[i]=phi[i]+2*PI;}
      }
      // Computation of feedback parameters
      if(n>0){
         S_DVQ=S_DVQ+DeltaVQ[n-1];
          DeltaPhiD[n]=g_d[n]*(PAv[n]-DeltaP[n])*beta[n]*beta[n]*gamma[n]*Mp/Vrf[n];
          DeltaVQ[n]=g_q[n]*(SigmaPhi[n]-SigmaPhi[n-1])-S_DVQ/TAU_Q;
         if(DeltaVQ[n]>V_SATUR)DeltaVQ[n]=V_SATUR;
         if(DeltaVQ[n]<-V_SATUR)DeltaVQ[n]=-V_SATUR;}
      else DeltaPhiD[n]=DeltaVQ[n]=0;
      // Computation particle longituidinal distribution (histogram)
      for(k=0; k<Nhh; k++){Ref[k]=0; Imf[k]=0;}
       for(i=0; i<Np; i++){
          if(ro[i]>0.){
            if(fabs(Dpp[i])>p_accept[n]){ro[i]=0.; continue;}
              j=floor((phi[i]+PI)/h);
            x=(phi[i]+PI-h*j)/h;
            while(j<0)j=j+Nhh;
            while(j>=Nhh)j=j-Nhh;
             Ref[j]=Ref[j]+(1.-x)*ro[i];
            if((j+1)==Nhh)Ref[j+1-Nhh]=Ref[j+1-Nhh]+x*ro[i];
            else Ref[j+1]=Ref[j+1]+x*ro[i];
          }
      }
        if(nc==0){
          fprintf(fp1, "%d", n);
          for(k=0; k<Nhh; k++) fprintf(fp1, " \t%le", Ref[k]);
            fprintf(fp1, "\n");
        }
      // Computation of impedance induced force
      FFT(1, m, Ref, Imf);
      Asc=beta[0]*gamma[0]*gamma[0] / (beta[n]*gamma[n]*gamma[n]);
      nh_max=(int)(AMPL_H*Nh/(16.*SigmaPhi[n]));
      if(nh_max<5)nh_max=5;
      if(nh_max>Nh)nh_max=Nh;
      for(k=0; k<nh_max; k++){
          ReZt=ReZ[k]+Asc*ReZsc[k];
            ImZt=-(ImZ[k]+Asc*ImZsc[k]);
         x=ReZt*Ref[k]-ImZt*Imf[k];
         Imf[k]=-(ReZt*Imf[k]+ImZt*Ref[k]);
         Ref[k]=-x;
      }
      if(nh_max<(Nh-1))for(k=nh_max; k<Nh; k++)Ref[k]=Imf[k]=0.;
      for(k=1; k<Nh-1; k++){
         Ref[Nhh-k]=Ref[k];
         Imf[Nhh-k]=-Imf[k];
      }
      FFT(-1, m, Ref, Imf);
        if(nc==0){
          fprintf(fp1, "%d ", n);
          for(k=0; k<Nhh; k++) fprintf(fp1, " \t%le", Ref[k]);
            fprintf(fp1, "\n");
         nc=OutPerDstr;
        }
      else {nc--;}
      // Tracking step over particle momentum (dpp)
      capa2=1.+alpha*(beta[n+1]*gamma[n+1]/(beta[n]*gamma[n])-1.);
      capa=capa2*(beta[n]*beta[n]*gamma[n] / (beta[n+1]*beta[n+1]*gamma[n+1]));
      Vdd=1./(Mp*beta[n+1]*beta[n+1]*gamma[n+1]);
      Vd=Vrf[n]*Vdd;
      s=sin(PhiAcc[n+1]);
      Vd_a=Vd*(1.+DeltaVQ[n]-DeltaPhiD[n]*tan(PhiAcc[n+1]));
      dphi_tot=DeltaPhiD[n]+DeltaVQ[n]*tan(PhiAcc[n+1]);
          for(i=0; i<Np; i++)if(ro[i]>0.){
            j=floor((phi[i]+PI)/h);
         x=(phi[i]+PI-h*j)/h;
         while(j<0)j=j+Nhh;
         while(j>=Nhh)j=j-Nhh;
         if((j+1)==Nhh)Vb=Ref[j+1-Nhh];
         else Vb=Ref[j+1];
         Vb=(1.-x)*Ref[j] + x*Vb;
         Dpp[i]=capa*Dpp[i]-Vd_a*sin(phi[i]+dphi_tot)-Vd*s+Vb*Vdd;
         phi[i]=phi[i]/capa2;
      }
   }
   // computation of beam parameters for the last turn
   n=Nt-2;
    if(last_turn){
       Ntot[n]=0;
       PhiAv[n]=0;
       PAv[n]=0;
       for(i=0; i<Np; i++){
           Ntot[n]=Ntot[n]+ro[i];
          PhiAv[n]=PhiAv[n]+phi[i]*ro[i];
          PAv[n]=PAv[n]+Dpp[i]*ro[i];}
      if(Ntot[n]>0){
           PhiAv[n]=PhiAv[n]/Ntot[n];
           PAv[n]=PAv[n]/Ntot[n];
           SigmaPhi[n]=0;
           SigmaP[n]=0;
           P_Phi[n]=0;
           for(i=0; i<Np; i++){
               x=phi[i]-PhiAv[n];
              SigmaPhi[n]=SigmaPhi[n]+x*x*ro[i];
              y=Dpp[i]-PAv[n];
              SigmaP[n]=SigmaP[n]+y*y*ro[i];
              P_Phi[n]=P_Phi[n]+x*y*ro[i];}
           SigmaPhi[n]=sqrt(SigmaPhi[n]/Ntot[n]);
           SigmaP[n]=sqrt(SigmaP[n]/Ntot[n]);
           P_Phi[n]=P_Phi[n]/Ntot[n];
      }
      else n=n-1;
   }
   //
   // writing final particle coordinates
   fp=fopen(fnameOut,"w");
   for(k=0; k<Np; k++){
      fprintf(fp,"%le\t %le\t %le\t %le\t %le\n",phi[k],Dpp[k],ro[k]);
   }
    // writing average beam parameters
   if(n!=Nt-2)Nt=n;
   fp=fopen(fnameCB,"w");
   for(n=0; n<Nt-1; n++){
      fprintf(fp,"%le\t %le\t %le\t %le\t %le\t %le\t %le\t %le\t %le\t %le\n",
      Ntot[n],PhiAv[n],PAv[n],SigmaPhi[n],SigmaP[n],Vrf[n],PhiAcc[n],DeltaPhiD[n],DeltaVQ[n],P_Phi[n]);
   }
    time(&t2);
    printf("Done in %d sec.\n", t2-t1);
}


/*
for(k=0; k<Nhh; k++){
      printf("k=%d\t %le\t %le\t \n", k, Ref[k],Imf[k]);
}
fp=fopen(fnameOut,"w");
for(k=0; k<Nhh; k++){
      fprintf(fp,"%d\t %le\t %le\t \n", k, Ref[k],Imf[k]);
      printf("%d\t %le\t %le\t \n", k, Ref[k],Imf[k]);
}
exit (1);
*/
