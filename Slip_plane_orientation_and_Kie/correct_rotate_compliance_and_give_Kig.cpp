#include <iostream>
#include <cmath>
#include<stdlib.h>

// The code is for PhD research purposes

// Author: Praveenkumar Hiremath
// Email: praveenkumar.hiremath@mek.lth.se (Email at the University)
//       praveenkumar.hiremath2911@gmail.com (Private email)


//Function to calculate compliance constants in rotated system with plane strain condition
//void DisplaceAtoms::rotation(double o,double p,double q,double u,double v,double w,double s11,double s12,double s44,double s11_New,double s12_New,double s22_New,double s16_New,double s26_New,double s66_New,double A[])
//{
using namespace std;
int main()
{
double o,p,q,u,v,w;
double s11,s12,s13,s14,s15,s16;
double s21,s22,s23,s24,s25,s26;
double s31,s32,s33,s34,s35,s36;
double s41,s42,s43,s44,s45,s46;
double s51,s52,s53,s54,s55,s56;
double s61,s62,s63,s64,s65,s66;
double S[36];

double S_New[3][3][3][3]={0},T[3][3]={0},S_Old[3][3][3][3]={0};
double a,b,c;
int g,h,m,n,s,t,k,l;

//Crack propagation direction
o=1; p=0;  q=0; 

//Crack plane direction
u=0; v=0 ; w=4; 


//compliance constants in original system
//double s22=s11,s33=s11;
//double s13=s12,s21=s12,s23=s12,s31=s12,s32=s12;
//double s55=s44,s66=s44;
//double s14=0,s15=0,s16=0,s24=0,s25=0,s26=0,s34=0,s35=0,s36=0,s41=0,s42=0,s43=0,s45=0,s46=0;
//double s51=0,s52=0,s53=0,s54=0,s56=0,s61=0,s62=0,s63=0,s64=0,s65=0;

//compliance constants in original system
n=36;
FILE *input;
input=fopen("compliance.dat","r");
int i;
for(i=0;i<n;i++)
 {
  //fscanf(input, "%lf \t %lf \n",&a[i],&b[i]);           
  fscanf(input, "%lf \n",&S[i]);                                            
  //printf("Sij= %lf \n",S[i]);
}

s11=S[0]; s12=S[1]; s13=S[2]; s14=S[3]; s15=S[4]; s16=S[5]; 
s21=S[6]; s22=S[7]; s23=S[8]; s24=S[9]; s25=S[10]; s26=S[11];
s31=S[12]; s32=S[13]; s33=S[14]; s34=S[15]; s35=S[16]; s36=S[17];
s41=S[18]; s42=S[19]; s43=S[20]; s44=S[21]; s45=S[22]; s46=S[23];
s51=S[24]; s52=S[25]; s53=S[26]; s54=S[27]; s55=S[28]; s56=S[29];
s61=S[30]; s62=S[31]; s63=S[32]; s64=S[33]; s65=S[34]; s66=S[35];

//Building 4th order compliance tensor in original system --> (*) condition in Report.
S_Old[0][0][0][0]=s11;
S_Old[0][0][1][1]=s12;
S_Old[0][0][2][2]=s13;
S_Old[0][0][0][1]=s16/2;     S_Old[0][0][1][0]=s16/2;
S_Old[0][1][0][1]=s66/4;     S_Old[0][1][1][0]=s66/4;    S_Old[1][0][0][1]=s66/4;    S_Old[1][0][1][0]=s66/4;
S_Old[0][1][1][1]=s62/2;     S_Old[1][0][1][1]=s62/2;
S_Old[0][1][2][2]=s63/2;     S_Old[1][0][2][2]=s63/2;
S_Old[1][1][1][1]=s22;
S_Old[1][1][2][2]=s23;
S_Old[1][1][0][1]=s26/2;     S_Old[1][1][1][0]=s26/2;
S_Old[2][2][2][2]=s33;
S_Old[1][1][0][0]=s21;
S_Old[2][2][0][0]=s31;
S_Old[2][2][1][1]=s32;
S_Old[1][2][1][2]=s44/4;     S_Old[2][1][2][1]=s44/4;    S_Old[1][2][2][1]=s44/4;    S_Old[2][1][1][2]=s44/4;
S_Old[0][2][0][2]=s55/4;     S_Old[2][0][2][0]=s55/4;    S_Old[0][2][2][0]=s55/4;    S_Old[2][0][0][2]=s55/4;



//to find the third axis(crack front) of the new coordinate system-> Cross Product:
//a=(p*w)-(q*v);  b=(u*q)-(w*o);  c=(o*v)-(p*u);
a=(p*w)-(q*v);  b=(u*q)-(w*o);  c=(o*v)-(p*u);
cout<<"a="<<a<<"\t"<<"b="<<b<<"\t"<<"c="<<c<<"\n";

//Normalization of the above vectors to form Orthonormal basis set:
double X1=o/sqrt(pow(o,2)+pow(p,2)+pow(q,2)); double Y1=p/sqrt(pow(o,2)+pow(p,2)+pow(q,2)); double Z1=q/sqrt(pow(o,2)+pow(p,2)+pow(q,2));
double X3=a/sqrt(pow(a,2)+pow(b,2)+pow(c,2)); double Y3=b/sqrt(pow(a,2)+pow(b,2)+pow(c,2)); double Z3=c/sqrt(pow(a,2)+pow(b,2)+pow(c,2));
double X2=u/sqrt(pow(u,2)+pow(v,2)+pow(w,2)); double Y2=v/sqrt(pow(u,2)+pow(v,2)+pow(w,2)); double Z2=w/sqrt(pow(u,2)+pow(v,2)+pow(w,2));

//Rotation matrix: T(ij)=x'(i).x(j) Here i is associated with the rotated system axes and j with original.
T[0][0]=X1;    T[0][1]=Y1;     T[0][2]=Z1;
T[1][0]=X2;    T[1][1]=Y2;     T[1][2]=Z2;
T[2][0]=X3;    T[2][1]=Y3;     T[2][2]=Z3;


s=0; t=0; k=0; l=0;
//compliance constants in rotated coordinate system
for(s=0;s<3;s++)
{ 
 for(t=0;t<3;t++)
 {
  for(k=0;k<3;k++)
  {
   for(l=0;l<3;l++)
   {
for(g=0;g<3;g++)
{ 
 for(h=0;h<3;h++)
 {
  for(m=0;m<3;m++)
  {
   for(n=0;n<3;n++)
   {

  
S_New[s][t][k][l]+= T[s][g]*T[t][h]*S_Old[g][h][m][n]*T[k][m]*T[l][n];  //Rotation


   }
  }
 }
}
   }
  }
 }
}


double sN11,sN12,sN13,sN14,sN15,sN16;
double sN21,sN22,sN23,sN24,sN25,sN26;
double sN31,sN32,sN33,sN34,sN35,sN36;
double sN41,sN42,sN43,sN44,sN45,sN46;
double sN51,sN52,sN53,sN54,sN55,sN56;
double sN61,sN62,sN63,sN64,sN65,sN66;

//No plane stress and No plane strain   (*) condition in Report.
sN11=S_New[0][0][0][0]; sN12=S_New[0][0][1][1]; sN13=S_New[0][0][2][2]; sN14=2*S_New[0][0][1][2]; sN15=2*S_New[0][0][0][2]; sN16=2*S_New[0][0][0][1];
sN21=S_New[1][1][0][0]; sN22=S_New[1][1][1][1]; sN23=S_New[1][1][2][2]; sN24=2*S_New[1][1][1][2]; sN25=2*S_New[1][1][0][2]; sN26=2*S_New[1][1][0][1];
sN31=S_New[2][2][0][0]; sN32=S_New[2][2][1][1]; sN33=S_New[2][2][2][2]; sN34=2*S_New[2][2][1][2]; sN35=2*S_New[2][2][0][2]; sN36=2*S_New[2][2][0][1];
sN41=2*S_New[1][2][0][0]; sN42=2*S_New[1][2][1][1]; sN43=2*S_New[1][2][2][2]; sN44=4*S_New[1][2][1][2]; sN45=4*S_New[1][2][0][2]; sN46=4*S_New[1][2][0][1];
sN51=2*S_New[0][2][0][0]; sN52=2*S_New[0][2][1][1]; sN53=2*S_New[0][2][2][2]; sN54=4*S_New[0][2][1][2]; sN55=4*S_New[0][2][0][2]; sN56=4*S_New[0][2][0][1];
sN61=2*S_New[0][1][0][0]; sN62=2*S_New[0][1][1][1]; sN63=2*S_New[0][1][2][2]; sN64=4*S_New[0][1][1][2]; sN65=4*S_New[0][1][0][2]; sN66=4*S_New[0][1][0][1];


//after rotation compliance constants to array
double Final_S[6][6]={0};
Final_S[0][0]=sN11; Final_S[0][1]=sN12; Final_S[0][2]=sN13; Final_S[0][3]=sN14; Final_S[0][4]=sN15; Final_S[0][5]=sN16; 
Final_S[1][0]=sN21; Final_S[1][1]=sN22; Final_S[1][2]=sN23; Final_S[1][3]=sN24; Final_S[1][4]=sN25; Final_S[1][5]=sN26; 
Final_S[2][0]=sN31; Final_S[2][1]=sN32; Final_S[2][2]=sN33; Final_S[2][3]=sN34; Final_S[2][4]=sN35; Final_S[2][5]=sN36; 
Final_S[3][0]=sN41; Final_S[3][1]=sN42; Final_S[3][2]=sN43; Final_S[3][3]=sN44; Final_S[3][4]=sN45; Final_S[3][5]=sN46; 
Final_S[4][0]=sN51; Final_S[4][1]=sN52; Final_S[4][2]=sN53; Final_S[4][3]=sN54; Final_S[4][4]=sN55; Final_S[4][5]=sN56; 
Final_S[5][0]=sN61; Final_S[5][1]=sN62; Final_S[5][2]=sN63; Final_S[5][3]=sN64; Final_S[5][4]=sN65; Final_S[5][5]=sN66; 


FILE *f = fopen("rotated_compliance.dat", "wb");
fwrite(Final_S, sizeof(double), sizeof(Final_S), f);
fclose(f);

double b11,b22,b12,b66;

//b11=((s11*s33)-(s13*s13))/s33;
//b22=((s22*s33)-(s23*s23))/s33;
//b12=((s12*s33)-(s13*s23))/s33;
//b66=((s66*s33)-(s26*s26))/s33;

b11=((sN11*sN33)-(sN13*sN13))/sN33;
b22=((sN22*sN33)-(sN23*sN23))/sN33;
b12=((sN12*sN33)-(sN13*sN23))/sN33;
b66=((sN66*sN33)-(sN26*sN26))/sN33;

double temp_1,temp_2, temp_3, temp_4;
temp_1=(b11*b22*0.5);
temp_2=sqrt(b22/b11);
temp_3=((2*b12)+b66);
temp_4=(2*b11);

double B_1=sqrt((temp_1)*(temp_2+(temp_3/temp_4)));

double gamma_s=3.241;
double G_G=2*gamma_s;
double K_G=sqrt(G_G/B_1);

double K_G_N=sqrt(G_G/B_1);
double K_G_N_MPa_root_m=K_G_N/1e6;

//print "Emit K (Pa.m^0.5) = ",K_G_N
printf("Emit K (Pa.m^0.5) = %lf \n",K_G_N);
//print "Emit K (MPa.m^0.5) = ",K_G_N_MPa_root_m
printf("Emit K (Pa.m^0.5) = %lf \n",K_G_N_MPa_root_m);

return 0;
}

