
#include <stdio.h>
#include <math.h>
#include <algorithm>

#define pi 3.14159265

void makjonswap(float hm0gew,float fp,float mainang,float rt,float scoeff,float gam,float *theta,int ntheta,float& TTrep,float * &Stt)
{
    

//float hm0gew=2.0f;//Significant wave height (m)
//float fp=1.0f/12.0f; //Wave peak frequency (Hz)
//float mainang=0; //wave mean direction (angle of incidence º)
//float rt=1000; //Boundary duration
//float scoeff=100;// spread coef n.u.
//float gam=3.3f;//: peak enhancement factor, optional parameter (DEFAULT 3.3)

//printf("hs=%f\nfp=%f\nang=%f\nscoeff=%f\ngam=%f\n",hm0gew,fp,mainang,scoeff,gam);	
	
	
//printf("fp=%f\n",fp);
//// 

float fnyq = 3.0f*fp;
float dfj= fp/20.0f;

// 
int nfreq=ceil((fnyq-dfj)/dfj);
float * f;
float * x;
float * y;

f=(float *)malloc(nfreq*sizeof(float));
//f= new float[nfreq];
//x= new float[nfreq];
//y= new float[nfreq];
x=(float *)malloc(nfreq*sizeof(float));
y=(float *)malloc(nfreq*sizeof(float));
//printf("Hello world!");
float xa,ymax,ysum;
float sigma,fac1,fac2,fac3;

ysum=0.0f;
ymax=0.0f;

for (int i=0; i<nfreq ;i++)
{
    f[i]=(i+1)*dfj;
    
    // x: nondimensional frequency, divided by the peak frequency
    x[i]=f[i]/fp;
    xa=abs(x[i]);
    if (xa==0)
    {
        xa=1e-7f;
    }
    if (xa<1)
    {
        sigma=0.07f;
    }
    else
    {
        sigma=0.09f;
    }

    fac1=pow(xa,-5);
    fac2=exp(-1.25*(pow(xa,-4.0f)));


    fac3=pow(gam,(exp(-(pow(xa-1.0f,2)/(2*pow(sigma,2))))));

    y[i]=fac1*fac2*fac3;
	//printf("fac1=%f\nfac2=%f\nfac3=%f\n",fac1,fac2,fac3);
    ysum=ysum+y[i];
    if (y[i]>=ymax)
    {
        ymax=y[i];
    }

}
ysum=ysum/ymax;
for (int i=0; i<nfreq;i++)
{
    y[i]=y[i]/ymax;
    
    y[i]=y[i]*pow((hm0gew/(4*sqrt(ysum*dfj))),2);
    //printf("y[%d]=%f\n",i,y[i]);
}


////
        
//float t1=-(pi)/2;
//int ntheta=101;
//float *theta;

float * Dd;
Dd= (float *)malloc(ntheta*sizeof(float));
//Dd= new float[ntheta];
//theta= (float *)malloc(ntheta*sizeof(float));//new float[ntheta];
float ddsum=0.0f;
//theta=(0:100)*((pi)/100)+t1;

for(int i=0; i<ntheta; i++)
{
    //theta[i]=i*((pi)/100.0f)+t1;
    Dd[i] = pow(cos((theta[i]-mainang)/2.0f),2.0f*scoeff);
    ddsum=ddsum+Dd[i];
    //printf("theta[%d]=%f\n",i,theta[i]);
}

float dang=theta[1]-theta[0];

//mainang=(1.5d0*p1-alfa)-mainang*atan(1.d0)/45.0d0;

for(int i=0; i<ntheta; i++)
{
    Dd[i] = Dd[i] / (ddsum*dang);
    //printf("Dd[%d]=%f\n",i,Dd[i]);
}


//float nang=ntheta;

float * S_array;
float * Sf;
float * St;
S_array= new float[nfreq*ntheta];
Sf= new float[nfreq];
//St= new float[ntheta];

float stsum;
float sfmax=0;


for (int i=0; i<ntheta; i++)                             //! Fill S_array
{
    stsum=0;
    for (int ii=0; ii<nfreq; ii++)
    {
        S_array[ii+i*nfreq]=y[ii]*Dd[i];// m2/Hz/rad ?
        //printf("S_array[%d,%d]=%f\n",ii+1,i+1,S_array[ii+i*nfreq]);
        if (i==0)
        {
            Sf[ii]=S_array[ii+i*nfreq]*dang;
        }
        else
        {
            Sf[ii]=Sf[ii]+S_array[ii+i*nfreq]*dang;
        }
        stsum=stsum+S_array[ii+i*nfreq];
        
    }
    Stt[i]=stsum*dfj;
    //printf("St[%d]=%f\n",i,Stt[i]);
    ///St?
    ///Sf?
}




//                                         ! integrate S_array over angles

//call tpDcalc(Sf,f,Trep)

//allocate(temp(size(Sf)))
float sumnom=0.0f;
float sumden=0.0f;

for (int i=0; i<nfreq; i++)
{
    float temp=0.0f;
    if (Sf[i]>0.8f*sfmax)
    {
        temp=1.0f;
    }
    sumnom=sumnom+temp*Sf[i]*f[i];
    sumden=sumden+temp*Sf[i];
}
//float Trep;
	
	//printf("sumnom=%f\n",sumnom);
	//printf("sumden=%f\n",sumden);

TTrep=sumnom/sumden;

TTrep=1.0f/TTrep;
	//printf("TTrep=%f\n",TTrep);

//return Trep;


} 
