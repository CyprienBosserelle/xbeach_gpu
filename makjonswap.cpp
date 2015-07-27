//////////////////////////////////////////////////////////////////////////////////
//XBeach_GPU                                                                    //
//Copyright (C) 2013 Bosserelle                                                 //
//                                                                              //
//This program is free software: you can redistribute it and/or modify          //
//it under the terms of the GNU General Public License as published by          //
//the Free Software Foundation.                                                 //
//                                                                              //
//This program is distributed in the hope that it will be useful,               //
//but WITHOUT ANY WARRANTY; without even the implied warranty of                //    
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//GNU General Public License for more details.                                  //
//                                                                              //
//You should have received a copy of the GNU General Public License             //
//along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <algorithm>

#define pi 3.14159265
using DECNUM = float;

void makjonswap(DECNUM hm0gew,DECNUM fp,DECNUM mainang,DECNUM rt,DECNUM scoeff,DECNUM gam,DECNUM *theta,int ntheta,DECNUM& TTrep,DECNUM * &Stt)
{
    

//DECNUM hm0gew=2.0f;//Significant wave height (m)
//DECNUM fp=1.0f/12.0f; //Wave peak frequency (Hz)
//DECNUM mainang=0; //wave mean direction (angle of incidence º)
//DECNUM rt=1000; //Boundary duration
//DECNUM scoeff=100;// spread coef n.u.
//DECNUM gam=3.3f;//: peak enhancement factor, optional parameter (DEFAULT 3.3)

//printf("hs=%f\nfp=%f\nang=%f\nscoeff=%f\ngam=%f\n",hm0gew,fp,mainang,scoeff,gam);	
	
	
//printf("fp=%f\n",fp);
//// 

DECNUM fnyq = 3.0f*fp;
DECNUM dfj= fp/20.0f;

// 
int nfreq=ceil((fnyq-dfj)/dfj);
DECNUM * f;
DECNUM * x;
DECNUM * y;

f=(DECNUM *)malloc(nfreq*sizeof(DECNUM));
//f= new DECNUM[nfreq];
//x= new DECNUM[nfreq];
//y= new DECNUM[nfreq];
x=(DECNUM *)malloc(nfreq*sizeof(DECNUM));
y=(DECNUM *)malloc(nfreq*sizeof(DECNUM));
//printf("Hello world!");
DECNUM xa,ymax,ysum;
DECNUM sigma,fac1,fac2,fac3;

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
        
//DECNUM t1=-(pi)/2;
//int ntheta=101;
//DECNUM *theta;

DECNUM * Dd;
Dd= (DECNUM *)malloc(ntheta*sizeof(DECNUM));
//Dd= new DECNUM[ntheta];
//theta= (DECNUM *)malloc(ntheta*sizeof(DECNUM));//new DECNUM[ntheta];
DECNUM ddsum=0.0f;
//theta=(0:100)*((pi)/100)+t1;

for(int i=0; i<ntheta; i++)
{
    //theta[i]=i*((pi)/100.0f)+t1;
    Dd[i] = pow(cos((theta[i]-mainang)/2.0f),2.0f*scoeff);
    ddsum=ddsum+Dd[i];
    //printf("theta[%d]=%f\n",i,theta[i]);
}

DECNUM dang=theta[1]-theta[0];

//mainang=(1.5d0*p1-alfa)-mainang*atan(1.d0)/45.0d0;

for(int i=0; i<ntheta; i++)
{
    Dd[i] = Dd[i] / (ddsum*dang);
    //printf("Dd[%d]=%f\n",i,Dd[i]);
}


//DECNUM nang=ntheta;

DECNUM * S_array;
DECNUM * Sf;
DECNUM * St;
S_array= new DECNUM[nfreq*ntheta];
Sf= new DECNUM[nfreq];
//St= new DECNUM[ntheta];

DECNUM stsum;
DECNUM sfmax=0;


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
DECNUM sumnom=0.0f;
DECNUM sumden=0.0f;

for (int i=0; i<nfreq; i++)
{
    DECNUM temp=0.0f;
    if (Sf[i]>0.8f*sfmax)
    {
        temp=1.0f;
    }
    sumnom=sumnom+temp*Sf[i]*f[i];
    sumden=sumden+temp*Sf[i];
}
//DECNUM Trep;
	
	//printf("sumnom=%f\n",sumnom);
	//printf("sumden=%f\n",sumden);

TTrep=sumnom/sumden;

TTrep=1.0f/TTrep;
	//printf("TTrep=%f\n",TTrep);

//return Trep;


} 
