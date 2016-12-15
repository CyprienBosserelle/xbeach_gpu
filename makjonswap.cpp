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

#include "XBeachGPU.h"

#define pi 3.14159265
using DECNUM = float;


void GenCstWave(XBGPUParam Param, std::vector<Wavebndparam> wavebnd, float * theta, double * &Stfile, double * &qfile, double * &Tpfile)
{
	int ny, ntheta,nwavbnd;
	double dtheta;
	ny = Param.ny;
	ntheta = Param.ntheta;
	dtheta = Param.dtheta;

	nwavbnd = wavebnd.size();

	for (int n = 0; n < nwavbnd; n++)
	{
		Tpfile[n] = wavebnd[n].Tp;
		double eetot = wavebnd[n].Hs*wavebnd[n].Hs*Param.rho*Param.g / (16.0*Param.dtheta);

		for (int j = 0; j < Param.ny; j++)
		{
			qfile[j + 0 * ny + n*ny * 4] = 0.0;
			qfile[j + 1 * ny + n*ny * 4] = 0.0;
			qfile[j + 2 * ny + n*ny * 4] = 0.0;
			qfile[j + 3 * ny + n*ny * 4] = 0.0;
		}

		double sumcos = 0.0;

		double * scaledir;

		scaledir = (double *)malloc(ntheta*sizeof(double));

		for (int t = 0; t < ntheta; t++)
		{
			scaledir[t] = pow(cos((theta[t] - wavebnd[n].Dp) / 2.0), 2.0*wavebnd[n].s);
			sumcos = sumcos + scaledir[t];
		}

		for (int t = 0; t < ntheta; t++)
		{
			double Stdir = scaledir[t] / sumcos*eetot;
			for (int j = 0; j < ny; j++)
			{
				Stfile[j + t*ny + n*ny*ntheta] = Stdir;
			}

		}
	}
}
	

void makjonswap(XBGPUParam Param, std::vector<Wavebndparam> wavebnd, int step, int &nfreq, int &ntheta, double * &HRfreq, double * &HRdir, double * &HRSpec)
{
    
	double * x, *y;

	double Hs = wavebnd[step].Hs;
	double Tp = wavebnd[step].Tp;
	double Dp = wavebnd[step].Dp; // converted to the main angle already
	double mainang = Dp;
	double fp = 1 / Tp;
	double gam = wavebnd[step].gamma;
	double scoeff = wavebnd[step].s;

	//printf("fp=%f\n",fp);
	//// 

	double fnyq = 3.0f*fp;
	double dfj= fp/20.0f;

	// 
	int nfreq=ceil((fnyq-dfj)/dfj);


	HRfreq=(double *)malloc(nfreq*sizeof(double));
	//f= new DECNUM[nfreq];
	//x= new DECNUM[nfreq];
	//y= new DECNUM[nfreq];
	x=(double *)malloc(nfreq*sizeof(double));
	y=(double *)malloc(nfreq*sizeof(double));
	//printf("Hello world!");
	double xa,ymax,ysum;
	double sigma,fac1,fac2,fac3;

	ysum=0.0f;
	ymax=0.0f;

	for (int i=0; i<nfreq ;i++)
	{
		HRfreq[i] = (i + 1)*dfj;
    
		// x: nondimensional frequency, divided by the peak frequency
		x[i] = HRfreq[i] / fp;
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
    
		y[i]=y[i]*pow((Hs/(4*sqrt(ysum*dfj))),2);
		//printf("y[%d]=%f\n",i,y[i]);
	}


	////
        
	//DECNUM t1=-(pi)/2;
	//int ntheta=101;
	//DECNUM *theta;
	double * HRtheta;
	int ntheta = 90;
	double dtheta = 2*pi / ntheta;
	double *theta;


	double * Dd;
	
	theta = (double *)malloc(ntheta*sizeof(double));
	Dd = (double *)malloc(ntheta*sizeof(double));
	HRdir = (double *)malloc(ntheta*sizeof(double));

	
	//Dd= new DECNUM[ntheta];
	//theta= (DECNUM *)malloc(ntheta*sizeof(DECNUM));//new DECNUM[ntheta];
	double ddsum = 0.0f;
	//theta=(0:100)*((pi)/100)+t1;

	for(int i=0; i<ntheta; i++)
	{
		theta[i]=i*dtheta-pi; // cover the full circle
		HRdir[i] = theta[i];
		Dd[i] = pow(cos((theta[i]-mainang)/2.0f),2.0f*scoeff);
		ddsum=ddsum+Dd[i];
		//printf("theta[%d]=%f\n",i,theta[i]);
	}

	double dang=theta[1]-theta[0];

	//mainang=(1.5d0*p1-alfa)-mainang*atan(1.d0)/45.0d0;
	
	for(int i=0; i<ntheta; i++)
	{
		Dd[i] = Dd[i] / (ddsum*dang);
		//printf("Dd[%d]=%f\n",i,Dd[i]);
	}


	//DECNUM nang=ntheta;

	HRSpec = (double *)malloc(nfreq*ntheta*sizeof(double));
	

	for (int i=0; i<ntheta; i++)                             //! Fill S_array
	{
		
		for (int ii=0; ii<nfreq; ii++)
		{
			HRSpec[ii + i*nfreq] = y[ii] * Dd[i];// m2/Hz/rad ?
			//printf("S_array[%d,%d]=%f\n",ii+1,i+1,S_array[ii+i*nfreq]);
			
        
		}
		
	}





} 

void GenWGnLBW(XBGPUParam Param, int nf, int ndir,double * HRfreq,double * HRdir, double * HRSpec, float Trep, double * qfile, double * Stfile)
{
	int ny = Param.ny;
	int K; //Number of wave components

	int nhf;
	int nhd;

	double * Sf;

	double fmax,Sfmax; // Should be in Param

	int Kmin = 200;

	double dtheta = HRdir[1] - HRdir[0];

	Sf = (double *)malloc(nf*sizeof(double));

	for (int n = 0; n < nf; n++)
	{
		Sf[n] = 0.0;
		for (int d = 0; d < ndir; d++)
		{
			//
			Sf[n] = Sf[n] + HRSpec[d+n*nf];
		}
		Sf[n] = Sf[n] * dtheta;
	}

	//////////////////////////////////////
	// Generate wave train component
	//////////////////////////////////////
	fmax = 0.0;
	Sfmax = 0.0;
	for (int n = 0; n < nf; n++)
	{
		fmax = max(fmax,HRfreq[n]);
		Sfmax = max(Sfmax, Sf[n]);
	}
	fmax = 2.0*fmax; //???? 
	
	int ind1 = 0;
	int ind2 = nf - 1;

	int first = 0;


	for (int n = 0; n < nf; n++)
	{
		if (Sf[n]>Sfmax*Param.sprdthr && HRfreq[n]<fmax )
		{
			if (first == 0)
			{
				ind1 = n;
				first = 1;
			}
			ind2 = n;
		}
	}

	//Calculate number of wave components to be included in determination of the
	//wave boundary conditions based on the wave record length and width of the
	//wave frequency range
	K = ceil(Param.dtbc*(HRfreq[ind2] - HRfreq[ind1]) + 1);
	//also include minimum number of components
	K = max(K, Kmin);





	//////////////////////////////////////
	// Generate wave time axis
	//////////////////////////////////////

	//////////////////////////////////////
	//Generate wave train variance
	//////////////////////////////////////

	//////////////////////////////////////
	//Generate wave train properties at each offshore points
	//////////////////////////////////////

	//////////////////////////////////////
	//Generate wave train fourier
	//////////////////////////////////////

	//////////////////////////////////////
	//Distribute wave train direction
	//////////////////////////////////////

	//////////////////////////////////////
	// Generate e (STfile)
	//////////////////////////////////////

	//////////////////////////////////////
	// Generate q (qfile)
	//////////////////////////////////////

	//////////////////////////////////////
	//Save to Netcdf file
	//////////////////////////////////////


}
