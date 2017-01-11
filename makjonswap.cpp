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
	nfreq=ceil((fnyq-dfj)/dfj);


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
	ntheta = 90;
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

	double * Sf; // size of nf

	double *fgen, *phigen, *thetagen, *kgen, *wgen, *vargen;
	int *Findex , *WDindex; // size of K
	//double *CompFn;//size of ny*K // now as a 2D vector of complex
	double *Sd,*pdf, *cdf; // size of ndir

	double fmax,Sfmax; // Should be in Param
	int nspr = 0; // Should be in Param (need test for nspr ==1)
	double * binedgeleft, * binedgeright; // size of ntheta

	int Kmin = 200;

	double dtheta = HRdir[1] - HRdir[0];
	double dfreq = HRfreq[1] - HRfreq[0];
	double Hm0=0.0;
	Sf = (double *)malloc(nf*sizeof(double));
	Sd = (double *)malloc(ndir*sizeof(double));
	pdf = (double *)malloc(ndir*sizeof(double));
	cdf = (double *)malloc(ndir*sizeof(double));

	for (int n = 0; n < nf; n++)
	{
		Sf[n] = 0.0;
		for (int d = 0; d < ndir; d++)
		{
			//
			Sf[n] = Sf[n] + HRSpec[n+d*nf];
		}
		Sf[n] = Sf[n] * dtheta;
		Hm0 = Hm0 + Sf[n];
	}
	Hm0 = 4.0*sqrt(Hm0*dfreq);

	for (int d = 0; d < ndir; d++)
	{
		Sd[d] = 0.0;
		for (int n = 0; n < nf; n++)
		{
			Sd[d] = Sd[d] + HRSpec[n + d*nf];
		}
		Sd[d] = Sd[d] * dfreq;
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
	K = ceil(Param.rtlength*(HRfreq[ind2] - HRfreq[ind1]) + 1);
	//also include minimum number of components
	K = (int)max(K*1.0, Kmin*1.0);// workaround because template for int not compiling for some reason

	fgen = (double *)malloc(K*sizeof(double));
	phigen = (double *)malloc(K*sizeof(double));
	thetagen = (double *)malloc(K*sizeof(double));
	kgen = (double *)malloc(K*sizeof(double));
	wgen = (double *)malloc(K*sizeof(double));
	vargen = (double *)malloc(K*sizeof(double));
	Findex = (int *)malloc(K*sizeof(int));
	WDindex = (int *)malloc(K*sizeof(int));
	//CompFn = (double *)malloc(K*Param.ny*sizeof(double));

	double dfgen = (HRfreq[ind2] - HRfreq[ind1]) / K;
	for (int i = 0; i < K; i++)
	{
		//
		fgen[i] = HRfreq[ind1] + i*dfgen;
	}

	unsigned seed;
	if (Param.random == 1)
	{
		seed = std::chrono::system_clock::now().time_since_epoch().count();
	}
	else
	{
		seed = 0;
	}


	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	//Determine a random phase for each wave train component between 0 and 2pi

	for (int i = 0; i < K; i++)
	{
		//
		phigen[i] = distribution(generator)*2.0*pi;
	}
	//double number = distribution(generator);

	//Determine random directions for each wave train component, based on the CDF of
	//the directional spectrum.For each wave train we will interpolate the directional
	//distribution at that frequency, generate a CDF, and then interpolate the wave
	//train direction from a random number draw and the CDF.
	for (int i = 0; i < K; i++)
	{
		//int fprev = ceil((fgen[i] - HRfreq[ind1]) / dfreq) + ind1;

		for (int d = 0; d < ndir; d++)
		{
			//pdf[d] = interptime(HRSpec[d + (fprev + 1)*nf], HRSpec[d + (fprev)*nf], dfreq, fgen[i] - HRfreq[fprev]);
			pdf[d] = Interp2(nf, ndir,HRfreq, HRdir,HRSpec, fgen[i],HRdir[d]);
		}

		//convert to pdf by ensuring total integral == 1, assuming constant directional bin size
		double sumpdf = 0;
		for (int d = 0; d < ndir; d++)
		{
			sumpdf = sumpdf + pdf[d];
		}

		for (int d = 0; d < ndir; d++)
		{
			pdf[d] = pdf[d] / sumpdf;
		}

		//convert to cdf by trapezoidal integration
		//Note: this only works if the directional
		//bins are constant in size.Assumed multiplication by one.
		cdf[0] = pdf[0];
		for (int d = 1; d < ndir; d++)
		{
			cdf[d] = cdf[d - 1] + 0.5*(pdf[d - 1] + pdf[d]);
		}
		double number = distribution(generator);

		//int dprev = ndir-1;
		//for (int d = 1; d < ndir; d++)
		//{
		//	double diff = number - cdf[d];
		//	if (diff < 0.0)
		//	{
		//		dprev = d - 1;
		//		break;
		//	}
		//}
		//interp1D(double *x, double *y, double xx)
		//thetagen[i] = interptime(HRdir[dprev + 1], HRdir[dprev], dtheta, dtheta*((number - cdf[dprev]) / (cdf[dprev + 1] - cdf[dprev])));
		thetagen[i] = interp1D(ndir, cdf, HRdir, number);
		printf("thetagen[i]=%f\n", thetagen[i]);
	}

	//determine wave number for each wave train component
	// RL Soulsby(2006) "Simplified calculation of wave orbital velocities"
	// HR Wallingford Report TR 155, February 2006
	// Eqns. 12a - 14

	// csherwood@usgs.gov
	// Sept 10, 2006
	for (int i = 0; i < K; i++)
	{
		double L, L0,w,x,y,h,t;
		h = Param.offdepth;
		L0 = Param.g*(1 / fgen[i])*(1 / fgen[i]) / 2 / pi;
		w = 2 * pi / ((1 / fgen[i]));//2pi/T
		x = w*w * h / Param.g;
		y = sqrt(x) * (x<1) + x* (x >= 1);
		t = tanh(y);
		y = y - ((y*t - x) / (t + y*(1 - t*t)));
		t = tanh(y);
		y = y - ((y*t - x) / (t + y*(1 - t*t)));
		t = tanh(y);
		y = y - ((y*t - x) / (t + y*(1 - t*t)));
		kgen[i] = y / h;
		wgen[i] = w;
	}

	


	//////////////////////////////////////
	// Generate wave time axis
	//////////////////////////////////////

	//First assume that internal and bc - writing time step is the same
	double dtin = Param.dtbc;
	int tslenbc = (int)(Param.rtlength / Param.dtbc) + 1;
	
	//Check whether the internal frequency is high enough to describe the highest frequency
	//wave train returned from frange(which can be used in the boundary conditions)
	if (dtin > 0.1 / fgen[K - 1])
	{
		dtin = 0.1 / fgen[K - 1];
	}
	//! The length of the internal time axis should be even (for Fourier transform) and
	//depends on the internal time step needed and the internal duration(~1 / dfgen) :
	int tslen = (int)(ceil(1 / dfgen / dtin) + 1);
	if (ceil(tslen/2)-tslen/2>0)
	{
		tslen = tslen + 1;
	}
	
	//Now we can make the internal time axis
	double rtin = tslen * dtin;
	double * tin, *taperf, *taperw;
	tin = (double *)malloc(tslen*sizeof(double));

	for (int n = 0; n < tslen; n++)
	{
		tin[n] = n*dtin;
	}

	//Make a taper function to slowly increase and decrease the boundary condition forcing
	//at the start and the end of the boundary condition file(including any time beyond
	//the external rtbc
	taperf = (double *)malloc(tslen*sizeof(double));
	taperw = (double *)malloc(tslen*sizeof(double));

	for (int n = 0; n < tslen; n++)
	{
		taperf[n] = 1;
		taperw[n] = 1;
	}

	double Tbc = 1 / fgen[0]; //Should be Trep or 1/fpeak...
	int ntaper = (int)((5.0*Tbc) / dtin);

	for (int n = 0; n < (int)min(1.0*ntaper, 1.0*tslen); n++)
	{
		taperf[n] = tanh(5.0*n / ntaper); //
		taperw[n] = tanh(5.0*n / ntaper); //
	}

	//We do not want to taperw the end anymore.Instead we pass the wave height at the end of rtbc to
	//the next wave generation iteration.
	//end taper by finding where tin = rtbc, taper before that and set everything to zero after
	//that.

	for (int n = tslen; n > 0; n--)
	{
		if (tin[n - 1] > (Param.rtlength - ntaper*dtin))
		{
			taperf[n - 1] = tanh(5.0*(Param.rtlength-tin[n-1])/dtin/ntaper);
		}
		
		
		if (tin[n - 1] > Param.rtlength)
		{
			taperf[n - 1] = 0;
		}

	}

	//////////////////////////////////////
	//Generate wave train variance
	//////////////////////////////////////
	//! Determine variance at each spectrum location
	double sumvargen=0.0;
	for (int i = 0; i < K; i++)
	{
		//interptime(double next, double prev, double timenext, double time)
		//vargen[i] = interp1D(HRfreq, Sf, fgen[i]);
		vargen[i] = Interp2( nf, ndir, HRfreq, HRdir, HRSpec, fgen[i], thetagen[i]);
		sumvargen = sumvargen + vargen[i];
	}

	// scale vargen so that 4*sqrt(sum(vargen)*dfgen)==4*sqrt(sum(Sf)*df)
	double scalefactor = pow(Hm0 / (4.0 * sqrt(sumvargen*dfgen)),2); // squared

	for (int i = 0; i < K; i++)
	{
		vargen[i] = vargen[i] * scalefactor;
	}

	//Not sure why this is done here (the prev values should alway be the min anyways)
	double dummy;
	for (int i = 0; i < K; i++)
	{
		dummy = interp1D(nf, HRfreq, Sf, fgen[i]);
		vargen[i] = min(vargen[i] ,dummy);
	}

	//////////////////////////////////////
	//Generate wave train properties at each offshore points
	//////////////////////////////////////
	//Skip A=sqrt(2*vargen[i]*dfgen) and Hmo=4.0 * sqrt(sumvargen*dfgen)

	//////////////////////////////////////
	//Generate wave train fourier
	//////////////////////////////////////
	//  ! Determine indices of wave train components in frequency axis and
	//  !Fourier transform result
	
	int tempi = floor(fgen[0] / dfgen);
	for (int i = 0; i < K; i++)
	{
		Findex[i] = tempi + i + 1; // why not simply i ??
	}
	
	//	! Determine first half of complex Fourier coefficients of wave train
	//	!components using random phase and amplitudes from sampled spectrum
	//	!until Nyquist frequency.The amplitudes are represented in a
	//	!two - sided spectrum, which results in the factor 1 / 2.
	//	!Unroll Fourier components along offshore boundary, assuming all wave trains
	//	!start at x(1, 1), y(1, 1).


	std::complex<double> par_compi (0.0,1.0);
	std::complex<double> tempcmplx;
	//= 0.0 + 1.0*I;
	TwoDee<std::complex<double>> CompFn(ny, tslen);


	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < Param.ny; j++)
		{
			CompFn(j,Findex[i]) = sqrt(2 * vargen[i] * dfgen) / 2 * exp(par_compi*phigen[i])*    //Bas: wp%Findex used in time dimension because t = j*dt in frequency space
				exp(-par_compi*kgen[i] * (sin(thetagen[i])*(j*Param.dx))); //dsin
				//+ cos(thetagen[i])*(xb[j] - x0))); //dcos

			//!Determine Fourier coefficients beyond Nyquist frequency(equal to
			//!coefficients at negative frequency axis) of relevant wave components for
			//!first y - coordinate by mirroring
			for (int n = 1; n < (tslen / 2); n++)
			{
				tempcmplx = std::conj(CompFn(j, n));
				CompFn(j, tslen - (tslen / 2 + n)) = tempcmplx; //Not sure this is right
			}


		}
	}
	//create2dnc(nfHR, ndHR, HRfreq[1] - HRfreq[0], HRdir[1] - HRdir[0], 0.0, HRfreq, HRdir, HRSpec);
	//////////////////////////////////////
	//Distribute wave train direction
	//////////////////////////////////////

	//!Calculate the bin edges of all the computational wave bins in the
	//!XBeach model(not the input spectrum)
	binedgeleft = (double *)malloc(Param.ntheta*sizeof(double));
	binedgeright = (double *)malloc(Param.ntheta*sizeof(double));
	
	for (int i = 0; i < Param.ntheta; i++)
	{
		binedgeleft[i] = fmod(i*Param.dtheta + Param.thetamin,2*pi);
		binedgeright[i] = fmod((i+1)*Param.dtheta + Param.thetamin,2*pi);
	}


	//All generated wave components are in the rang 0 <= theta<2pi.
	//	!We link wave components to a wave direction bin if the direction falls
	//	!within the bin boundaries.Note the >= and <= , ranther than >= and <.This
	//	!is not necessarily a problem, but solves having to make an exception for the
	//	!highest wave direction bin, in which >= and <= should be applicable.
	//	!In the case of a central bin and a wave direction exactly(!) on the bin
	//	!interface, the wave will be linked to the first wave bin in the ordering,
	//	!rather than the higher of the two bins.
	//	!
	//	!Initially set WDindex to zero.This marks a wave direction outside the
	//	!computational wave bins.In case it does not fit in a directional wave
	//	!bin, it remains zero at the end of the loops.
	//	!Note: this does not ensure all energy is included in the wave bins,
	//	!as wave energy may still fall outside the computational domain.

	
	for (int i = 0; i < K; i++)
	{
		WDindex[i] = 0;
		for (int itheta = 0; itheta < Param.ntheta; i++)
		{
			// special case if this bin spans 0 degrees
			if (binedgeleft[itheta]>binedgeright[itheta])
			{
				if ((thetagen[i] >= binedgeleft[itheta] && thetagen[i] <= (2 * pi)) || (thetagen[i] >= 0.0 && thetagen[i] <= binedgeright[itheta]))
				{
					WDindex[i] = itheta;
					// We now have the correct wave bin, move to next wave component K
					break;
				}
			}
			else
			{
				if (thetagen[i] >= binedgeleft[itheta] && thetagen[i] <= binedgeright[itheta])
				{
					WDindex[i] = itheta;
					// We now have the correct wave bin, move to next wave component K
					break;
				
				}
			}

		}
	}

	// If the user has set nspr == 1 then the randomly drawn wave directions
	// should be set to the centres of the wave directional bins.
	// Also move all wave energy falling outside the computational bins, into
	// the computational domain(in the outer wave direction bins)

	if (nspr == 1)
	{
		for (int i = 0; i < K; i++)
		{
			if (WDindex[i]>0)
			{
				thetagen[i] = binedgeleft[WDindex[i]] + Param.dtheta*0.5;
			}

		}
	}

	// Check the amount of energy lost to wave trains falling outside the computational domain

	double lostvar = 0.0;
	double keptvar = 0.0;
	for (int i = 0; i < K; i++)
	{
		if (WDindex[i] == 0)
		{
			lostvar = lostvar + sqrt(2 * vargen[i] * dfgen);
		}
		else
		{ 
			keptvar = keptvar + sqrt(2 * vargen[i] * dfgen);
		}
			
	}

	double perclost = 100 * (lostvar / (lostvar + keptvar));
	if (perclost > 10.0)
	{
		write_text_to_log_file("Large amounts of energy (" + std::to_string(perclost) + "%) fall outside computational domain at the offshore boundary");
	}
	else
	{
		write_text_to_log_file("Wave energy outside computational domain at offshore boundary: " + std::to_string(perclost) + "%");
	}



	//////////////////////////////////////
	// Generate e (STfile)
	//-------- - Calculate energy envelope time series from-------- -
	//--------Fourier components, and write to output file--------
	//////////////////////////////////////

	// Calculate wave energy for each y - coordinate along seaside boundary for
	/// current computational directional bin

	for (int itheta = 0; itheta < Param.ntheta; i++)
	{

	}

	//////////////////////////////////////
	// Generate q (qfile)
	//////////////////////////////////////

	//////////////////////////////////////
	//Save to Netcdf file
	//////////////////////////////////////


}
