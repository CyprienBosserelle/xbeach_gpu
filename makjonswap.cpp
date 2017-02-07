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
		free(scaledir);
	}
}
	

void makjonswap(XBGPUParam Param, std::vector<Wavebndparam> wavebnd, int step, int &nfreq, int &ntheta, double * &HRfreq, double * &HRdir, double * &HRSpec)
{
    
	double * x, *y;

	double Hs = wavebnd[step].Hs;
	double Tp = max(wavebnd[step].Tp,1.5); // for very small Tp the wave group generator will request a Giant amount of memory so we need to cap it here
	double Dp = wavebnd[step].Dp; // converted to the main angle already
	double mainang = Dp;
	double fp = 1 / Tp;
	double gam = wavebnd[step].gamma;
	double scoeff = max(round(min(wavebnd[step].s,1000.0)),10.0);

	
	printf("Generating JONSWAP spectrum: Hs=%f, Tp=%f, Dp=%f, gam=%f, s=%f\n",Hs,Tp,Dp,gam,scoeff);
	write_text_to_log_file("Generating JONSWAP spectrum: Hs=" + std::to_string(Hs) + " Tp=" + std::to_string(Tp) + " Dp=" + std::to_string(Dp) + " gam=" + std::to_string(gam) + " s=" + std::to_string(scoeff) + ";");
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
		Dd[i] = pow(cos((theta[i]-mainang)/2.0f),2.0*scoeff);
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
			HRSpec[ii + i*nfreq] = y[ii] * Dd[i];// m2/Hz/rad
			//printf("S_array[%d,%d]=%f\n",ii+1,i+1,S_array[ii+i*nfreq]);


			if (HRSpec[ii + i*nfreq] != HRSpec[ii + i*nfreq])
			{
				printf("Error in generating JONSWAP Spectrum: #NAN or #IND detected");
				write_text_to_log_file("Error in generating JONSWAP Spectrum: #NAN or #IND detected");
				exit(EXIT_FAILURE);
			}
        
		}
		
	}


	free(x);
	free(y);
	free(theta);
	free(Dd);




} 

void GenWGnLBW(XBGPUParam Param, int nf, int ndir,double * HRfreq,double * HRdir, double * HRSpec, float &Trep, double * &qfile, double * &Stfile)
{
	// Generating Boundary condition: Energy from wave group and Long bound waves

	printf("Generating Boundary condition: Energy from wave group and Long bound waves.\n");
	write_text_to_log_file("Generating Boundary condition: Energy from wave group and Long bound waves");

	int ny = Param.ny;
	int K; //Number of wave components

	//int nhf;
	//int nhd;

	double * Sf; // size of nf

	double trepfac = 0.01;
	double temptrep;

	double *fgen, *phigen, *thetagen, *kgen, *wgen, *vargen,*vargenq;
	int *Findex , *WDindex; // size of K
	//double *CompFn;//size of ny*K // now as a 2D vector of complex
	double *Sd,*pdf, *cdf; // size of ndir

	double fmax,Sfmax; // Should be in Param
	//int nspr = 0; // Should be in Param (need test for nspr ==1)
	double * binedgeleft, * binedgeright; // size of ntheta
	double * zeta, *Ampzeta; //water elevation ny*ntheta*tslen
	double *eta, *Amp; //water elevation integrated over directions ny*tslen
	double *stdzeta; //size of ntheta
	double *E_tdir; // tslen
	double *qx, *qy, *qtot, *qtempx, *qtempy;

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
	//printf("Hm0=%f\n", Hm0); 

	//Need a sanity check here!

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
	vargenq = (double *)malloc(K*sizeof(double));
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


	std::default_random_engine generator (seed);
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
		
		//Cannot use interp1DMono here because cdf is not monotonic
		thetagen[i] = interp1D(ndir, cdf, HRdir, number);

		//printf("thetagen[i]=%f\n", thetagen[i]);
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
	// Calculate Trep
	//////////////////////////////////////
	
	temptrep = 0.0;
	double tempf = 0.0;
	for (int n = 0; n < nf; n++)
	{

		if (Sf[n] >= Sfmax*trepfac)
		{
			temptrep += (Sf[n] / max(HRfreq[n], 0.001));
			tempf += Sf[n];

		}
	}
	Trep = (float) (temptrep / tempf);



	//////////////////////////////////////
	// Generate wave time axis
	//////////////////////////////////////

	//First assume that internal and bc - writing time step is the same
	double dtin = Param.dtbc;
	int tslenbc = ceil(Param.rtlength / Param.dtbc) + 1;// (int)(Param.rtlength / Param.dtbc) + 1;
	
	//Check whether the internal frequency is high enough to describe the highest frequency
	//wave train returned from frange(which can be used in the boundary conditions)
	if (dtin > 0.5 / fgen[K - 1])
	{
		dtin = 0.5 / fgen[K - 1];
	}
	//! The length of the internal time axis should be even (for Fourier transform) and
	//depends on the internal time step needed and the internal duration(~1 / dfgen) :
	int tslen = (int)(ceil(1 / dfgen / dtin) + 1);
	if ((ceil(tslen/2)-tslen/2)>0)
	{
		tslen = tslen + 1;
	}
	
	//int tsfft = pow(2,ceil(ceil(log(tslen) / log(2)))); /// for fft

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
	double Hm0post = 4.0 * sqrt(sumvargen*dfgen);
	double scalefactor = pow(Hm0 / Hm0post, 2); // squared

	for (int i = 0; i < K; i++)
	{
		vargen[i] = vargen[i] * scalefactor;
	}

	//Not sure why this is done here (the prev values should alway be the min anyways)
	double dummy;
	for (int i = 0; i < K; i++)
	{
		dummy = interp1DMono(nf, HRfreq, Sf, fgen[i]);
		vargenq[i] = min(vargen[i] ,dummy);
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
		Findex[i] = tempi + i; 
	}
	
	//	! Determine first half of complex Fourier coefficients of wave train
	//	!components using random phase and amplitudes from sampled spectrum
	//	!until Nyquist frequency.The amplitudes are represented in a
	//	!two - sided spectrum, which results in the factor 1 / 2.
	//	!Unroll Fourier components along offshore boundary, assuming all wave trains
	//	!start at x(1, 1), y(1, 1).


	std::complex<double> par_compi (0.0,1.0);
	//std::complex<double> tempcmplx;
	CArray tempcmplx(tslen / 2-1);
	//= 0.0 + 1.0*I;
	TwoDee<std::complex<double>> CompFn(ny, tslen);

	zeta = (double *)malloc(ny*Param.ntheta*tslen*sizeof(double));
	//Ampzeta = (double *)malloc(ny*Param.ntheta*tslen*sizeof(double));
	eta = (double *)malloc(ny*tslen*sizeof(double));
	Amp = (double *)malloc(ny*tslen*sizeof(double));
	stdzeta = (double *)malloc(Param.ntheta*sizeof(double));
	E_tdir = (double *)malloc(tslen*sizeof(double));
	qx = (double *)malloc(ny*tslen*sizeof(double));
	qy = (double *)malloc(ny*tslen*sizeof(double));
	qtot = (double *)malloc(ny*tslen*sizeof(double));
	qtempx = (double *)malloc(tslen*sizeof(double));
	qtempy = (double *)malloc(tslen*sizeof(double));

	//initialise the variables
	for (int n = 0; n < ny*Param.ntheta*tslen;n++)
	{
		zeta[n] = 0.0;
		//Ampzeta[n] = 0.0;
	}
	for (int n = 0; n < ny*tslen; n++)
	{
		eta[n] = 0.0;
		Amp[n] = 0.0;
	}

	for (int n = 0; n < Param.ntheta; n++)
	{
		stdzeta[n] = 0.0;
	}

	for (int n = 0; n <tslen; n++)
	{
		E_tdir[n] = 0.0;
	}





	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < Param.ny; j++)
		{
			CompFn(j,Findex[i]) = sqrt(2.0 * vargen[i] * dfgen) / 2 * exp(par_compi*phigen[i])*    //Bas: wp%Findex used in time dimension because t = j*dt in frequency space
				exp(-par_compi*kgen[i] * (sin(thetagen[i])*(j*Param.dx))); //dsin
				//+ cos(thetagen[i])*(xb[j] - x0))); //dcos

			//!Determine Fourier coefficients beyond Nyquist frequency(equal to
			//!coefficients at negative frequency axis) of relevant wave components for
			//!first y - coordinate by mirroring
			
			
			for (int n = 1; n < (tslen / 2); n++)
			{
				tempcmplx[n-1] = std::conj(CompFn(j, n));
			}

			flipiv(tempcmplx);

			for (int n = (tslen/2+1); n < tslen; n++)
			{
				CompFn(j, n) = tempcmplx[n - (tslen / 2 + 1)]; //Not sure this is right
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
		binedgeleft[i] = fmod((i*Param.dtheta + Param.thetamin),2*pi);
		binedgeright[i] = fmod(((i+1)*Param.dtheta + Param.thetamin),2*pi);
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
		WDindex[i] = -1;
		for (int itheta = 0; itheta < Param.ntheta; itheta++)
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

	if (Param.nspr == 1)
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
		if (WDindex[i] < 0)
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
	
	std::vector<int> wcompindx;
	CArray tempcplxarr(tslen*0.5 - 1);

	//std::valarray<double> zeta(tslen);
	CArray Gn(tslen);
	fftw_complex *out, *in;
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * tslen);
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * tslen);
	fftw_plan p;
	p = fftw_plan_dft_1d(tslen, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	
	for (int itheta = 0; itheta < Param.ntheta; itheta++)
	{
		//Check whether any wave components are in the current computational
		//directional bin
		for (int i = 0; i < K; i++)
		{
			if (WDindex[i] == itheta)
			{
				wcompindx.push_back(i);
			}
		}

		if (wcompindx.empty())//(waveinbin == 0)
		{
			// no wave component in bin so nothing to do!
			continue;
		}
		// else There are some wave component in the bin
		


		for (int j = 0; j < Param.ny; j++)
		{
			
			Gn = 0;// Reset the whole array
			//GnForFFT = 0;
			// Determine Fourier coefficients of all wave components for current
			// y - coordinate in the current computational directional bin
			for (int n = 0; n < wcompindx.size(); n++)
			{
				//printf("wcompindx[%d]=%d; Findex[%d]=%d\n", n, wcompindx[n], n, Findex[n]);
				Gn[Findex[wcompindx[n]]] = (CompFn(j, Findex[wcompindx[n]]));
			}

			tempcplxarr = Gn[std::slice(1, tslen*0.5 - 1, 1)];
			for (int jcplx = 0; jcplx < tempcplxarr.size(); jcplx++)
			{
				tempcplxarr[jcplx] = std::conj(tempcplxarr[jcplx]);
			}
			flipiv(tempcplxarr);
			Gn[std::slice(tslen / 2 + 1, tempcplxarr.size(), 1)] = tempcplxarr;
						
			for (int n = 0; n < tslen; n++)
			{
				in[n][0] = std::real(Gn[n]);
				in[n][1] = std::imag(Gn[n]);
			}

			
			// Inverse Discrete Fourier transformation to transform back to time
			// domain from frequency domain
			fftw_execute(p);
			
			//store the results in zeta
			for (int n = 0; n < tslen; n++)
			{
				zeta[j + itheta*ny + n*ny*Param.ntheta] = out[n][0]* taperw[n];
				
			}
			

			
		}

		wcompindx.clear();
	}

	//Temporarily output results for debugging
	
	/*double * yyfx, *thetafx;

	yyfx=(double *)malloc(ny*sizeof(double));
	thetafx=(double *)malloc(Param.ntheta*sizeof(double));

	for (int j = 0; j < Param.ny; j++)
	{
		yyfx[j] = j*Param.dx;
	}

	for (int itheta = 0; itheta < Param.ntheta; itheta++)
	{
		thetafx[itheta] = itheta*(Param.dtheta) + Param.thetamin + 0.5f*Param.dtheta;
	}
	*/
	
	
	// Calculate energy envelope amplitude
	
	//Integrate instantaneous water level excitation of wave
	//components over directions
	double temp = 0.0;
	CArray tmpcplx(tslen);
	for (int j = 0; j < Param.ny; j++)
	{
		tmpcplx = 0.0; //need to reset 
		for (int n = 0; n < tslen; n++)
		{
			temp = 0.0;
			for (int itheta = 0; itheta < Param.ntheta; itheta++)
			{
				temp = temp + zeta[j + itheta*ny + n*ny*Param.ntheta];
			}
			eta[j + n*ny] = temp;
			tmpcplx[n] = temp;
		}

		hilbert(tmpcplx);
		
		for (int n = 0; n < tslen; n++)
		{
			Amp[j + n*ny] = abs(tmpcplx[n]);
		}

		double stdeta = 0.0;
		for (int n = 0; n < tslen; n++)
		{
			stdeta = stdeta + eta[j + n*ny] * eta[j + n*ny];
		}

		stdeta = sqrt(stdeta / (tslen - 1));
		for (int itheta = 0; itheta < Param.ntheta; itheta++)
		{
			temp = 0.0;

			for (int n = 0; n < tslen; n++)
			{
				temp = temp + zeta[j + itheta*ny + n*ny*Param.ntheta] * zeta[j + itheta*ny + n*ny*Param.ntheta];

			}
			stdzeta[itheta] = sqrt(temp / (tslen - 1));

			for (int n = 0; n < tslen; n++)
			{
				zeta[j + itheta*ny + n*ny*Param.ntheta] = Amp[j + n*ny] * stdzeta[itheta] / stdeta;
			}

		}

		//Calculate energy and interpolate to the output time step
		// (Maybe dtin==dtbc or maybe not)
		for (int itheta = 0; itheta < Param.ntheta; itheta++)
		{
			for (int n = 0; n < tslen; n++)
			{
				E_tdir[n] = zeta[j + itheta*ny + n*ny*Param.ntheta] * zeta[j + itheta*ny + n*ny*Param.ntheta] * 0.5*Param.rho*Param.g / (Param.dtheta);
				if (E_tdir[n] != E_tdir[n])
				{
					printf("Error in generating Wave component eebc: #NAN or #IND detected");
					write_text_to_log_file("Error in generating Wave component eebc: #NAN or #IND detected");
					exit(EXIT_FAILURE);
				}
			}
			//interpolate to boundary timeseries
			for (int m = 0; m < tslenbc; m++)
			{
				Stfile[j + itheta*ny + m*ny*Param.ntheta] = interp1DMono(tslen, tin, E_tdir, m*Param.dtbc);
			}

		}


		


	}

	//Temporarily output results for debugging
	/*
	double * yyfx, *thetafx, *bctimin;

	yyfx = (double *)malloc(ny*sizeof(double));
	thetafx = (double *)malloc(Param.ntheta*sizeof(double));
	bctimin = (double *)malloc(tslenbc*sizeof(double));

	for (int j = 0; j < Param.ny; j++)
	{
	yyfx[j] = j*Param.dx;
	}

	for (int itheta = 0; itheta < Param.ntheta; itheta++)
	{
	thetafx[itheta] = itheta*(Param.dtheta) + Param.thetamin + 0.5f*Param.dtheta;
	}

	for (int m = 0; m < tslenbc; m++)
	{
		bctimin[m] = m*Param.dtbc;
	}

	create3dnc(ny, Param.ntheta, tslenbc, Param.dx, Param.dtheta, Param.dtbc, 0.0, yyfx, thetafx, bctimin, Stfile);
	*/

	//////////////////////////////////////
	// Bound long waves
	//////////////////////////////////////
	double deltaf, deltatheta, k3, Eforc;
	double *KKx, *KKy, *dphi3, *cg3, *theta3, *D, *Abnd;

	//deltatheta = (double *)malloc(K*(K-1)*sizeof(double));
	KKx = (double *)malloc(K*(K - 1)*sizeof(double));
	KKy = (double *)malloc(K*(K - 1)*sizeof(double));
	//k3 = (double *)malloc(K*(K - 1)*sizeof(double));
	cg3 = (double *)malloc(K*(K - 1)*sizeof(double));
	D = (double *)malloc(K*(K - 1)*sizeof(double));
	Abnd = (double *)malloc(K*(K - 1)*sizeof(double));
	theta3 = (double *)malloc(K*(K - 1)*sizeof(double));
	dphi3 = (double *)malloc(K*(K - 1)*sizeof(double));
	double t1, t2, t2n, dift, chk1, chk2;
	std::complex<double> Comptemp, Comptemp2;


	//initiallise all variables 
	for (int n = 0; n < K*(K - 1); n++)
	{
		KKx[n] = 0.0;
		KKy[n] = 0.0;
		cg3[n] = 0.0;
		D[n] = 0.0;
		Abnd[n] = 0.0;
		theta3[n] = 0.0;
		dphi3[n] = 0.0;

	}




	//Run loop over wave-wave interaction components
	for (int i = 1; i < (K); i++)
	{
		// Determine difference frequency
		deltaf = (i+1)*dfgen;

		for (int m = 0; m < (K-i); m++)
		{
			int mi = i + m; // Cahnged from m+i+1
			//! Determine difference frequency
			deltatheta = abs(thetagen[mi] - thetagen[m]) + pi;

			//Determine x- and y-components of wave numbers of difference waves
			KKy[i + m*(K - 1)] = kgen[mi] * sin(thetagen[mi]) - kgen[m] * sin(thetagen[m]);
			KKx[i + m*(K - 1)] = kgen[mi] * cos(thetagen[mi]) - kgen[m] * cos(thetagen[m]);

			// Determine difference wave numbers according to Van Dongeren et al. 2003
			//	eq. 19
			k3 = sqrt(kgen[m] * kgen[m] + kgen[mi] * kgen[mi] + 2 * kgen[m] * kgen[mi] * cos(deltatheta));

			

			//Determine group velocity of difference waves

			cg3[i + m*(K - 1)] = 2 * pi*deltaf / k3;
			//Modification Robert + Jaap: make sure that the bound long wave amplitude does not
			//	!explode when offshore boundary is too close to shore,
			//	!by limiting the interaction group velocity
			cg3[i + m*(K - 1)] = min(cg3[i + m*(K - 1)], Param.nmax*sqrt(Param.g / k3*tanh(k3*Param.offdepth)));

			//Determine difference - interaction coefficient according to Herbers 1994
			//	!eq.A5
			t1 = (-1.0*wgen[m])*wgen[mi];
			t2 = (-1.0*wgen[m])+wgen[mi];
			t2n = cg3[i + m*(K - 1)] * k3;
			dift = abs(t2 - t2n);

			chk1 = cosh(kgen[m] * Param.offdepth);
			chk2 = cosh(kgen[mi] * Param.offdepth);

			D[i + m*(K - 1)] = -1.0*Param.g*kgen[m] * kgen[mi] * cos(deltatheta)*0.5 / t1 +
				Param.g*t2*(chk1*chk2) / ((Param.g*k3*tanh(k3*Param.offdepth) - t2n*t2n)*t1*cosh(k3*Param.offdepth))*
				(t2*(t1*t1 / (Param.g*Param.g) - kgen[m] * kgen[mi] * cos(deltatheta))
				- 0.5*((-1.0*wgen[m])*kgen[mi] * kgen[mi] / (chk2*chk2) + wgen[mi] * kgen[m] * kgen[m] / (chk1*chk1)));

			//Correct for surface elevation input and output instead of bottom pressure
			//	!so it is consistent with Van Dongeren et al 2003 eq. 18
			D[i + m*(K - 1)] = D[i + m*(K - 1)] * cosh(k3*Param.offdepth) / (cosh(kgen[m] * Param.offdepth)*cosh(kgen[mi] * Param.offdepth));

			// Exclude interactions with components smaller than or equal to current
			// component according to lower limit Herbers 1994 eq. 1
			if (fgen[m] <= deltaf)
			{
				D[i + m*(K - 1)] = 0.0;
			}

			// Exclude interactions with components that are cut - off by the fcutoff
			// parameter
			if (deltaf<=Param.fcutoff)
			{
				D[i + m*(K - 1)] = 0.0;
			}

			// Determine phase of bound long wave assuming a local equilibrium with
			// forcing of interacting primary waves according to Van Dongeren et al.
			// 2003 eq. 21 (the angle is the imaginary part of the natural log of a
			// complex number as long as the complex number is not zero)
			Comptemp = std::conj(CompFn(0, Findex[0] + mi  - 1)); //i or i+1 or m+i+1 or m+i ??
			Comptemp2 = std::conj(CompFn(0, Findex[0] + m - 1)); //m or m-1
			dphi3[i + m*(K - 1)] = pi + std::imag(log(Comptemp)) - std::imag(log(Comptemp2));

			// Determine angle of bound long wave according to Van Dongeren et al. 2003 eq. 22
			theta3[i + m*(K - 1)] = atan2(KKy[i + m*(K - 1)], KKx[i + m*(K - 1)]);

			// Determine energy of bound long wave according to Herbers 1994 eq. 1 based
			// on difference - interaction coefficient and energy density spectra of
			// primary waves
			// Robert: E = 2 * D**2 * S**2 * dtheta**2 * df can be rewritten as
			// E = 2 * D**2 * Sf**2 * df
			Eforc = 2 * D[i + m*(K - 1)] * D[i + m*(K - 1)] * vargenq[m] * vargenq[mi] * dfgen;
			Abnd[i + m*(K - 1)] = sqrt(2.0 * Eforc*dfgen);
			//if (!(Abnd[i + m*(K - 1)] == Abnd[i + m*(K - 1)])) // is nan
			//{
			//	printf("Arggg!");
			//}

		}
		
	}

	TwoDee<std::complex<double>> Ftempx(K-1, K);
	TwoDee<std::complex<double>> Ftempy(K-1, K);
	TwoDee<std::complex<double>> Ftemptot(K-1, K);

	// Run a loop over the offshore boundary
	for (int j = 0; j < Param.ny; j++)
	{
		//Reset Ftemp. necessary?
		for (int i = 1; i < K; i++)
		{

			for (int m = 0; m < (K); m++)
			{
				Ftempx(i-1, m) = 0.0;
				Ftempy(i-1, m) = 0.0;
				Ftemptot(i-1, m) = 0.0;
			}
		}
		
		for (int i = 1; i < K; i++)
		{

			for (int m = 0; m < (K - i); m++)
			{
				//qx
				Ftempx(i-1, m) = Abnd[i + m*(K - 1)] * 0.50 * exp(-1.0 * par_compi* dphi3[i + m*(K - 1)])*cg3[i + m*(K - 1)] * cos(theta3[i + m*(K - 1)]);
				
				if (Ftempx(i-1, m) != Ftempx(i-1, m))
				{
					printf("Error in generating Wave component q: #NAN or #IND detected");
					write_text_to_log_file("Error in generating  Wave component q: #NAN or #IND detected");
					exit(EXIT_FAILURE);
				}
				//qy
				Ftempy(i-1, m) = Abnd[i + m*(K - 1)] * 0.50 * exp(-1.0 * par_compi* dphi3[i + m*(K - 1)])*cg3[i + m*(K - 1)] * sin(theta3[i + m*(K - 1)]);
				//eta
				Ftemptot(i-1, m) = Abnd[i + m*(K - 1)] * 0.50 * exp(-1.0 * par_compi* dphi3[i + m*(K - 1)]);
			}
		}

		Gn = 0.0;
		tempcplxarr = 0.0;

		//! Unroll wave component to correct place along the offshore boundary
		for (int i = 1; i < K; i++)
		{

			for (int m = 0; m < (K); m++)
			{
				//qx
				Ftempx(i - 1, m) = Ftempx(i - 1, m) * exp(-1.0*par_compi*(KKy[i + m*(K - 1)] * (j*Param.dx) + KKx[i + m*(K - 1)] * (0.0*Param.dx)));
			}
		}

		for (int i = 1; i < K; i++)
		{

			std::complex<double> sume = 0.0;

			for (int m = 0; m < (K); m++)
			{
				sume = sume + Ftempx(i-1, m);
			}
			
			Gn[i] = sume; //[i] ot [i-1]
		}

		tempcplxarr = Gn[std::slice(1, tslen*0.5-1, 1)];
		for (int jcplx = 0; jcplx < tempcplxarr.size(); jcplx++)
		{
			tempcplxarr[jcplx] = std::conj(tempcplxarr[jcplx]);
		}
		flipiv(tempcplxarr);
		Gn[std::slice(tslen / 2 + 1, tempcplxarr.size(), 1)] = tempcplxarr;

		for (int n = 0; n < tslen; n++)
		{
			in[n][0] = std::real(Gn[n]);
			in[n][1] = std::imag(Gn[n]);
		}


		// Inverse Discrete Fourier transformation to transform back to time
		// domain from frequency domain
		fftw_execute(p);

		//fft function takes care of the scaling here?
		//store the results in qx
		for (int n = 0; n < tslen; n++)
		{
			qx[j + n*ny] = out[n][0] * taperf[n];

		}
		



		Gn = 0.0;
		tempcplxarr = 0.0;

		//! Unroll wave component to correct place along the offshore boundary
		for (int i = 1; i < K; i++)
		{

			for (int m = 0; m < (K); m++)
			{
				//qx
				Ftempy(i-1, m) = Ftempy(i-1, m) * exp(-1.0*par_compi*(KKy[i + m*(K-1)] * (j*Param.dx) + KKx[i + m*(K-1)] * (0.0*Param.dx)));
			}
		}
		for (int i = 1; i < K; i++)
		{

			std::complex<double> sume = 0.0;
			for (int m = 0; m < (K); m++)
			{
				sume = sume + Ftempy(i-1, m);
			}
			Gn[i] = sume;
		}

		tempcplxarr = Gn[std::slice(1, tslen*0.5 - 1, 1)];
		for (int jcplx = 0; jcplx < tempcplxarr.size(); jcplx++)
		{
			tempcplxarr[jcplx] = std::conj(tempcplxarr[jcplx]);
		}
		flipiv(tempcplxarr);
		Gn[std::slice(tslen / 2 + 1, tempcplxarr.size(), 1)] = tempcplxarr;

		for (int n = 0; n < tslen; n++)
		{
			in[n][0] = std::real(Gn[n]);
			in[n][1] = std::imag(Gn[n]);
		}


		// Inverse Discrete Fourier transformation to transform back to time
		// domain from frequency domain
		fftw_execute(p);

		//fft function takes care of the scaling here?
		//store the results in qi
		for (int n = 0; n < tslen; n++)
		{
			qy[j + n*ny] = out[n][0] * taperf[n];

		}

		Gn = 0.0;
		tempcplxarr = 0.0;

		//! Unroll wave component to correct place along the offshore boundary
		for (int i = 1; i < K; i++)
		{

			for (int m = 0; m < (K); m++)
			{
				//qx
				Ftemptot(i-1, m) = Ftemptot(i-1, m) * exp(-1.0*par_compi*(KKy[i + m*(K-1)] * (j*Param.dx) + KKx[i + m*(K-1)] * (0.0*Param.dx)));
			}
		}
		for (int i = 1; i < K; i++)
		{

			std::complex<double> sume = 0.0;
			for (int m = 0; m < (K); m++)
			{
				sume = sume + Ftemptot(i-1, m);
			}
			Gn[i] = sume;
		}

		tempcplxarr = Gn[std::slice(1, tslen*0.5 - 1, 1)];
		for (int jcplx = 0; jcplx < tempcplxarr.size(); jcplx++)
		{
			tempcplxarr[jcplx] = std::conj(tempcplxarr[jcplx]);
		}
		flipiv(tempcplxarr);
		Gn[std::slice(tslen / 2 + 1, tempcplxarr.size(), 1)] = tempcplxarr;

		for (int n = 0; n < tslen; n++)
		{
			in[n][0] = std::real(Gn[n]);
			in[n][1] = std::imag(Gn[n]);
		}


		// Inverse Discrete Fourier transformation to transform back to time
		// domain from frequency domain
		fftw_execute(p);

		//fft function takes care of the scaling here?
		//store the results in qx
		for (int n = 0; n < tslen; n++)
		{
			qtot[j + n*ny] = out[n][0] * taperf[n];

		}

	}

	
	
	//////////////////////////////////////
	// Generate q (qfile)
	//////////////////////////////////////

	//Until ee is fully tested leave as qfile=0.0
	for (int j = 0; j < Param.ny; j++)
	{

		for (int n = 0; n < tslen; n++)
		{
			qtempx[n] = qx[j + n*ny];
			qtempy[n] = qy[j + n*ny];
		}

		for (int m = 0; m < tslenbc; m++)
		{
			qfile[j + 0 * ny + m*ny * 4] = interp1DMono(tslen, tin, qtempx, m*Param.dtbc); 
			qfile[j + 1 * ny + m*ny * 4] = interp1DMono(tslen, tin, qtempy, m*Param.dtbc);
			qfile[j + 2 * ny + m*ny * 4] = 0.0;
			qfile[j + 3 * ny + m*ny * 4] = 0.0;// qtot[j + m*ny];
		}
		
	}

	

	//////////////////////////////////////
	//Clean up and desallocate all that memory
	//////////////////////////////////////

	free(Sf);
	free(Sd);
	free(pdf);
	free(cdf);

	free(fgen);
	free(phigen);
	free(thetagen);
	free(kgen);
	free(wgen);
	free(vargen);
	free(vargenq);
	free(Findex);
	free(WDindex);

	free(tin);
	free(zeta);
	//free(Ampzeta);//Reused zeta
	free(eta);
	free(Amp);
	free(stdzeta);
	free(E_tdir);
	free(qx);
	free(qy);
	free(qtot);

	free(taperf);
	free(taperw);

	free(binedgeleft);
	free(binedgeright);
	free(KKx);
	free(KKy);
	free(cg3);
	free(D);
	free(Abnd);
	free(theta3);
	free(dphi3);
	free(qtempx);
	free(qtempy);
	
	fftw_free(in);
	fftw_free(out);

}
