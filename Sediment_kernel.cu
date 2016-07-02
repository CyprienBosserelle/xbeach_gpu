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

#define pi 3.14159265

// declare texture reference for 2D DECNUM texture
//texture<DECNUM, 2, cudaReadModeElementType> texU;
//texture<DECNUM, 2, cudaReadModeElementType> texV;
//texture<DECNUM, 2, cudaReadModeElementType> texZ;

__global__ void longturb(int nx, int ny,DECNUM dx, DECNUM rho,DECNUM g,DECNUM dt,DECNUM beta,DECNUM * c,DECNUM *kturb,DECNUM * rolthick,DECNUM *dzsdt,DECNUM * uu,DECNUM *vv, DECNUM *hu, DECNUM *hv,int * wetu, int * wetv,DECNUM *h)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;


	
	
	__shared__ DECNUM  uui[16][16];
	__shared__ DECNUM  uul[16][16];
	__shared__ DECNUM  vvi[16][16];
	__shared__ DECNUM  vvb[16][16];
	__shared__ DECNUM kturbi[16][16];
	__shared__ DECNUM kturbl[16][16];
	__shared__ DECNUM kturbr[16][16];
	__shared__ DECNUM kturbb[16][16];
	__shared__ DECNUM kturbt[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		// use lagrangian velocities
		DECNUM kturbu = 0.0f;
		DECNUM kturbv = 0.0f;
		DECNUM dzsdt_cr = beta*c[i];
		DECNUM kkturb;
		DECNUM kturbumin, kturbvmin;
		DECNUM Sturbu, Sturbv, Sturbumin, Sturbvmin;
		DECNUM ksource, rolth;
		DECNUM betad = 1.0f;

		DECNUM hold = h[i] - dzsdt[i] * dt;

		kturbi[tx][ty] = kturb[i];
		kturbr[tx][ty] = kturb[xplus + iy*nx];
		kturbl[tx][ty] = kturb[xminus + iy*nx];
		kturbb[tx][ty] = kturb[ix + yminus*nx];
		kturbt[tx][ty] = kturb[ix + yplus*nx];

		uui[tx][ty] = uu[i];
		uul[tx][ty] = uu[xminus + iy*nx];

		vvi[tx][ty] = vv[i];
		vvb[tx][ty] = vv[ix + yminus*nx];


		// Update roller thickness
		rolth = rolthick[i] + dt*(abs(dzsdt[i]) - dzsdt_cr);
		rolthick[i] = max(rolth, 0.0f);

		//  X-direction
		kturbu = kturbi[tx][ty] * max(uui[tx][ty], 0.0f) + kturbr[tx][ty] * min(uui[tx][ty], 0.0f);
		/*if(uu[i]>0.0f)
		{
		kturbu=kturb[i];
		}
		else
		{
		if (uu[i]<0.0f)
		{     kturbu=kturb[xplus+iy*nx];}
		else
		{     kturbu=0.5f*(kturb[i]+kturb[xplus+iy*nx]);}
		}*/
		kturbumin = kturbl[tx][ty] * max(uul[tx][ty], 0.0f) + kturbi[tx][ty] * min(uul[tx][ty], 0.0f);
		/*if(uu[xminus+iy*nx]>0.0f)
		{
		kturbumin=kturb[xminus+iy*nx];
		}
		else
		{
		if(uu[xminus+iy*nx]<0.0f)
		{
		kturbumin=kturb[i];
		}
		else
		{
		kturbumin=0.5f*(kturb[xminus+iy*nx]+kturb[i]);
		}
		}*/


		Sturbu = kturbu*hu[i] * wetu[i];
		Sturbumin = kturbumin*hu[xminus + iy*nx] * wetu[xminus + iy*nx];


		// Y-direction
		kturbv = kturbi[tx][ty] * max(vvi[tx][ty], 0.0f) + kturbt[tx][ty] * min(vvi[tx][ty], 0.0f);

		/*if(vv[i]>0.0f)
		{
		kturbv=kturb[i];
		}
		else
		{
		if(vv[i]<0)
		{
		kturbv=kturb[ix+yplus*nx];
		}
		else
		{
		kturbv=0.5f*(kturb[i]+kturb[ix+yplus*nx]);
		}
		}*/
		kturbvmin = kturbb[tx][ty] * max(vvb[tx][ty], 0.0f) + kturbi[tx][ty] * min(vvb[tx][ty], 0.0f);
		/*if(vv[ix+yminus*nx]>0.0f)
		{
		kturbvmin=kturb[ix+yminus*nx];
		}
		else
		{
		if(vv[ix+yminus*nx]<0.0f)
		{
		kturbvmin=kturb[i];
		}
		else
		{
		kturbvmin=0.5f*(kturb[ix+yminus*nx]+kturb[i]);
		}
		}*/

		Sturbv = kturbv*hv[i] * wetv[i];
		Sturbvmin = kturbvmin*hv[ix + yminus*nx] * wetv[ix + yminus*nx];

		ksource = g*rolthick[i] * beta*c[i];     // only important in shallow water, where c=sqrt(gh)  

		kkturb = hold*kturb[i] - dt*((Sturbu - Sturbumin) / dx + (Sturbv - Sturbvmin) / dx - (ksource - betad*powf(kturb[i], 1.5f)));

		kturb[i] = 0.0f;//max(kkturb,0.0f);
	}


}

__global__ void Sbvr(int nx, int ny, DECNUM rho,DECNUM g,DECNUM eps, DECNUM Trep,DECNUM D50, DECNUM D90, DECNUM rhosed,DECNUM ws,DECNUM nuhfac,DECNUM * ueu, DECNUM * vev,DECNUM *H,DECNUM * DR,DECNUM * R, DECNUM * c,DECNUM * hh,DECNUM *urms,DECNUM * ceqsg,DECNUM * ceqbg, DECNUM *Tsg, DECNUM *zom, DECNUM * kturb)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;




	

	__shared__ DECNUM  hhi[16][16];
	//__shared__ DECNUM  Hi[16][16];
	__shared__ DECNUM  ueui[16][16];
	__shared__ DECNUM  ueul[16][16];
	__shared__ DECNUM  vevi[16][16];
	__shared__ DECNUM  vevb[16][16];
	
	DECNUM ue,ve;
	
	DECNUM vmags,vmag,ML,Tbore,dcfin,dcf,kb,Urms2;
	DECNUM B2,T1,Ucrc,Ucrw,Ucr,Ass,Asb,Cd,ceqb,ceqs;
	//DECNUM D50=0.0038;
	//DECNUM D90=0.0053;
	DECNUM zo=0.006f;//zom[i];
	DECNUM sedcal=1.0f;
	int wetz;
	DECNUM bulk=1.0f;//1.0f;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		//DECNUM rhosed=2500; //Sediment density
		DECNUM drho = (rhosed - rho) / rho;
		DECNUM dester = powf(drho*g, 1.0f / 3.0f) / 0.0001f*D50; //1.19e-4 comes from (Kb^2)^1/3 with Kb = 1.3e-6 m2s-2 kinematic viscosity of water
		// Calc euler velocities at cell center
		hhi[tx][ty] = max(hh[i], 0.01f);;
		ueui[tx][ty] = ueu[i];
		vevi[tx][ty] = vev[i];
		//Hi[tx][ty]=H[i];



		ueul[tx][ty] = ueu[xminus + iy*nx];
		vevb[tx][ty] = vev[ix + yminus*nx];
		__syncthreads;

		ue = 0.5*(ueul[tx][ty] + ueui[tx][ty]);
		ve = 0.5*(vevb[tx][ty] + vevi[tx][ty]);

		//need to check this...
		vmags = ue*ue + ve*ve;
		vmag = sqrt(vmags);

		//Mixing length
		ML = max(min(sqrt(2.0f*R[i] * Trep / (rho*c[i])), hhi[tx][ty]), 0.01f);

		//Bore period
		Tbore = Trep / 4.0f;// should be more complex //to improve later

		//Exponential decay of turbulence over time
		dcfin = expf(min(hhi[tx][ty], 100.0f) / ML);
		dcf = min(1 / (dcfin - 1), 1.0f);

		//Short wave turbulence (Breaking):
		kb = nuhfac*powf(DR[i] / rho, 0.66666666667f)*dcf;

		Urms2 = urms[i] * urms[i] + 1.45f*(kb + kturb[i]);//not been tested yet!!!



		DECNUM tsfac = 0.1f;
		//DECNUM ws=0.0509f;
		DECNUM Tsmin = 0.5f;
		Tsg[i] = max(tsfac*hhi[tx][ty] / ws, Tsmin); //should be different for each sediment class


		//critical U due to current
		//Ucrc=8.5f*pow(D50,0.6f)*log(4.0f*hhi[tx][ty]/D90)/log(10.0f);//Shields

		//Critical U due to Waves
		//Ucrw=0.95f*pow(1.65f*g,0.57f)*pow(D50,0.43f)*pow(Trep,0.14f);//Komar and Miller 1975

		//Critical velocity
		//Ucr=B2*Ucrc+(1.0f-B2)*Ucrw;
		if (D50 <= 0.0005f)
		{
			Ucr = 0.19f*powf(D50, 0.1f)*log10f(4.0f*hhi[tx][ty] / D90);
		}
		else
		{
			Ucr = 8.50f*powf(D50, 0.6f)*log10f(4.0f*hhi[tx][ty] / D90);

		}

		//drag coeff
		DECNUM hdrag = max(hhi[tx][ty], 10.0f*zo);
		Cd = 0.4f / (logf(hdrag / zo) - 1.0f);
		Cd = Cd*Cd;

		//Bottom sediment
		//Asb=0.005f*hhi[tx][ty]*powf(D50/hhi[tx][ty]/(drho*g*D50),1.2f);
		Asb = 0.005f*hhi[tx][ty] * powf(1 / hhi[tx][ty] / (drho*g), 1.2f);//simplified from above to limit the propagation of round of error with D50

		//Suspended Sediment
		Ass = 0.012f*D50*powf(dester, -0.6f) / (powf(drho*g*D50, 1.2f));



		//
		T1 = vmags + 0.018f / Cd*Urms2;
		T1 = min(T1, 100000.0f*g / zom[i] * D50*drho);
		T1 = sqrtf(T1);




		// Calculate Cequilibrium


		if (hhi[tx][ty] > eps)
		{
			wetz = 1;
		}
		else
		{
			wetz = 0;
		}


		DECNUM T2;
		T2 = 0.0f;

		if (T1 > Ucr && hhi[tx][ty] > eps)
		{
			T2 = powf(T1 - Ucr, 2.4f);
		}





		ceqb = Asb*T2;
		ceqb = min(ceqb / hhi[tx][ty], 0.05f);             //maximum equilibrium bed concentration

		// This should be different for each sediment fraction
		ceqbg[i] = (1 - bulk)*ceqb*sedcal*wetz;

		ceqs = min(Ass*T2 / hhi[tx][ty], 0.05f);// maximum equilibrium suspended concentration		      

		ceqsg[i] = (ceqs + bulk*ceqb)*sedcal*wetz;
	}
}


__global__ void Sednew(int nx, int ny, DECNUM rho,DECNUM g,DECNUM eps, DECNUM Trep,DECNUM D50, DECNUM D90, DECNUM rhosed,DECNUM ws,DECNUM nuhfac,DECNUM * ueu, DECNUM * vev,DECNUM *H,DECNUM * DR,DECNUM * R, DECNUM * c,DECNUM * hh,DECNUM *urms,DECNUM * ceqsg,DECNUM * ceqbg, DECNUM *Tsg, DECNUM *zom, DECNUM * kturb)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;




	

	__shared__ DECNUM  hhi[16][16];
	//__shared__ DECNUM  Hi[16][16];
	__shared__ DECNUM  ueui[16][16];
	__shared__ DECNUM  ueul[16][16];
	__shared__ DECNUM  vevi[16][16];
	__shared__ DECNUM  vevb[16][16];
	
	DECNUM ue,ve;
	
	DECNUM vmags,vmag,ML,Tbore,dcfin,dcf,kb,Urms2;
	DECNUM B2,T1,Ucrc,Ucrw,Ucr,Ass,Asb,Cd,ceqb,ceqs;
	//DECNUM D50=0.0038;
	//DECNUM D90=0.0053;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		DECNUM zo = zom[i];
		DECNUM sedcal = 1.0f;
		int wetz;
		DECNUM bulk = 1.0f;

		//DECNUM rhosed=25000; //Sediment density
		DECNUM drho = (rhosed - rho) / rho;
		DECNUM dester = powf(drho*g, 1.0f / 3.0f) / 0.0001*D50;
		// Calc euler velocities at cell center
		hhi[tx][ty] = hh[i];
		ueui[tx][ty] = ueu[i];
		vevi[tx][ty] = vev[i];
		//Hi[tx][ty]=H[i];



		ueul[tx][ty] = ueu[xminus + iy*nx];
		vevb[tx][ty] = vev[ix + yminus*nx];
		__syncthreads;

		ue = 0.5*(ueul[tx][ty] + ueui[tx][ty]);
		ve = 0.5*(vevb[tx][ty] + vevi[tx][ty]);

		//need to check this...
		vmags = ue*ue + ve*ve;
		vmag = sqrt(vmags);

		//Mixing length
		ML = max(min(sqrt(2 * R[i] * Trep / (rho*c[i])), hhi[tx][ty]), 0.01f);

		//Bore period
		Tbore = Trep / 4.0f;// should be more complex //to improve later

		//Exponential decay of turbulence over time
		dcfin = exp(min(hhi[tx][ty], 100.0f) / ML);
		dcf = min(1 / (dcfin - 1), 1.0f);

		//Short wave turbulence (Breaking):
		kb = nuhfac*powf(DR[i] / rho, 0.6666667f)*dcf;

		Urms2 = urms[i] * urms[i] + 1.45f*(kb + kturb[i]);//not been tested yet!!!
		//DECNUM dester=rhosed*D50;//dester=25296*D50;
		//DECNUM dster=(drho*g/1.0f-12)**onethird*s%D50(jg) 
		DECNUM tsfac = 0.1f;
		//DECNUM ws=0.043f;
		DECNUM Tsmin = 0.5f;
		Tsg[i] = max(tsfac*hhi[tx][ty] / ws, Tsmin); //should be different for each sediment class

		//DECNUM Ucrc,Ucrw;
		if (D50 <= 0.0005)
		{
			Ucrc = powf(0.19f*D50, 0.10f)*log10f(4.0f*hhi[tx][ty] / D90);
			Ucrw = powf(0.24f*drho*g, 0.66f)*powf(D50*Trep, 0.33);
		}
		if (D50<0.002 && D50>0.0005)
		{
			//critical U due to current
			Ucrc = 8.5f*powf(D50, 0.6f)*log(4.0f*hhi[tx][ty] / D90) / log(10.0f);//Shields

			//Critical U due to Waves
			Ucrw = 0.95f*powf(1.65f*g, 0.57f)*powf(D50, 0.43f)*powf(Trep, 0.14f);//Komar and Miller 1975
		}


		B2 = vmag / max(vmag + sqrtf(Urms2), eps);
		//Critical velocity
		Ucr = B2*Ucrc + (1.0f - B2)*Ucrw;

		//Bottom sediment
		Asb = 0.015f*hhi[tx][ty] * powf(D50 / hhi[tx][ty], 1.2f) / powf(drho*g*D50, 0.75f);

		//Suspended Sediment
		Ass = 0.012f*D50*powf(dester, -0.6f) / (powf(drho*g*D50, 1.2f));

		//
		T1 = vmags + 0.64f*Urms2;

		T1 = min(T1, 100000 * g / zom[i] * D50*drho);
		T1 = sqrtf(T1);


		// Calculate Cequilibrium


		if (hhi[tx][ty] > eps)
		{
			wetz = 1;
		}
		else
		{
			wetz = 0;
		}


		DECNUM T2;
		T2 = 0.0f;

		if (T1 > Ucr && hhi[tx][ty] > eps)
		{
			T2 = powf((T1 - Ucr), 1.5f);

			ceqb = Asb*T2;
			ceqb = min(ceqb / hhi[tx][ty], 0.05f);             //maximum equilibrium bed concentration
			// This should be different for each sediment fraction
			T2 = powf((T1 - Ucr), 2.4f);
			ceqs = min(Ass*T2 / hhi[tx][ty], 0.05f);// maximum equilibrium suspended concentration		      
		}

		ceqbg[i] = (1 - bulk)*ceqb*sedcal*wetz;
		ceqsg[i] = (ceqs + bulk*ceqb)*sedcal*wetz;
	}
}



__global__ void Rvr(int nx, int ny,DECNUM Trep,DECNUM facsk,DECNUM facas,DECNUM * H, DECNUM * hh, DECNUM * urms, DECNUM * c, DECNUM *ua)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;


	if (ix < nx && iy < ny)
	{
		// time averaged flows due to wave asymmetry

		DECNUM m1 = 0.0f;       // a = 0
		DECNUM m2 = 0.7939f;  // b = 0.79 +/- 0.023
		DECNUM m3 = -0.6065f; // c = -0.61 +/- 0.041
		DECNUM m4 = 0.3539f;  // d = -0.35 +/- 0.032 
		DECNUM m5 = 0.6373f;  // e = 0.64 +/- 0.025
		DECNUM m6 = 0.5995f;  // f = 0.60 +/- 0.043
		DECNUM alpha = -1.0f*log10(exp(1.0f)) / m4;
		DECNUM beta = exp(m3 / m4);
		DECNUM k = 2 * pi / (c[i] * Trep);

		DECNUM Ur, Bm, B1, Sk, As;
		//	if (abs(facua)>0.d0f)
		//	{
		Ur = 3.0f / 8.0f*sqrt(2.0f)*H[i] * k / powf(k*hh[i], 3.0f);                  //Ursell number
		Ur = max(Ur, 0.00001f);
		Bm = m1 + (m2 - m1) / (1.0f + beta*powf(Ur, alpha));                    //Boltzmann sigmoid (eq 6)         
		B1 = (-90.0f + 90.0f*tanh(m5 / powf(Ur, m6)))*pi / 180.0f;
		Sk = Bm*cos(B1);                                            //Skewness (eq 8)
		As = Bm*sin(B1);                                            //Asymmetry(eq 9)
		ua[i] = (facsk*Sk - facas*As)*urms[i];
		//	}
	}
}

__global__ void Erosus(int nx, int ny, DECNUM dt,DECNUM morfac,DECNUM por ,DECNUM * hh,DECNUM * ceqsg,DECNUM * ceqbg, DECNUM *Tsg, DECNUM * facero, DECNUM * structdepth)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	

	if (ix < nx && iy < ny)
	{
		//DECNUM morfac=0.0f; //Morphological factor 0.0= no bed update 1.0= normal bed update >1.0= enhance bed update
		//DECNUM por=0.4f;
		DECNUM pbbed = 1.0f;// sand fraction everywhere
		DECNUM exp_ero;
		//to be done for each sediment class


		exp_ero = morfac*dt / (1.0f - por)*hh[i] * (ceqsg[i] * pbbed / Tsg[i] + ceqbg[i] * pbbed / dt);
		facero[i] = min(1.0f, structdepth[i] * pbbed / max(0.000001f, exp_ero));        // limit erosion to available sediment on top 
	}
}


__global__ void Susp(int nx, int ny,DECNUM dx, DECNUM eps, DECNUM nuh,DECNUM nuhfac, DECNUM rho,DECNUM sus,DECNUM bed,DECNUM * ueu,DECNUM * vev,DECNUM * uu,DECNUM * uvg,DECNUM * hug,DECNUM * vv,DECNUM *vug,DECNUM *hvg,DECNUM * zb,DECNUM *h,DECNUM * DR, DECNUM * C,DECNUM * ceqbg,DECNUM * Sus, DECNUM * Svs,DECNUM * Sub, DECNUM * Svb,DECNUM * thetamean,DECNUM * ua)
{

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;


	DECNUM cu,cv,Dc,dcsdx,dcsdy,hu,hv;
	DECNUM dzbdx,dzbdy;
	DECNUM wetu,wetv;
	
	DECNUM vmagu, vmagv;
	DECNUM uau,uav;
	DECNUM uv,vu;
	DECNUM cub,cvb;
	
	//DECNUM sus=1.0f;
	//DECNUM bed=1.0f;
	
	DECNUM pbbed=1.0f; // WARNING sand fraction every where
	
	DECNUM facsl=1.6f; // between 0 and 1.6 tke into account the bed slope in bed load calculations
	DECNUM urep,vrep;



	
	__shared__ DECNUM  cci[16][16];
	__shared__ DECNUM  ccr[16][16];
	__shared__ DECNUM  cct[16][16];


	__shared__ DECNUM  cbi[16][16];
	__shared__ DECNUM  cbr[16][16];
	__shared__ DECNUM  cbt[16][16];

	__shared__ DECNUM hhi[16][16];
	
	__shared__ DECNUM zbi[16][16];
	__shared__ DECNUM zbr[16][16];
	__shared__ DECNUM zbt[16][16];
	
	__shared__ DECNUM uui[16][16];

	
	__shared__ DECNUM vvi[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		hhi[tx][ty] = h[i];


		zbi[tx][ty] = zb[i];
		zbr[tx][ty] = zb[xplus + iy*nx];
		zbt[tx][ty] = zb[ix + yplus*nx];


		cci[tx][ty] = C[i];
		ccr[tx][ty] = C[xplus + iy*nx];

		cct[tx][ty] = C[ix + yplus*nx];

		cbi[tx][ty] = pbbed*ceqbg[i];
		cbr[tx][ty] = pbbed*ceqbg[xplus + iy*nx];

		cbt[tx][ty] = pbbed*ceqbg[ix + yplus*nx];

		uui[tx][ty] = uu[i];


		vvi[tx][ty] = vv[i];



		__syncthreads;

		uau = 0.5*cosf(thetamean[i])*(ua[i] + ua[xplus + iy*nx]);
		uav = 0.5*sinf(thetamean[i])*(ua[xplus + iy*nx] + ua[i]);

		uv = uvg[i];//0.25f*(uul[tx][ty]+uui[tx][ty]+uutl[tx][ty]+uut[tx][ty]);
		vu = vug[i];//0.25f*(vvb[tx][ty]+vvi[tx][ty]+vvbr[tx][ty]+vvr[tx][ty]);
		urep = ueu[i] + uau;
		vmagu = sqrtf(powf((uui[tx][ty] + uau), 2.0f) + powf((vu + uav), 2.0f));

		uau = 0.5*cosf(thetamean[i])*(ua[i] + ua[ix + yplus*nx]);
		uav = 0.5*sinf(thetamean[i])*(ua[ix + yplus*nx] + ua[i]);

		vrep = vev[i] + uav;
		vmagv = sqrtf(powf(uv + uau, 2.0f) + powf(vvi[tx][ty] + uav, 2.0f));

		dzbdx = -1.0f*(zbr[tx][ty] - zbi[tx][ty]) / dx;
		dzbdy = -1.0f*(zbt[tx][ty] - zbi[tx][ty]) / dx;


		hu = hug[i];//0.50f*(hhi[tx][ty]+hhr[tx][ty]);
		hv = hvg[i];//0.50f*(hhi[tx][ty]+hht[tx][ty]);

		wetu = 0.0f;
		wetv = 0.0f;

		if (hu > eps)
		{
			wetu = 1.0f;
		}
		if (hv > eps)
		{
			wetv = 1.0f;
		}






		if (urep > 0.0f)
		{

			cu = cci[tx][ty];
			cub = cbi[tx][ty];
		}
		else
		{
			if (urep < 0.0f)
			{
				cu = ccr[tx][ty];
				cub = cbr[tx][ty];
			}
			else
			{
				cu = 0.50f*(cci[tx][ty] + ccr[tx][ty]);
				cub = 0.50f*(cbi[tx][ty] + cbr[tx][ty]);
			}
		}
		dcsdx = (ccr[tx][ty] - cci[tx][ty]) / dx;

		if (vrep > 0.0f)
		{

			cv = cci[tx][ty];
			cvb = cbi[tx][ty];
			//cvb(i,j)=par%thetanum*pbbed(i,j,1,jg)*ceqbg(i,j,jg)+(1.d0-par%thetanum)*pbbed(i,min(j+1,ny),1,jg)*ceqbg(i,min(j+1,ny),jg)
		}
		else
		{
			if (vrep < 0.0f)
			{
				cv = cct[tx][ty];
				cvb = cbt[tx][ty];
				//cvb(i,j)=par%thetanum*pbbed(i,j+1,1,jg)*ceqbg(i,j+1,jg)+(1.d0-par%thetanum)*pbbed(i,max(j,2),1,jg)*ceqbg(i,max(j,2),jg)
			}
			else
			{
				cv = 0.50f*(cci[tx][ty] + cct[tx][ty]);
				cvb = 0.50f*(cbi[tx][ty] + cbt[tx][ty]);
				//cvb(i,j)=0.5d0*(pbbed(i,j,1,jg)*ceqbg(i,j,jg)+pbbed(i,j+1,1,jg)*ceqbg(i,j+1,jg))

			}
		}

		dcsdy = (cct[tx][ty] - cci[tx][ty]) / dx;

		Dc = nuh + nuhfac*hhi[tx][ty] * powf(DR[i] / rho, 1.0f / 3.0f);

		Sus[i] = sus*(cu*urep*hu - Dc*hu*dcsdx - facsl*cu*vmagu*hu*dzbdx)*wetu;

		Svs[i] = sus*(cv*vrep*hv - Dc*hv*dcsdy - facsl*cv*vmagv*hv*dzbdy)*wetv;

		Sub[i] = bed*(cub*urep*hu - facsl*cub*vmagu*hu*dzbdx)*wetu;

		Svb[i] = bed*(cvb*vrep*hv - facsl*cvb*vmagv*hv*dzbdy)*wetv;
	}

}



__global__ void Conc(int nx, int ny, DECNUM dx, DECNUM dt,DECNUM eps,DECNUM * hh,DECNUM * C, DECNUM * ceqsg, DECNUM *Tsg,DECNUM *facero,DECNUM * ero,DECNUM * depo,DECNUM * Sus,DECNUM *Svs)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;



	__shared__ DECNUM  Susi[16][16];
	__shared__ DECNUM  Susl[16][16];
	__shared__ DECNUM  Svsi[16][16];
	__shared__ DECNUM  Svsb[16][16];
	__shared__ DECNUM  hhi[16][16];

	DECNUM cs,dsusdx,dsvsdy,wetz;
	DECNUM Pbed=1.0f;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		hhi[tx][ty] = hh[i];

		Susi[tx][ty] = Sus[i];
		Susl[tx][ty] = Sus[xminus + iy*nx];
		Svsi[tx][ty] = Svs[i];
		Svsb[tx][ty] = Svs[ix + yminus*nx];

		__syncthreads;

		wetz = 0.0f;
		if (hhi[tx][ty] > eps)
		{
			wetz = 1.0f;
		}


		ero[i] = facero[i] * hhi[tx][ty] * ceqsg[i] * Pbed / Tsg[i];

		dsusdx = (Susi[tx][ty] - Susl[tx][ty]) / dx;

		dsvsdy = (Svsi[tx][ty] - Svsb[tx][ty]) / dx;

		cs = (dt*Tsg[i]) / (dt + Tsg[i])*(hhi[tx][ty] * C[i] / dt - (dsusdx + dsvsdy - ero[i]))*wetz;
		cs = max(cs, 0.0f);
		cs = min(cs, 0.1f*hhi[tx][ty]);

		depo[i] = cs / Tsg[i];
		//cs=cs/hh[i];



		C[i] = cs / hh[i];
	}
}


__global__ void CClatbnd(int nx, int ny, DECNUM eps, DECNUM * hh, DECNUM * C)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;


	__shared__ DECNUM cci[16][16];
	__shared__ DECNUM cct[16][16];
	__shared__ DECNUM ccb[16][16];
	__shared__ DECNUM ccr[16][16];

	//cci[tx][ty]=C[i];
	//cct[tx][ty]=C[ix+yplus*nx];
	//ccb[tx][ty]=C[ix+yminus*nx];
	//ccr[tx][ty]=C[xplus+iy*nx];
	//__syncthreads;
	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		if (iy == 0)
		{
			C[i] = C[ix + yplus*nx];
		}

		if (iy == ny - 1)
		{
			C[i] = C[ix + yminus*nx];
		}
		if (ix == 0)
		{
			C[i] = 0.0f;//ccr[tx][ty];
		}
	}
}

__global__ void hardlayer(int nx, int ny,DECNUM dx,DECNUM dt,DECNUM * Sub, DECNUM * Svb, DECNUM * Sout, int * indSub,int * indSvb)
{

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;


	
	__shared__ DECNUM  Subi[16][16];
	__shared__ DECNUM  Subl[16][16];
	__shared__ DECNUM  Svbi[16][16];
	__shared__ DECNUM  Svbb[16][16];
	__shared__ DECNUM  Souti[16][16];
	
	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		Subi[tx][ty] = Sub[i];
		Subl[tx][ty] = Sub[xminus + iy*nx];
		Svbi[tx][ty] = Svb[i];
		Svbb[tx][ty] = Svb[ix + yminus*nx];
		Souti[tx][ty] = 0.0f;
		indSub[i] = 0;
		indSvb[i] = 0;

		if (Subi[tx][ty] > 0.0f) //      ! bed load u-direction
		{
			indSub[i] = 1;
			Souti[tx][ty] = Souti[tx][ty] + Subi[tx][ty] * dx;
		}
		if (Svbi[tx][ty] > 0.0f) //     ! bed load v-direction
		{
			indSvb[i] = 1;
			Souti[tx][ty] = Souti[tx][ty] + Svbi[tx][ty] * dx;
		}
		// fluxes at i-1,j
		if (Subl[tx][ty] < 0.0f) //   ! bed load u-direction
		{
			Souti[tx][ty] = Souti[tx][ty] - Subl[tx][ty] * dx;
		}
		// fluxes at i,j-1
		if (Svbb[tx][ty] < 0.0f) //   ! bed load v-direction
		{
			Souti[tx][ty] = Souti[tx][ty] - Svbb[tx][ty] * dx;
		}


		Sout[i] = Souti[tx][ty];
	}
}

__global__ void bedupdate(int nx, int ny, DECNUM eps, DECNUM dx, DECNUM dt, DECNUM morfac, DECNUM por, DECNUM * hh, DECNUM * ero, DECNUM * depo, DECNUM * Sub, DECNUM * Svb, DECNUM * Sout, int * indSub, int * indSvb, DECNUM * zb, DECNUM *ddzb, DECNUM * structdepth)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;





	__shared__ DECNUM  Subi[16][16];
	__shared__ DECNUM  Subl[16][16];
	__shared__ DECNUM  Svbi[16][16];
	__shared__ DECNUM  Svbb[16][16];

	__shared__ int indSubi[16][16];
	__shared__ int indSubl[16][16];
	__shared__ int indSvbi[16][16];
	__shared__ int indSvbb[16][16];

	if (ix < nx && iy < ny)
	{

		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		DECNUM oldzb = zb[i];
		DECNUM fac;
		DECNUM Savailable;
		DECNUM pbbed = 1.0f;

		Subi[tx][ty] = Sub[i];
		Subl[tx][ty] = Sub[xminus + iy*nx];
		Svbi[tx][ty] = Svb[i];
		Svbb[tx][ty] = Svb[ix + yminus*nx];

		indSubi[tx][ty] = indSub[i];
		indSubl[tx][ty] = indSub[xminus + iy*nx];
		indSvbi[tx][ty] = indSvb[i];
		indSvbb[tx][ty] = indSvb[ix + yminus*nx];


		__syncthreads;


		Savailable = structdepth[i] * pbbed / morfac / dt*(1.0f - por)*dx*dx;

		//	 ! reduction factor for cell outgoing sediment transports

		fac = 1.0f;
		if (Sout[i] > 0.0f)
		{
			fac = min(1.0f, Savailable / Sout[i]);
		}

		if (fac < 1.0f)
		{
			Subi[tx][ty] = fac*indSubi[tx][ty] * Subi[tx][ty] + (1 - indSubi[tx][ty])*Subi[tx][ty];
			Subl[tx][ty] = fac*(1 - indSubl[tx][ty])*Subl[tx][ty] + indSubl[tx][ty] * Subl[tx][ty];
			Svbi[tx][ty] = fac*indSvbi[tx][ty] * Svbi[tx][ty] + (1 - indSvbi[tx][ty])*Svbi[tx][ty];
			Svbb[tx][ty] = fac*(1 - indSvbb[tx][ty])*Svbb[tx][ty] + indSvbb[tx][ty] * Svbb[tx][ty];
		}




		DECNUM dzg;


		dzg = morfac*dt / (1.0f - por)*(ero[i] - depo[i] /*+ (Subi[tx][ty]-Subl[tx][ty])/dx + (Svbi[tx][ty]-Svbb[tx][ty])/dx*/);


		zb[i] = zb[i] + dzg;
		hh[i] = hh[i] + dzg;
		ddzb[i] = -1 * dzg;
		structdepth[i] = structdepth[i] - dzg;


	}
}





__global__ void zblatbnd(int nx,int ny,DECNUM * F)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;

	
	
	
	//__shared__ DECNUM Fi[16][16];
	//__shared__ DECNUM Ft[16][16];
	//__shared__ DECNUM Fb[16][16];
	//__shared__ DECNUM Fr[16][16];
	
	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		//Fi[tx][ty]=F[i];
		//Ft[tx][ty]=
		//Fb[tx][ty]=F[ix+yminus*nx];
		//Fr[tx][ty]=F[xplus+iy*nx];
		__syncthreads;

		//F[i]=Fi[tx][ty];
		if (iy == 0)
		{
			F[i] = F[ix + yplus*nx];
		}
		if (iy == ny - 1)
		{
			F[i] = F[ix + yminus*nx];
		}
		if (ix == 0)
		{
			F[i] = F[xplus + iy*nx];
		}

	}
			



}

__global__ void avalanching(int nx, int ny,DECNUM eps,DECNUM dx,DECNUM dt,DECNUM por,DECNUM drydzmax,DECNUM wetdzmax,DECNUM maxslpchg,DECNUM * hh,DECNUM * zb,DECNUM * dzb,DECNUM * structdepth)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;

	
	__shared__ DECNUM Zbi[16][16];
	__shared__ DECNUM Zbt[16][16];
	__shared__ DECNUM Zbr[16][16];
	
	__shared__ DECNUM dzbi[16][16];
	__shared__ DECNUM dzbt[16][16];
	__shared__ DECNUM dzbr[16][16];
	
	__shared__ DECNUM stdepi[16][16];
	__shared__ DECNUM stdept[16][16];
	__shared__ DECNUM stdepr[16][16];
	

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		Zbi[tx][ty] = zb[i];
		Zbt[tx][ty] = zb[ix + yplus*nx];
		Zbr[tx][ty] = zb[xplus + iy*nx];

		dzbi[tx][ty] = dzb[i];
		dzbt[tx][ty] = dzb[ix + yplus*nx];
		dzbr[tx][ty] = dzb[xplus + iy*nx];

		stdepi[tx][ty] = structdepth[i];
		stdept[tx][ty] = structdepth[ix + yplus*nx];
		stdepr[tx][ty] = structdepth[xplus + iy*nx];
		__syncthreads;

		DECNUM dzmaxdry = drydzmax;//1.0; // critical avalanching slope above water (dzbdx)
		DECNUM dzmaxwet = wetdzmax;//0.3; // critical avalanching slope under water
		DECNUM maxchg = maxslpchg;//0.05; // 0.05max bedlavel change due to Avalanching in m/s/m This avoid generatng tsunamis from avalanching
		DECNUM dzbdx, dzbdy, dzmax, dzbdxsign, dzbdysign;
		DECNUM dzbx = 0.0f;
		DECNUM dzby = 0.0f;

		if (hh[i] > eps)
		{
			dzmax = dzmaxwet;
		}
		else
		{
			dzmax = dzmaxdry;
		}


		// X direction
		dzbdx = (Zbr[tx][ty] - Zbi[tx][ty]) / dx;
		if (fabs(dzbdx) > dzmax)
		{
			dzbdxsign = dzbdx / fabs(dzbdx);
			dzbx = dzbdxsign*(fabs(dzbdx) - dzmax)*dx;
			if (dzbdxsign > 0)
			{
				dzbx = min(dzbx, maxchg*dt / dx);
				dzbx = min(dzbx, stdepi[tx][ty]);
			}
			else
			{
				dzbx = max(dzbx, -1.0f*maxchg*dt / dx);
				dzbx = max(dzbx, -1.0f*stdepr[tx][ty]);

			}
		}






		// Y direction
		dzbdy = (Zbt[tx][ty] - Zbi[tx][ty]) / dx;
		if (abs(dzbdy) > dzmax)
		{
			dzbdysign = dzbdy / abs(dzbdy);
			dzby = dzbdysign*(abs(dzbdy) - dzmax)*dx;
			if (dzbdysign > 0)
			{
				dzby = min(dzby, maxchg*dt / dx);
				dzby = min(dzby, stdepi[tx][ty]);
			}
			else
			{
				dzby = max(dzby, -1 * maxchg*dt / dx);
				dzby = max(dzby, -1 * stdept[tx][ty]);

			}
		}

		dzbi[tx][ty] = dzbi[tx][ty] + dzbx + dzby;
		dzbr[tx][ty] = dzbr[tx][ty] - dzbx;
		dzbt[tx][ty] = dzbt[tx][ty] - dzby;
		//__syncthreads;

		//Zb[i]=zbi[tx][ty]+dzbi[tx][ty]+dzbr[tx-1][ty]+dzbt[tx][ty-1];
		dzb[i] = dzb[i] + dzbx + dzby;
		dzb[xplus + iy*nx] = dzb[xplus + iy*nx] - dzbx;

		dzb[ix + yplus*nx] = dzb[ix + yplus*nx] - dzby;

	}


}

__global__ void updatezb(int nx,int ny,DECNUM dx,DECNUM dt,DECNUM * zb,DECNUM * ddzb,DECNUM * dzb,DECNUM * zs,DECNUM *hh, DECNUM * structdepth)
						
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	
	int tx =threadIdx.x;
	int ty= threadIdx.y;


	if (ix < nx && iy < ny)
	{

		zb[i] = zb[i] + ddzb[i];
		dzb[i] = dzb[i] - ddzb[i];
		zs[i] = zs[i] - ddzb[i];
		structdepth[i] = structdepth[i] - ddzb[i];
		hh[i] = hh[i] - ddzb[i];
	}
}


__global__ void updatezom(int nx, int ny,DECNUM cf,DECNUM cf2,DECNUM fw,DECNUM fw2,DECNUM * structdepth, DECNUM * cfm,DECNUM * fwm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	
	if (ix < nx && iy < ny)
	{

		if (structdepth[i] < 0.05f)
		{
			cfm[i] = cf2;
			fwm[i] = fw2;
		}
		else
		{
			cfm[i] = cf;
			fwm[i] = fw;
		}
	}
	
	
}



//__global__ void SedEnt(int nx, int ny,DECNUM dx,DECNUM dt, DECNUM rho,DECNUM g,DECNUM eps, DECNUM Trep,DECNUM * ueu, DECNUM * vev,DECNUM *H, DECNUM *DR,DECNUM *R,DECNUM *c,DECNUM * hh,DECNUM *urms,DECNUM *Sedup)
//{
//	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
//	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
//	unsigned int i=ix+iy*nx;
//	int tx =threadIdx.x;
//	int ty= threadIdx.y;
//
//
//	DECNUM Urms2,kb,teta,tau,Ab,fw,abkb,fent;
//
//	DECNUM D50=0.0038;
//	DECNUM D90=0.0090;
//	DECNUM zo=0.05;
//	DECNUM ggm=0.05;
//
//    DECNUM ws=0.043f;
//	DECNUM wp=0.00001f;
//	
//	DECNUM rhosed=25000; //Sediment density
//	DECNUM drho = (rhosed-rho)/rho;
//
//
//	//Short wave turbulence (Breaking):
//	kb=pow(DR[i]/rho,0.6666667f)*1/(exp(min(hh[i]/max(H[i],0.1f),100.0f))-1);
//	
//	Urms2=urms[i]*urms[i]+0.5f*kb;
//	Ab=urms[i]*Trep/(2*pi);
//	kb=30*zo;
//	abkb=Ab/kb;
//
//	if(abkb<=0.2f)
//	{
//		fw=0.3f;
//	}
//	else
//	{
//		if(abkb<=100)
//		{
//			fw=exp(-8.82+7.02*pow(abkb,-0.078f));
//		}
//		else
//		{
//			fw=exp(-7.30+5.61*pow(abkb,-0.109f));
//		}
//	}
//
//
//	tau=max(0.5*fw*rho*Urms2,(ueu[i]*ueu[i]+vev[i]*vev[i])*0.41f*0.41f/pow(log(hh[i]/(3*zo)),2));
//
//	teta=tau/((rhosed-rho)*g*D50);
//
//	fent=ggm*pow(teta,3)*ws;
//	
//	Sedup[i]=round((fent*dt*dx*dx)/wp);
//
//
//
//
//
//}
//
//__global__ void up3ddGPU(int nx, DECNUM *uu,DECNUM *vv,DECNUM *xx, DECNUM *yy, DECNUM *zz, DECNUM *dd_rand, DECNUM dx, DECNUM dt)
//{
//      DECNUM Dpx=10.0;
//      DECNUM Ux=0.05;
//      DECNUM Vx=0.05;
//     
//      //__shared__ DECNUM xxx[256];
//      //__shared__ DECNUM yyy[256];
//      DECNUM xxx,yyy;
//	  
// 
//      
// 
// 
//      
//      DECNUM Eh=0.001;//m2/s
//      DECNUM Ev=0.001;//m2/s
//      DECNUM ws=-0.043f;
//	  //DECNUM ws=-0.1;//m/s
//      int a=0;//abitrary number
//      //DECNUM dt=1;
//      DECNUM zo=0.001;//m roughness length
//      DECNUM ttc=0.01;// critical resuspension velocity m/s
//     
//      int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;
//      int tx = threadIdx.x;
//      int idx = threadIdx.x + blockIdx.x*blockDim.x;
// 
//      xxx=xx[i];
//      yyy=yy[i];
//
//      //DECNUM xp=xxx/dx;
//      //DECNUM yp=yyy/dx;
//      
//      //int x1=floor(xxx/dx);
//      //int y1=floor(yyy/dx);
//      
//      //int x2=x1+1;
//      //int y2=y1+1;
//
//      //DECNUM den=(x2-x1)*(y2-y1);
//      //DECNUM U11,U12,U21,U22,V11,V12,V21,V22;
//      
//      //U11=uu[x1+y1*nx];
//      //U21=uu[x2+y1*nx];
//      //U12=uu[x1+y2*nx];
//      //U22=uu[x2+y2*nx];
//
//      //V11=vv[x1+y1*nx];
//      //V21=vv[x2+y1*nx];
//      //V12=vv[x1+y2*nx];
//      //V22=vv[x2+y2*nx];
//
//
//
// 
//           
//      //Interpolate wter depth, Uvel Vvel at the particle position
//     
//      //Dpx=tex2D(texZ, xxx[tx]/dx, yyy[tx]/dx);
//     
//      Ux=tex2D(texU, xxx/dx, yyy/dx);
//      Vx=tex2D(texV, xxx/dx, yyy/dx);
//	      
//	//Ux=U11;//den*(x2-xp)*(y2-yp)+U21/den*(xp-x1)*(y2-yp)+U12/den*(x2-xp)*(yp-y1)+U22/den*(xp-x1)*(yp-y1);
//	//Vx=V11;//den*(x2-xp)*(y2-yp)+V21/den*(xp-x1)*(y2-yp)+V12/den*(x2-xp)*(yp-y1)+V22/den*(xp-x1)*(yp-y1);
//      //DECNUM T=9.81*1021*(Ux*Ux+Vx*Vx)/(18*log10(0.37*Dpx/zo));//Warning doggy equation!!!!!
//      
//     
//     
// 
//      //update the particle position
//     
//      //if (zz[i]*Dpx>zo)
//      //{
//         
//           
//      //Ux=Ux*(log10f(Dpx*zz[i]/zo)/log10f(0.37*Dpx/zo));
//      //Vx=Vx*(log10f(Dpx*zz[i]/zo)/log10f(0.37*Dpx/zo));
//     
// 
//      xx[i]=xxx+Ux*dt/*+(dd_rand[i]-0.5)*2*sqrtf(6*Eh*dt)*/;
//      yy[i]=yyy+Vx*dt/*+(dd_rand[i+a]-0.5)*2*sqrtf(6*Eh*dt)*/;
//      zz[i] =(zz[i]*Dpx+ws*dt/*+(dd_rand[i+2*a]-0.5)*2*sqrtf(6*Ev*dt)*/)/Dpx;
//
//      
//      
//      
//      
//     
// 
//            
// 
// 
//     
//
//
//}
//
