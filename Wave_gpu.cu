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

// includes, system

using DECNUM = float;

//#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
//#include <cutil.h>
#include <string>
#include <netcdf.h>
#include <time.h>



using namespace std;






//global variables

DECNUM Trep, Trepold, Trepnew;
DECNUM * St, *Stnew, *Stold;
double * Stfile;
double * qfile;
double * Tpfile;
int nwbndstep = 0;
int wavebndtype;
int nstep = 0;
int breakmod = 1;

DECNUM * hh;//=ones(20,40).*10; //Depth+SL
DECNUM * zb, *qbndold, *qbndnew;
DECNUM * qbndold_g, *qbndnew_g;
DECNUM * umeanbnd_g;
DECNUM * vmeanbnd_g;
DECNUM * umeanbnd;
DECNUM * hh_g, *uu_g, *vv_g, *zs_g, *zb_g, *hhold_g;
DECNUM * ueu_g, *vev_g;
DECNUM * vmageu_g, *vmagev_g;
DECNUM * uu;
DECNUM  *vv;
DECNUM *zs;
DECNUM *dummy;

DECNUM *xadvec_g, *yadvec_g, *thetaadvec_g;

DECNUM *hum_g, *hu_g;
DECNUM *hvm_g, *hv_g;
int *wetu_g, *wetv_g;
DECNUM *ududx_g, *vdudy_g;
DECNUM *udvdx_g, *vdvdy_g;
DECNUM *vu_g, *uv_g;
DECNUM  nuh, nuhfac;
DECNUM *nuh_g;
DECNUM * viscu_g, *viscv_g;

DECNUM uumin = 0.0f;


int nx, ny;
DECNUM dx, dt, eps;
DECNUM *arrmin, *arrmax;
DECNUM *arrmin_g, *arrmax_g;
DECNUM grdalpha;
double totaltime;
int nstpw, nwstp;//nb of hd step between wave step and next step for calculating waves and
DECNUM wdt;// wave model time step
DECNUM wavbndtime;
DECNUM slbndtime;
DECNUM windtime;
DECNUM Cd; //Wind drag
DECNUM fp, hm0gew, mainang, rt, scoeff, gam;
int nwavbnd, nwavfile;
DECNUM dtwavbnd;
int roller;
DECNUM wci;
DECNUM *wci_g;
DECNUM gammax, hwci, gammaa, n, alpha, beta, t1;
DECNUM fw, fw2;
DECNUM *fwm_g;//Wave dissipation factor map

DECNUM phi = (1.0f + sqrt(5.0f)) / 2;
DECNUM aphi = 1 / (phi + 1);
DECNUM bphi = phi / (phi + 1);
DECNUM twopi = 8 * atan(1.0f);

DECNUM g = 9.81f;
DECNUM rho = 1025.0f;
DECNUM zo;
DECNUM cf, cf2;//friction
DECNUM *cfm, *cfm_g; //friction map

DECNUM lat; //lattitude 
DECNUM fc; //coriolis

int ntheta;
DECNUM thetamin, thetamax;
DECNUM *theta;//=(0:100)*((pi)/100)+t1;see below
DECNUM *theta_g;

DECNUM * cgx, *cgy, *cx, *cy, *ctheta, *cxsth, *sxnth;//
DECNUM * cgx_g, *cgy_g, *cx_g, *cy_g, *ctheta_g, *cxsth_g, *sxnth_g, *eect_g;//

int var2plot = 2;// 1: wave height 2: eta 3: u 4: v
int colorindx;

DECNUM dang;
DECNUM dtheta;

DECNUM * ee;//=zeros(nx+1,ny+1,ntheta);
DECNUM * dd;
DECNUM * wete;

DECNUM * ee_g, *St_g;

DECNUM * tm_g;

DECNUM * rr;//;rr=zeros(nx+1,ny+1,ntheta);

DECNUM * drr;//=zeros(nx+1,ny+1,ntheta);
DECNUM * usd;//=zeros(nx+1,ny+1);
DECNUM * D;//=zeros(nx+1,ny+1);
DECNUM * D_g;
DECNUM * E, *H;
DECNUM * E_g, *H_g;
DECNUM * drr_g, *rr_g;
DECNUM * DR_g, *R_g;

DECNUM * DR, *R;

DECNUM * Sxx_g, *Syy_g, *Sxy_g, *Fx_g, *Fy_g;
DECNUM /** Sxx,* Syy,* Sxy,*/* Fx, *Fy;
DECNUM * thetamean;
DECNUM * thetamean_g;

DECNUM * urms_g;
DECNUM * ust_g;
DECNUM omega;// = 2*pi/Trep;

DECNUM D50, D90, rhosed;

DECNUM * kturb_g, *rolthick_g, *dzsdt_g;
DECNUM * Ceq_g, *ceqsg_g, *ceqbg_g, *Tsg_g, *facero_g;
DECNUM * C, *Cc_g, *stdep, *stdep_g;
DECNUM * Sus_g, *Svs_g, *Sub_g, *Svb_g;
DECNUM morfac, por;
DECNUM * ero_g, *depo_g;
DECNUM *dzb, *dzb_g, *ddzb_g;
DECNUM * Sout_g;
int * indSub_g, *indSvb_g;
DECNUM sus = 1.0f;
DECNUM bed = 1.0f;
DECNUM facsk, facas;

DECNUM * zsbnd;
DECNUM rtsl;
DECNUM zsbndnew, zsbndold;

DECNUM windth, windthold, windv, windvold, windvnew, windthnew, rtwind;
FILE * fsl;
FILE * fwav;
FILE * fwind;

FILE * Tsout;
//FILE * fXq,* fXE;

char tsoutfile[256];
int iout, jout;
int imodel; //1: Wave only; 2: Current Only; 3: Wave + current; 4: Wave + current + Sediment transport
int nstepout = 0;//nb of step between outputs. 0= no output
int istepout = 0;//
int nstepplot = 0;//nb of step between plots. 0= no plotting
int istepplot = 1;
int displayon = 0;
float endtime;



DECNUM wws;
DECNUM drydzmax;
DECNUM wetdzmax;
DECNUM maxslpchg;



DECNUM Hplotmax;

DECNUM * sigm, *thet, *costhet, *sinthet;
DECNUM * sigm_g, *thet_g, *costhet_g, *sinthet_g;
DECNUM * ueu, *vev, *urms;
DECNUM * ua_g;

DECNUM *k, *c, *kh, *cg, *sinh2kh;
DECNUM *k_g, *c_g, *kh_g, *cg_g, *sinh2kh_g;

DECNUM * dhdx, *dhdy, *dudx, *dudy, *dvdx, *dvdy;
DECNUM * dhdx_g, *dhdy_g, *dudx_g, *dudy_g, *dvdx_g, *dvdy_g;
DECNUM *dzsdx_g, *dzsdy_g;
DECNUM *zeros;

DECNUM * Hmean_g, *uumean_g, *vvmean_g, *hhmean_g, *zsmean_g, *Cmean_g;
DECNUM *dtflow_g;
DECNUM *Hmean, *uumean, *vvmean, *hhmean, *zsmean, *Cmean;

int GPUDEVICE;

int startflowstep;
int usesmago;

// Particle stuff
int npart;
int sedstart;
int wxstep = 1;


char wavebndfile[256];


#include "Wave_kernel.cu"
#include "Flow_kernel.cu"
#include "Sediment_kernel.cu"

//#include "read_input.cpp"

#include "XBeachGPU.h"
#include "Wavestep.cu"


template <class T> const T& min(const T& a, const T& b) {
	return !(b < a) ? a : b;     // or: return !comp(b,a)?a:b; for version (2)
}






/*void CUDA_CHECK(cudaError CUDerr)
{


if( cudaSuccess != CUDerr) {

fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \

__FILE__, __LINE__, cudaGetErrorString( CUDerr) );

exit(EXIT_FAILURE);

}
}*/


// Main loop that actually runs the model
void mainloopGPU(void)
{
	

	while (totaltime <= endtime)
	{
		dim3 blockDim(16, 16, 1);// This means that the grid has to be a factor of 16 on both x and y
		dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);

		dim3 blockDimLine(32, 1, 1);
		dim3 gridDimLine(ceil((nx*ny*1.0f) / blockDimLine.x), 1, 1);

		nstep++;
		wdt = dt; // previous timestep

		//Calculate timestep
		FLOWDT << <gridDim, blockDim, 0 >> >(nx, ny, dx, 0.8f, dtflow_g, hh_g);
		CUDA_CHECK(cudaThreadSynchronize());


		minmaxKernel << <gridDimLine, blockDimLine, 0 >> >(nx*ny, arrmax_g, arrmin_g, dtflow_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		finalminmaxKernel << <1, blockDimLine, 0 >> >(arrmax_g, arrmin_g);
		CUDA_CHECK(cudaThreadSynchronize());

		//CUDA_CHECK(cudaMemcpy(arrmax, arrmax_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(arrmin, arrmin_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));

		dt = arrmin[0]*0.5f;

		
		
		if ((imodel == 1 || imodel > 2) && totaltime>0.0)
		{
			float dtwave;
			// Make sure the CFL condition for flow do not violate CFL condition for Waves
			WAVEDT << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, 0.8f, dtheta, dtflow_g, ctheta_g);
			CUDA_CHECK(cudaThreadSynchronize());


			minmaxKernel << <gridDimLine, blockDimLine, 0 >> >(nx*ny, arrmax_g, arrmin_g, dtflow_g);
			//CUT_CHECK_ERROR("UpdateZom execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			finalminmaxKernel << <1, blockDimLine, 0 >> >(arrmax_g, arrmin_g);
			CUDA_CHECK(cudaThreadSynchronize());

			//CUDA_CHECK(cudaMemcpy(arrmax, arrmax_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(arrmin, arrmin_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));

			dtwave = arrmin[0]*0.5f;

			dt = min(dt, dtwave);
		}

		float diffdt = dt - wdt;
		// prevent time step to change too quickly (less than 10%)
		if (abs(diffdt) / wdt > 0.1)
		{
			dt = wdt*(1 + (diffdt) / abs(diffdt)*0.1);
		}


		//dt=0.18;

		//dt = wdt;

		



		//printf("Timestep: %f\n", dt);


		totaltime = totaltime + (double) dt;	//total run time acheived until now in s

		

		if (imodel == 1 || imodel > 2)
		{
			wavebnd(); // Calculate the boundary condition for this step
		}

		if (imodel >= 2)
		{
			flowbnd();// Calculate the flow boundary for this step
		}

		if (imodel == 1 || imodel > 2)
		{

			wavestep(); // Calculate the wave action ballance for this step
		}



		if (imodel >= 2)
		{
			flowstep();// solve the shallow water and continuity for this step
		}
		if (imodel >= 4 && nstep >= sedstart)
		{
			//Sediment step
			sedimentstep();//solve the sediment dispersion, and morphology
		}

		//add last value for avg calc
		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Hmean_g, H_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, uumean_g, uu_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, vvmean_g, vv_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, hhmean_g, hh_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, zsmean_g, zs_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Cmean_g, Cc_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());


		if (nstep == istepout && nstepout > 0)
		{
			istepout = istepout + nstepout;

			//Avg mean variables

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstepout, Hmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstepout, uumean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstepout, vvmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstepout, hhmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstepout, zsmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstepout, Cmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			// Download mean vars
			CUDA_CHECK(cudaMemcpy(Hmean, Hmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(uumean, uumean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(vvmean, vvmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(hhmean, hhmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(zsmean, zsmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(Cmean, Cmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));


			CUDA_CHECK(cudaMemcpy(H, H_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(uu, uu_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(vv, vv_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(zs, zs_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(Fx, Fx_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(Fy, Fy_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(thetamean, thetamean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(D, D_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(urms, urms_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(ueu, ueu_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(vev, vev_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK( cudaMemcpy(C, ceqsg_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			CUDA_CHECK(cudaMemcpy(C, hum_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK( cudaMemcpy(C,k_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			//CUDA_CHECK( cudaMemcpy(ctheta,ee_g, nx*ny*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			CUDA_CHECK(cudaMemcpy(hh, hh_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			if (imodel == 4)// If moprhology is on
			{
				CUDA_CHECK(cudaMemcpy(zb, zb_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
				CUDA_CHECK(cudaMemcpy(dzb, dzb_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			}
			//CUDA_CHECK( cudaMemcpy(xxp, xxp_g, npart*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			//CUDA_CHECK( cudaMemcpy(yyp, yyp_g, npart*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			printf("Writing output, totaltime:%f s\n", totaltime);
			writestep2nc(tsoutfile, nx, ny,/*npart,*/(float) totaltime, imodel,/*xxp,yyp,*/zb, zs, uu, vv, H, H, thetamean, D, urms, ueu, vev, C, dzb, Fx, Fy, hh, Hmean, uumean, vvmean, hhmean, zsmean, Cmean);

			//write3dvarnc(nx,ny,ntheta,totaltime,ctheta);
			//outfile[],nx,ny,npart,totaltime,xxp,yyp,zs,uu, vv, H,Tp,Dp,      D,Urms,ueu,vev)
			//fprintf(Tsout,"%f\t%f\t%f\t%f\t%f\t%f\n",totaltime,hh[iout+jout*nx],zs[iout+jout*nx],uu[iout+jout*nx],vv[iout+jout*nx],H[iout+jout*nx]);


			//Clear avg vars
			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Hmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, uumean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, vvmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, hhmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, zsmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Cmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());
		}
	}
}


void mainloopCPU(void)
{
	printf("Computing CPU mode\n");
	while (totaltime <= endtime)
	{

		nstep++;
		wdt = dt; // Sometinmes in stationary wave run one can have a larger wdt (wave time step)
		totaltime = totaltime + dt;	//total run time acheived until now in s



		if (imodel == 1 || imodel > 2)
		{
			wavebnd(); // Calculate the boundary condition for this step
		}

		if (imodel >= 2)
		{
			flowbnd();// Calculate the flow boundary for this step
		}
		if (imodel == 1 || imodel > 2)
		{
			wavestepCPU(); // Calculate the wave action ballance for this step
		}
		if (imodel >= 2)
		{
			//flowstepCPU();// solve the shallow water and continuity for this step
		}
		if (imodel >= 4 && nstep >= sedstart)
		{
			//Sediment step
			//sedimentstepCPU();//solve the sediment dispersion, and morphology
		}

		//add last value for avg calc
		addavg_varCPU(nx, ny, Hmean_g, H_g);
		addavg_varCPU(nx, ny, uumean_g, uu_g);
		addavg_varCPU(nx, ny, vvmean_g, vv_g);
		addavg_varCPU(nx, ny, hhmean_g, hh_g);
		addavg_varCPU(nx, ny, zsmean_g, zs_g);
		addavg_varCPU(nx, ny, Cmean_g, Cc_g);

		if (nstep == istepout && nstepout > 0)
		{
			istepout = istepout + nstepout;

			//Avg mean variables

			divavg_varCPU(nx, ny, nstepout, Hmean_g);
			divavg_varCPU(nx, ny, nstepout, uumean_g);
			divavg_varCPU(nx, ny, nstepout, vvmean_g);
			divavg_varCPU(nx, ny, nstepout, hhmean_g);
			divavg_varCPU(nx, ny, nstepout, zsmean_g);
			divavg_varCPU(nx, ny, nstepout, Cmean_g);

			printf("Writing output, totaltime:%d s\n", totaltime);
			//printf("test Hs: %f\n",H_g[0+16*nx]);
			writestep2nc(tsoutfile, nx, ny, (float) totaltime, imodel, zb_g, zs_g, uu_g, vv_g, H_g, xadvec_g, thetamean_g, D_g, urms_g, ueu_g, vev_g, Cc_g, dzb_g, Fx_g, Fy_g, hh_g, Hmean_g, uumean_g, vvmean_g, hhmean_g, zsmean_g, Cmean_g);

			//Clear avg vars
			resetavg_varCPU(nx, ny, Hmean_g);
			resetavg_varCPU(nx, ny, uumean_g);
			resetavg_varCPU(nx, ny, vvmean_g);
			resetavg_varCPU(nx, ny, hhmean_g);
			resetavg_varCPU(nx, ny, zsmean_g);
			resetavg_varCPU(nx, ny, Cmean_g);
		}

	}

}




void flowbnd(void)
{
	//update sl bnd

	if (totaltime >= slbndtime)
	{

		zsbndold = zsbndnew;
		rtsl = slbndtime;
		fscanf(fsl, "%f\t%f", &slbndtime, &zsbndnew);
		//slbndtime=+rtsl;
		//zsbnd=zsbndold+(t-slbndtime+rtsl)*(zsbndnew-zsbndold)/rtsl;
	}





	if (wavebndtype == 1)
	{
		for (int ni = 0; ni < ny; ni++)
		{
			zsbnd[ni] = zsbndold + ((float) totaltime - rtsl)*(zsbndnew - zsbndold) / (slbndtime - rtsl);
		}
	}

	if (wavebndtype == 2)
	{
		if (GPUDEVICE >= 0)
		{
			dim3 blockDim(16, 16, 1);
			dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);
			// FLow abs_2d should be here not at the flow step		
			// Set weakly reflective offshore boundary
			ubnd1D << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, g, rho, (float) totaltime, wavbndtime, dtwavbnd, slbndtime, rtsl, zsbndold, zsbndnew, Trep, qbndold_g, qbndnew_g, zs_g, uu_g, vv_g, vu_g, umeanbnd_g, vmeanbnd_g, zb_g, cg_g, hum_g, cfm_g, Fx_g, hh_g);
			//CUT_CHECK_ERROR("ubnd execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());
		}
		else
		{
			ubndCPU(nx, ny, dx, dt, g, rho, (float) totaltime, wavbndtime, dtwavbnd, slbndtime, rtsl, zsbndold, zsbndnew, Trep, qbndold_g, qbndnew_g, zs_g, uu_g, vv_g, vu_g, umeanbnd_g, vmeanbnd_g, zb_g, cg_g, hum_g, cfm_g, Fx_g, hh_g);

		}
	}




	if (totaltime >= windtime)
	{
		windthold = windthnew;
		windvold = windvnew;
		rtwind = windtime;
		fscanf(fwind, "%f\t%f\t%f", &windtime, &windvnew, &windthnew);
		//windtime=windtime+rtwind;
		//printf("windthold=%f\n",windthold);
		//printf("windthnew=%f\n",windthnew);
	}
	windth = windthold + (totaltime - rtwind)*(windthnew - windthold) / (windtime - rtwind);
	windv = windvold + (totaltime - rtwind)*(windvnew - windvold) / (windtime - rtwind);
	//printf("windv=%f\n",windv);

	windth = (1.5*pi - grdalpha) - windth*pi / 180;
	//printf("windv=%f\twindth=%f\n",windv,windth);





}


void flowstep(void)
{
	// Flow model timestep
	dim3 blockDim(16, 16, 1);
	dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);

	dim3 blockDim4(4, 4, 1);
	dim3 gridDim4(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);


	//////////////////////////////////////////
	// BELOW IS FOR DEBUGGING ONLY
	/////////////////////////////////////////
	//read2Dnc(nx,ny,"Hfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(H_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"Fxfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(Fx_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"Fyfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(Fy_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"urmsfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(urms_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"ustfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(ust_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"thetameanfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(thetamean_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"uufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(uu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vvfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vv_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"zsfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(zs_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hhfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hh_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"ueufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(ueu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vevfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vev_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );



	//read2Dnc(nx,ny,"hufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"humfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hum_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hvfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hv_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hvmfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hvm_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vmageufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vmageu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vmagevfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vmagev_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );


	// Set weakly reflective offshore boundary ! MOVED TO FLOW BND SUBROUTINE!!
	//ubnd<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,g,rho,totaltime,wavbndtime,dtwavbnd,slbndtime,rtsl,zsbndold,zsbndnew,Trep,qbndold_g,qbndnew_g,zs_g,uu_g,vv_g,vu_g,umeanbnd_g,vmeanbnd_g,zb_g,cg_g,hum_g,cfm_g,Fx_g,hh_g);
	//CUT_CHECK_ERROR("ubnd execution failed\n");
	//CUDA_CHECK( cudaThreadSynchronize() );


	//
	// Water level slopes
	//
	wlevslopes << <gridDim, blockDim, 0 >> >(nx, ny, dx, eps, zs_g, dzsdx_g, dzsdy_g, hh_g);
	//CUT_CHECK_ERROR("wlevslopes execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//
	// Water depth at u pts for momentum and continuity eq (hum hu)
	//

	udepthmomcont << <gridDim, blockDim, 0 >> >(nx, ny, dx, eps, uumin, wetu_g, zs_g, uu_g, hh_g, hum_g, hu_g, zb_g);
	//CUT_CHECK_ERROR("udepthmomcont execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//
	// Water depth at v pts for momentum and continuity eq (hvm hv)
	//

	vdepthmomcont << <gridDim, blockDim, 0 >> >(nx, ny, dx, eps, uumin, wetv_g, zs_g, vv_g, hh_g, hvm_g, hv_g, zb_g);
	//CUT_CHECK_ERROR("vdepthmomcont execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());






	//
	// Advection in the x direction using 2n order finite difference
	//

	ududx_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, dx, hu_g, hum_g, uu_g, ududx_g);
	//CUT_CHECK_ERROR("uadvec execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//vdudy
	vdudy_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, dx, hv_g, hum_g, uu_g, vv_g, vdudy_g);
	//CUT_CHECK_ERROR("uadvec execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());




	//
	// Smagorinsky formulation or Normal eddy viscosity
	//
	CUDA_CHECK(cudaMalloc((void **)&nuh_g, nx*ny*sizeof(DECNUM)));
	smago << <gridDim, blockDim, 0 >> >(nx, ny, dx, uu_g, vv_g, nuh, nuh_g, usesmago);
	//CUT_CHECK_ERROR("uadvec execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//
	// increase eddy viscosity by wave induced breaking as in Reniers 2004 & Set viscu = 0.0 near water line
	//
	CUDA_CHECK(cudaMalloc((void **)&viscu_g, nx*ny*sizeof(DECNUM)));
	viscou << <gridDim, blockDim, 0 >> >(nx, ny, dx, rho, eps, nuhfac, nuh_g, hh_g, hum_g, hvm_g, DR_g, uu_g, wetu_g, viscu_g);
	//CUT_CHECK_ERROR("visco execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Explicit Euler step momentum u-direction
	//

	eulerustep << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, g, rho, cfm_g, fc, windth, windv, Cd, uu_g, urms_g, ududx_g, vdudy_g, viscu_g, dzsdx_g, hu_g, hum_g, Fx_g, vu_g, ueu_g, vmageu_g, wetu_g);
	//CUT_CHECK_ERROR("eulerustep execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd << <gridDim, blockDim, 0 >> >(nx, ny, uu_g, vv_g, zs_g);
	//CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Advection in the y direction using 2n order finite difference
	//
	//vdvdy
	vdvdy_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, dx, hv_g, hvm_g, vv_g, vdvdy_g);
	//CUT_CHECK_ERROR("vadvec for v execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());
	//udvdx

	udvdx_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, dx, hu_g, hvm_g, uu_g, vv_g, udvdx_g);
	//CUT_CHECK_ERROR("vadvec for v execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//
	// increase eddy viscosity by wave induced breaking as in Reniers 2004 & Set viscv = 0.0 near water line
	//
	CUDA_CHECK(cudaMalloc((void **)&viscv_g, nx*ny*sizeof(DECNUM)));
	viscov << <gridDim, blockDim, 0 >> >(nx, ny, dx, rho, eps, nuhfac, nuh_g, hh_g, hum_g, hvm_g, DR_g, vv_g, wetv_g, viscv_g);
	//CUT_CHECK_ERROR("visco v execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaFree(nuh_g));


	viscovbnd << <gridDim, blockDim, 0 >> >(nx, ny, viscv_g);
	//CUT_CHECK_ERROR("visco v execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Explicit Euler step momentum v-direction
	//

	eulervstep << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, g, rho, cfm_g, fc, windth, windv, Cd, vv_g, urms_g, udvdx_g, vdvdy_g, viscv_g, dzsdy_g, hv_g, hvm_g, Fy_g, uv_g, vev_g, vmagev_g, wetv_g);
	//CUT_CHECK_ERROR("eulervstep execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd << <gridDim, blockDim, 0 >> >(nx, ny, uu_g, vv_g, zs_g);
	//CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//
	//v velocities at u pts and u velocities at v pts
	//

	calcuvvu << <gridDim, blockDim, 0 >> >(nx, ny, dx, uu_g, vv_g, vu_g, uv_g, ust_g, thetamean_g, ueu_g, vev_g, vmageu_g, vmagev_g, wetu_g, wetv_g);
	//CUT_CHECK_ERROR("calcuvvu execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	uvlatbnd << <gridDim, blockDim, 0 >> >(nx, ny, vu_g, uv_g, ueu_g, vev_g, vmageu_g, vmagev_g);
	//fix side bnd for vu
	//twodimbndnoix<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,vu_g);
	//CUT_CHECK_ERROR("wave force X bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());









	//
	//Calculate hu
	//
	depthhu << <gridDim, blockDim, 0 >> >(nx, ny, dx, uumin, eps, hh_g, uu_g, hu_g, zs_g, zb_g);
	//CUT_CHECK_ERROR("depthhu execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//
	//Calculate hv
	//
	depthhv << <gridDim, blockDim, 0 >> >(nx, ny, dx, uumin, eps, hh_g, vv_g, hv_g, zs_g, zb_g);
	//CUT_CHECK_ERROR("depthhv execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Update water level using continuity eq.
	//
	continuity << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, eps, uu_g, hu_g, vv_g, hv_g, zs_g, hh_g, zb_g, dzsdt_g);
	//CUT_CHECK_ERROR("continuity execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd << <gridDim, blockDim, 0 >> >(nx, ny, uu_g, vv_g, zs_g);
	//CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	hsbnd << <gridDim, blockDim, 0 >> >(nx, ny, eps, hh_g, zb_g, zs_g);
	//CUT_CHECK_ERROR("hh lateral bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	CUDA_CHECK(cudaFree(viscu_g));
	CUDA_CHECK(cudaFree(viscv_g));


}

void sedimentstep(void)
{
	// suspended sediment timestep
	dim3 blockDim(16, 16, 1);
	dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);

	dim3 blockDim4(4, 4, 1);
	dim3 gridDim4(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);

	/////////////////////////////////////////////////////
	// BELOW IS FOR DEBUGGING ONLY
	/////////////////////////////////////////////////////

	//read2Dnc(nx,ny,"uufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(uu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vvfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vv_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hvfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hv_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );


	//read2Dnc(nx,ny,"ueufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(ueu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );


	//read2Dnc(nx,ny,"vevfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vev_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"Hfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(H_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hhfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hh_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"urmsfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(urms_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"uvfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(uv_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );


	//
	// Compute long wave turbulence due to breaking
	//
	longturb << <gridDim, blockDim, 0 >> >(nx, ny, dx, rho, g, dt, beta, c_g, kturb_g, rolthick_g, dzsdt_g, uu_g, vv_g, hu_g, hv_g, wetu_g, wetv_g, hh_g);
	//CUT_CHECK_ERROR("longturb execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//
	// Calculate Equilibrium concentration Ceq
	//
	//CUDA_CHECK( cudaMalloc((void **)&ceqsg_g, nx*ny*sizeof(DECNUM )) );
	CUDA_CHECK(cudaMalloc((void **)&ceqbg_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Tsg_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&ua_g, nx*ny*sizeof(DECNUM)));
	//BEWARE BELOW SHOULD BE hh_old_g
	//Sbvr or Sednew
	Sbvr << <gridDim, blockDim, 0 >> >(nx, ny, rho, g, eps, Trep, D50, D90, rhosed, wws, nuhfac, ueu_g, vev_g, H_g, DR_g, R_g, c_g, hh_g, urms_g, ceqsg_g, ceqbg_g, Tsg_g, cfm_g, kturb_g);
	//CUT_CHECK_ERROR("CalcCeq execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	Rvr << <gridDim, blockDim, 0 >> >(nx, ny, Trep, facsk, facas, H_g, hh_g, urms_g, c_g, ua_g);
	//CUT_CHECK_ERROR("Rvr execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	// CUDA_CHECK( cudaMemcpy(uumean, ua_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );




	//
	// Limit erosion to available sediment on top of hard layer
	//
	CUDA_CHECK(cudaMalloc((void **)&facero_g, nx*ny*sizeof(DECNUM)));
	Erosus << <gridDim, blockDim, 0 >> >(nx, ny, dt, morfac, por, hh_g, ceqsg_g, ceqbg_g, Tsg_g, facero_g, stdep_g);
	//CUT_CHECK_ERROR("Erosus execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//////////////////
	// LOOP for sediment fractions should come here
	/////////////////////

	//
	// suspended load 
	// 
	CUDA_CHECK(cudaMalloc((void **)&Sus_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Svs_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Sub_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Svb_g, nx*ny*sizeof(DECNUM)));
	Susp << <gridDim, blockDim, 0 >> >(nx, ny, dx, eps, nuh, nuhfac, rho, sus, bed, ueu_g, vev_g, uu_g, uv_g, hu_g, vv_g, vu_g, hv_g, zb_g, hh_g, DR_g, Cc_g, ceqbg_g, Sus_g, Svs_g, Sub_g, Svb_g, thetamean_g, ua_g);
	//CUT_CHECK_ERROR("Susp execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	//Calculate suspended concentration
	//

	CUDA_CHECK(cudaMalloc((void **)&ero_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&depo_g, nx*ny*sizeof(DECNUM)));
	Conc << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, eps, hh_g, Cc_g, ceqsg_g, Tsg_g, facero_g, ero_g, depo_g, Sus_g, Svs_g);
	//CUT_CHECK_ERROR("Conc execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//
	// Update global variables and fix bnds
	//
	CClatbnd << <gridDim, blockDim, 0 >> >(nx, ny, eps, hh_g, Cc_g);
	//CUT_CHECK_ERROR("CClatbnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	if (morfac > 0.0f)// Only if morphology is need i.e. if imodel=4 and morphac >0.0
	{
		//
		// Adjust sediment fluxes for rocklayer
		//

		CUDA_CHECK(cudaMalloc((void **)&Sout_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&indSub_g, nx*ny*sizeof(int)));
		CUDA_CHECK(cudaMalloc((void **)&indSvb_g, nx*ny*sizeof(int)));
		hardlayer << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, Sub_g, Svb_g, Sout_g, indSub_g, indSvb_g);
		//CUT_CHECK_ERROR("hardlayer execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, Sout_g);
		//CUT_CHECK_ERROR("Sout twodimbnd execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());


		//
		// Bed update
		//
		bedupdate << <gridDim, blockDim, 0 >> >(nx, ny, eps, dx, dt, morfac, por, hh_g, ero_g, depo_g, Sub_g, Svb_g, Sout_g, indSub_g, indSvb_g, zb_g, dzb_g, stdep_g);
		//CUT_CHECK_ERROR("bedupdate execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());


		//
		// Update lateral bnd	
		//	
		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, zb_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());


		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, stdep_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());


		//
		// Avalanching
		//
		CUDA_CHECK(cudaMalloc((void **)&ddzb_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMemcpy(ddzb_g, zeros, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		avalanching << <gridDim, blockDim, 0 >> >(nx, ny, eps, dx, dt, por, drydzmax, wetdzmax, maxslpchg, hh_g, zb_g, ddzb_g, stdep_g);
		//CUT_CHECK_ERROR("avalanching execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		//
		// Update Zb for avalanching
		//

		updatezb << <gridDim, blockDim, 0 >> >(nx, ny, dx, dt, zb_g, ddzb_g, dzb_g, zs_g, hh_g, stdep_g);
		//CUT_CHECK_ERROR("avalanching execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		//
		// Update lateral bnd	
		//	
		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, zb_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		//zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,dzb_g);
		////CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		//CUDA_CHECK( cudaThreadSynchronize() );

		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, stdep_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		updatezom << <gridDim, blockDim, 0 >> >(nx, ny, cf, cf2, fw, fw2, stdep_g, cfm_g, fwm_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());
	}




	///////////////////
	// END LOOP sediment fraction

	CUDA_CHECK(cudaFree(ddzb_g));
	CUDA_CHECK(cudaFree(Sout_g));
	CUDA_CHECK(cudaFree(indSub_g));
	CUDA_CHECK(cudaFree(indSvb_g));
	CUDA_CHECK(cudaFree(ua_g));


	CUDA_CHECK(cudaFree(ero_g));
	CUDA_CHECK(cudaFree(depo_g));
	CUDA_CHECK(cudaFree(Sus_g));
	CUDA_CHECK(cudaFree(Svs_g));
	CUDA_CHECK(cudaFree(Sub_g));
	CUDA_CHECK(cudaFree(Svb_g));
	CUDA_CHECK(cudaFree(facero_g));
	//CUDA_CHECK( cudaFree(ceqsg_g));
	CUDA_CHECK(cudaFree(ceqbg_g));
	CUDA_CHECK(cudaFree(Tsg_g));



}



int main(int argc, char **argv)
{

	// Start timer to keep track of time 
	clock_t startcputime, endcputime;


	startcputime = clock();

	// Initialise totaltime
	totaltime = 0.0f;

	//////////////////////////////////////////////////////
	/////             Read Operational file          /////
	//////////////////////////////////////////////////////


	char opfile[] = "opfile.dat"; // Compulsory input file

	char filename[256];

	char slbnd[256];
	char windfile[256];
	char zofile[256];
	char HLfile[256];




	FILE * fop;
	fop = fopen(opfile, "r");
	fscanf(fop, "%*s");//Dummy string
	fscanf(fop, "%s\t%*s", &filename);// Bathy file name needs to be md format
	fscanf(fop, "%d\t%*s", &imodel);// Type of model: 1: wave only; 2: currents only 3: waves+currents 4:waves+currents+sediment(+ morphology if morfac>0)
	fscanf(fop, "%d\t%*s", &GPUDEVICE);// What GPU device to use 
	fscanf(fop, "%f\t%*s", &dt);// Model time step in s. //This should be calculated by the model
	fscanf(fop, "%d\t%*s", &nstpw);// Number of flow and sediment step between wave step needs to be 1 for unsteady runs 
	fscanf(fop, "%f\t%*s", &eps);//drying height in m
	fscanf(fop, "%f,%f\t%*s", &cf, &cf2);// bottom friction for flow model cf is for sand and cf2 fro reef area (Reef and sand discrimination is done based on structure file if none is present cf2 should not be used )
	//fscanf(fop,"%f\t%*s",&cf);
	//fscanf(fop,"%s\t%*s",&zofile);
	fscanf(fop, "%f\t%*s", &nuh); // Viscosity coeff or samgo coeff depennding on usesmago
	fscanf(fop, "%f\t%*s", &nuhfac);//nuhfac=1.0f;//0.001f; //viscosity coefficient for roller induced turbulent horizontal viscosity// it should be small contrary to what XBeach recommend as default
	fscanf(fop, "%d\t%*s", &usesmago);// Uses smagorynsky formulation to calculate viscosity 0: No 1: Yes
	fscanf(fop, "%f\t%*s", &lat);// Latitude of the grid use negative for south hemisphere (this implies the grid is small on earth scale)
	fscanf(fop, "%f\t%*s", &Cd);// Wind drag coeff
	fscanf(fop, "%f\t%*s", &wci); // Wave current interaction switch (can also be used as a number between 0 and 1 to reduce the interaction if unstable)
	fscanf(fop, "%f\t%*s", &hwci);// hwci=0.010f;//min depth for wci
	fscanf(fop, "%d\t%*s", &breakmod); // Wave dissipation model 1: roelvink 2: Baldock. use 1 for unsteady runs (i.e. with wave group) and use 2 for steady runs
	fscanf(fop, "%f\t%*s", &gammaa);// Wave breaking gamma param 
	fscanf(fop, "%f\t%*s", &n);// exponential; in Roelving breaking model
	fscanf(fop, "%f\t%*s", &alpha);// calibration for wave dissipation (should be 1)
	fscanf(fop, "%f\t%*s", &gammax);//gammax=2.0f; //maximum ratio Hrms/hh
	fscanf(fop, "%f\t%*s", &beta);// Roller slope dissipation param
	fscanf(fop, "%f,%f\t%*s", &fw, &fw2);// Wave bottom dissipation parameters fw is for sand fw2 is for reefs. see cf comments
	//fscanf(fop,"%f\t%*s",&fw);
	fscanf(fop, "%f,%f\t%*s", &D50, &D90);// sand grain size in m
	//fscanf(fop,"%f\t%*s",&D50);
	//fscanf(fop,"%f\t%*s",&D90);
	fscanf(fop, "%f\t%*s", &rhosed);// sand density
	fscanf(fop, "%f\t%*s", &wws);// sand fall velocity (should be calculated)
	fscanf(fop, "%f,%f\t%*s", &drydzmax, &wetdzmax);// max slope in avalannching model
	//fscanf(fop,"%f\t%*s",&drydzmax);
	//fscanf(fop,"%f\t%*s",&wetdzmax);
	fscanf(fop, "%f\t%*s", &maxslpchg);// max change within a step to avoid avalanching tsunami
	fscanf(fop, "%f\t%*s", &por);// sand porosity (should not be constant)
	fscanf(fop, "%f\t%*s", &morfac);// morphological factor 0 no changes in morphology 1 normal changes in morpho >1 accelerated morphological changes (beware this doesn't accelerate the bnd you have to do this manually)
	fscanf(fop, "%f,%f\t%*s", &sus, &bed);// calibration coeff for suspended load and bed load
	//fscanf(fop,"%f\t%*s",&sus);
	//fscanf(fop,"%f\t%*s",&bed);
	fscanf(fop, "%f,%f\t%*s", &facsk, &facas);// calibration factor for wave skewness and Asymetry
	//fscanf(fop,"%f\t%*s",&facsk);
	//fscanf(fop,"%f\t%*s",&facas);
	fscanf(fop, "%s\t%*s", &HLfile);// Structure file write down "none" if none present
	fscanf(fop, "%d\t%*s", &wavebndtype); // 1 is quasistationary wave spectrum; 2 is for infrgravity and long bound waves Xbeach type
	fscanf(fop, "%s\t%*s", &wavebndfile);// wave bnd file see wiki for details
	fscanf(fop, "%s\t%*s", &slbnd); // tide/surge bnd file
	fscanf(fop, "%s\t%*s", &windfile); // Wind forcing file
	//fscanf(fop,"%d\t%*s",&npart);
	fscanf(fop, "%d\t%*s", &sedstart);// which step to start sediment transport and morpho
	//fscanf(fop,"%f\t%*s",&Hplotmax);
	//fscanf(fop,"%d\t%*s",&nstepplot);
	fscanf(fop, "%d\t%*s", &nstepout); // output step
	fscanf(fop, "%f\t%*s", &endtime);// end step
	//fscanf(fop,"%d\t%d\t%*s",&iout,&jout);
	fscanf(fop, "%s\t%*s", &tsoutfile);// output file
	fclose(fop);

	printf("bathy file: %s\n", filename);
	printf("Imodel: %d\n", imodel);
	//printf("nstepplot: %d\n",nstepplot);
	printf("smago?: %d\n", usesmago);
	printf("bed: %f\n", bed);

	printf("facsk=%f facas=%f\n", facsk, facsk);



	wdt = 0.0;

	FILE * fid;
	FILE * fiz;


	//read input data:
	printf("bathy: %s\n", filename);

	fid = fopen(filename, "r");
	fscanf(fid, "%u\t%u\t%f\t%*f\t%f", &nx, &ny, &dx, &grdalpha);
	printf("nx=%d\tny=%d\tdx=%f\talpha=%f\n", nx, ny, dx, grdalpha);


	//READ INITIAL ZS CONDITION
	//fiz=fopen("zsinit.md","r");
	//fscanf(fiz,"%u\t%u\t%f\t%*f\t%f",&nx,&ny,&dx,&grdalpha);

	grdalpha = grdalpha*pi / 180; // grid rotation

	//if(imodel>=2)
	//{
	printf("Opening sea level bnd...");
	fsl = fopen(slbnd, "r");

	fscanf(fsl, "%f\t%f", &rtsl, &zsbndold);

	//Note: the first rtsl should be 0 
	fscanf(fsl, "%f\t%f", &slbndtime, &zsbndnew);
	printf("done\n");

	//zsbnd sea level in bnd file
	//rtsl bnd time
	//slbndtime=rtsl;
	//}
	//else
	//{
	//	zsbndold=0;
	//	zsbndnew=0;
	//}

	// Allocate CPU memory

	hh = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	uu = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	vv = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	zs = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	zb = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	cfm = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	dzb = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	stdep = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	zeros = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	umeanbnd = (DECNUM *)malloc(ny*sizeof(DECNUM));
	



	// set initital condition and read bathy file
	printf("Set initial condition...");

	int jread;
	//int jreadzs;
	for (int fnod = ny; fnod >= 1; fnod--)
	{

		fscanf(fid, "%u", &jread);
		//fscanf(fiz,"%u",&jreadzs);
		umeanbnd[(jread - 1)] = 0.0f;
		for (int inod = 0; inod < nx; inod++)
		{
			fscanf(fid, "%f", &zb[inod + (jread - 1)*nx]);
			uu[inod + (jread - 1)*nx] = 0.0f;
			vv[inod + (jread - 1)*nx] = 0.0f;
			dzb[inod + (jread - 1)*nx] = 0.0f;
			cfm[inod + (jread - 1)*nx] = cf;
			stdep[inod + (jread - 1)*nx] = 0.0f;
			zeros[inod + (jread - 1)*nx] = 0.0f;
			//fscanf(fiz,"%f",&zs[inod+(jreadzs-1)*nx]);

			//hh[inod+(jread-1)*nx]=max(zb[inod+(jread-1)*nx]+zs[inod+(jreadzs-1)*nx],eps);
			//zs[inod+(jread-1)*nx]=max(zs[inod+(jreadzs-1)*nx],-1*zb[inod+(jread-1)*nx]);

			zs[inod + (jread - 1)*nx] = max(zsbndold, -1 * zb[inod + (jread - 1)*nx]);
			hh[inod + (jread - 1)*nx] = max(zb[inod + (jread - 1)*nx] + zsbndold, eps);



		}
	}

	fclose(fid);
	printf("...done\n");
	//fclose(fiz);
	char nofrictionfile[] = "none";

	// Friction file is now obsolete but may be needed agin in the future 
	/*	////Friction file


		int testfrictfile=strcmp(zofile,nofrictionfile);
		if (testfrictfile!=0)
		{
		printf("Friction file found\n");
		fid=fopen(zofile,"r");
		int jx,jy;
		fscanf(fid,"%u\t%u\t%*f\t%*f\t%*f",&jx,&jy);

		if (jx!=nx || jy!=ny)
		{
		printf("Error friction file dimension mismatch. Model will run with constant friction.\n");
		}
		}
		else
		{
		printf("No friction file found\n");
		}

		for (int fnod=ny; fnod>=1;fnod--)
		{
		if(testfrictfile!=0)
		{fscanf(fid,"%u",&jread);}
		for(int inod=0; inod<nx; inod++)
		{
		if(testfrictfile!=0)
		{
		//fscanf(fid,"%f",&zom[inod+(jread-1)*nx]);
		cfm[inod+(fnod-1)*nx]=zo;
		}
		else
		{
		cfm[inod+(fnod-1)*nx]=zo;
		}

		}
		}
		if(testfrictfile!=0)
		{
		fclose(fid);
		}

		*/
	//// read Hard layer file


	int testfrictfile = strcmp(HLfile, nofrictionfile);
	if (testfrictfile != 0)
	{
		printf("Hard layer file found\n");
		fid = fopen(HLfile, "r");
		int jx, jy;
		fscanf(fid, "%u\t%u\t%*f\t%*f\t%*f", &jx, &jy, &dx, &grdalpha);

		if (jx != nx || jy != ny)
		{
			printf("Error Hard layer file dimension mismatch. Model will run with constant friction.\n");
		}
	}
	else
	{
		printf("No hard layer file found\n");
	}

	for (int fnod = ny; fnod >= 1; fnod--)
	{
		if (testfrictfile != 0)
		{
			fscanf(fid, "%u", &jread);
		}
		for (int inod = 0; inod < nx; inod++)
		{
			if (testfrictfile != 0)
			{
				fscanf(fid, "%f", &stdep[inod + (jread - 1)*nx]);
			}
			else
			{
				stdep[inod + (fnod - 1)*nx] = 5.0f;

			}

		}
	}
	if (testfrictfile != 0)
	{
		fclose(fid);
	}





	// Read Wind forcing
	printf("Opening wind bnd\n");
	fwind = fopen(windfile, "r");
	fscanf(fwind, "%f\t%f\t%f", &rtwind, &windvold, &windthold);
	fscanf(fwind, "%f\t%f\t%f", &windtime, &windvnew, &windthnew);





	windv = windvold;
	windth = (1.5f*pi - grdalpha) - windthold*pi / 180.0f;



	//zo=0.1;//roughness length
	//cf=zo;//zo;
	//lat=-35.0;
	//calculate coriolis force
	lat = lat*pi / 180.0f;
	DECNUM wearth = pi*(1.0f / 24.0f) / 1800.0f;
	fc = 2.0f*wearth*sin(lat);


	//gammax=2.0f; //maximum ratio Hrms/hh
	//wci=1.0f; //switch for wave/current interaction.
	//hwci=0.010f;//min depth for wci
	//gammaa = 0.55; //breaker parameter in Baldock or Roelvink formulation
	//n=10.0;// power in roelvink dissipation model
	//alpha = 1.;//! wave dissipation coefficient
	//beta=0.15f;//! breaker slope coefficient in roller model
	roller = 1; // option to turn off/on roller model (0/1)

	//nuh=0.05; // Eddy viscosity [m2/s]
	//nuhfac=1.0f;//0.001f; //viscosity coefficient for roller induced turbulent horizontal viscosity

	t1 = -(pi) / 2;
	//thetamin=-60;
	//thetamax=60;

	//dtheta=10;

	//Allocate More array on CPU

	Fx = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	Fy = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	zsbnd = (DECNUM *)malloc(ny*sizeof(DECNUM));

	cgx = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cgy = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cx = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cy = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	ctheta = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));







	//Sxx=(DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	//Syy=(DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	//Sxy=(DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	usd = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	D = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	E = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	H = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	urms = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	ueu = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	vev = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	thetamean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));


	Hmean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	uumean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	vvmean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	hhmean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	zsmean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	Cmean = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	arrmax = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	arrmin = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));



	// Set initial water level on offshore bnd
	//Unnecessary
	//for (int ii = 0; ii < nx; ii++)
	//{
	//	for (int jj = 0; jj < ny; jj++)
	//	{
	//		zs[ii + jj*nx] = max(zsbndold, -1.0f*zb[ii + jj*nx]+eps);
	//		hh[ii + jj*nx] = zb[ii + jj*nx] + zs[ii + jj*nx];
	//	}
//
	//}
	// Allocate more CPU memory


	omega = 2 * pi / Trep;

	sigm = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	//sigt= (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	thet = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	//costhet=(DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	//sinthet=(DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));





	k = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	c = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	kh = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	cg = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	sinh2kh = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));

	dhdx = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	dhdy = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	dudx = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	dudy = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	dvdx = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	dvdy = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));

	C = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	R = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	DR = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));



	//ONLY FOR GPU BELOW CPU


	if (GPUDEVICE >= 0)
	{
		// Init GPU

		CUDA_CHECK(cudaSetDevice(GPUDEVICE));


		if (imodel == 1 || imodel > 2)
		{
			waveinitGPU();
		}

		//CUT_DEVICE_INIT(argc, argv);
		printf("Allocating GPU memory\n");




		//Allocate GPU memory
		CUDA_CHECK(cudaMalloc((void **)&hh_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&uu_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vv_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&wci_g, nx*ny*sizeof(DECNUM)));

		CUDA_CHECK(cudaMalloc((void **)&ueu_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vev_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vmageu_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vmagev_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&zs_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&zb_g, nx*ny*sizeof(DECNUM)));
		//CUDA_CHECK( cudaMalloc((void **)&Ceq_g, nx*ny*sizeof(DECNUM )) );
		CUDA_CHECK(cudaMalloc((void **)&dzsdx_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&dzsdy_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&cfm_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&fwm_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&dzsdt_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&hu_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&hv_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&hum_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&hvm_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vu_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&uv_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&wetu_g, nx*ny*sizeof(int)));
		CUDA_CHECK(cudaMalloc((void **)&wetv_g, nx*ny*sizeof(int)));
		CUDA_CHECK(cudaMalloc((void **)&ududx_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vdudy_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vdvdy_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&udvdx_g, nx*ny*sizeof(DECNUM)));



		CUDA_CHECK(cudaMalloc((void **)&D_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&urms_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&ust_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&Fx_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&Fy_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&k_g, nx*ny*sizeof(DECNUM)));


		CUDA_CHECK(cudaMalloc((void **)&ee_g, nx*ny*ntheta*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&rr_g, nx*ny*ntheta*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&St_g, ny*ntheta*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&sigm_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&DR_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&R_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&H_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&qbndold_g, 4 * ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&qbndnew_g, 4 * ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&umeanbnd_g, ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vmeanbnd_g, ny*sizeof(DECNUM)));
		//CUDA_CHECK( cudaMalloc((void **)&sigt_g, nx*ny*ntheta*sizeof(DECNUM )) );
		CUDA_CHECK(cudaMalloc((void **)&c_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&cg_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&theta_g, ntheta*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&cxsth_g, ntheta*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&sxnth_g, ntheta*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&thetamean_g, nx*ny*sizeof(DECNUM)));
		//CUDA_CHECK( cudaMalloc((void **)&Sxx_g, nx*ny*sizeof(DECNUM )) );
		//CUDA_CHECK( cudaMalloc((void **)&Syy_g, nx*ny*sizeof(DECNUM )) );
		//CUDA_CHECK( cudaMalloc((void **)&Sxy_g, nx*ny*sizeof(DECNUM )) );
		CUDA_CHECK(cudaMalloc((void **)&ctheta_g, nx*ny*ntheta*sizeof(DECNUM)));



		CUDA_CHECK(cudaMalloc((void **)&stdep_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&Cc_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&dzb_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&kturb_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&rolthick_g, nx*ny*sizeof(DECNUM)));





		// Averaged variables
		CUDA_CHECK(cudaMalloc((void **)&Hmean_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&uumean_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&vvmean_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&hhmean_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&zsmean_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&Cmean_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&ceqsg_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&arrmin_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&arrmax_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&dtflow_g, nx*ny*sizeof(DECNUM)));


		

	}
	else
	{
		// Waveinit CPU!!

		if (imodel == 1 || imodel > 2)
		{
			waveinitGPU();
		}

		//Allocate GPU memory

		//wci_g=(DECNUM *)malloc(nx*ny*sizeof(DECNUM));;

		hh_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		uu_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vv_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		wci_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));

		ueu_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vev_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vmageu_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vmagev_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		zs_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		zb_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));


		//CUDA_CHECK( cudaMalloc((void **)&Ceq_g, nx*ny*sizeof(DECNUM )) );
		dzsdx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dzsdy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		cfm_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		fwm_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dzsdt_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		hu_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		hv_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		hum_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		hvm_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vu_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		uv_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		wetu_g = (int *)malloc(nx*ny*sizeof(int));
		wetv_g = (int *)malloc(nx*ny*sizeof(int));
		ududx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vdudy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vdvdy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		udvdx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));



		D_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		urms_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		ust_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Fx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Fy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		k_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));


		ee_g = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
		rr_g = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
		St_g = (DECNUM *)malloc(ny*ntheta*sizeof(DECNUM));
		sigm_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		DR_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		R_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		H_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		qbndold_g = (DECNUM *)malloc(4 * ny*sizeof(DECNUM));
		qbndnew_g = (DECNUM *)malloc(4 * ny*sizeof(DECNUM));
		umeanbnd_g = (DECNUM *)malloc(ny*sizeof(DECNUM));
		vmeanbnd_g = (DECNUM *)malloc(ny*sizeof(DECNUM));
		//CUDA_CHECK( cudaMalloc((void **)&sigt_g, nx*ny*ntheta*sizeof(DECNUM )) );
		c_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		cg_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		theta_g = (DECNUM *)malloc(ntheta*sizeof(DECNUM));
		cxsth_g = (DECNUM *)malloc(ntheta*sizeof(DECNUM));
		sxnth_g = (DECNUM *)malloc(ntheta*sizeof(DECNUM));
		thetamean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		//CUDA_CHECK( cudaMalloc((void **)&Sxx_g, nx*ny*sizeof(DECNUM )) );
		//CUDA_CHECK( cudaMalloc((void **)&Syy_g, nx*ny*sizeof(DECNUM )) );
		//CUDA_CHECK( cudaMalloc((void **)&Sxy_g, nx*ny*sizeof(DECNUM )) );
		ctheta_g = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));



		stdep_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Cc_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dzb_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		kturb_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		rolthick_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));





		// Averaged variables
		Hmean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		uumean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		vvmean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		hhmean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		zsmean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Cmean_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		ceqsg_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		kh_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		sinh2kh_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dhdx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dhdy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dudx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dudy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dvdx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		dvdy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		xadvec_g = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
		yadvec_g = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
		thetaadvec_g = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
		E_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Sxx_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Sxy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
		Syy_g = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));


	}

	if (GPUDEVICE >= 0)
	{

		////////////////////////////////////////////////////////////////////////////////////////////
		// Copy CPU array to the GPU                                                       /////////
		////////////////////////////////////////////////////////////////////////////////////////////



		CUDA_CHECK(cudaMemcpy(hh_g, hh, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(hu_g, hh, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(hv_g, hh, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(hum_g, hh, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(hvm_g, hh, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(zb_g, zb, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(zs_g, zs, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(uu_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vv_g, vv, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		//CUDA_CHECK( cudaMemcpy(zom_g, zom, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );
		CUDA_CHECK(cudaMemcpy(vu_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(uv_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vev_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(ueu_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vmageu_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vmagev_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(ududx_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vdudy_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vdvdy_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(udvdx_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));


		CUDA_CHECK(cudaMemcpy(H_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Fx_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Fy_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(urms_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(ust_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(D_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));



		CUDA_CHECK(cudaMemcpy(ee_g, ee, nx*ny*ntheta*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(rr_g, rr, nx*ny*ntheta*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(cxsth_g, cxsth, ntheta*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(sxnth_g, sxnth, ntheta*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(theta_g, theta, ntheta*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(thetamean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(R_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(DR_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(c_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(cg_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));



		CUDA_CHECK(cudaMemcpy(stdep_g, stdep, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Cc_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(kturb_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(rolthick_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(dzb_g, dzb, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));

		CUDA_CHECK(cudaMemcpy(umeanbnd_g, umeanbnd, ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vmeanbnd_g, umeanbnd, ny*sizeof(DECNUM), cudaMemcpyHostToDevice));

		//Averaged variables
		CUDA_CHECK(cudaMemcpy(Hmean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(uumean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(vvmean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(hhmean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(zsmean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Cmean_g, uu, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));


		/*CUDA_CHECK( cudaMemcpy(xxp_g, xxp, npart*sizeof(DECNUM ), cudaMemcpyHostToDevice) );
		CUDA_CHECK( cudaMemcpy(yyp_g, yyp, npart*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

		*/

	}
	else
	{
		for (int jj = 0; jj < ny; jj++)
		{
			umeanbnd_g[jj] = umeanbnd[jj];
			vmeanbnd_g[jj] = umeanbnd[jj];
			for (int ii = 0; ii < nx; ii++)
			{

				hh_g[ii + jj*nx] = hh[ii + jj*nx];
				hu_g[ii + jj*nx] = hh[ii + jj*nx];
				hv_g[ii + jj*nx] = hh[ii + jj*nx];
				hum_g[ii + jj*nx] = hh[ii + jj*nx];
				hvm_g[ii + jj*nx] = hh[ii + jj*nx];
				zb_g[ii + jj*nx] = zb[ii + jj*nx];
				zs_g[ii + jj*nx] = zs[ii + jj*nx];
				uu_g[ii + jj*nx] = uu[ii + jj*nx];
				vv_g[ii + jj*nx] = vv[ii + jj*nx];
				//CUDA_CHECK( cudaMemcpy(zom_g, zom, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );
				vu_g[ii + jj*nx] = uu[ii + jj*nx];
				uv_g[ii + jj*nx] = uu[ii + jj*nx];
				vev_g[ii + jj*nx] = uu[ii + jj*nx];
				ueu_g[ii + jj*nx] = uu[ii + jj*nx];
				vmageu_g[ii + jj*nx] = uu[ii + jj*nx];
				vmagev_g[ii + jj*nx] = uu[ii + jj*nx];
				ududx_g[ii + jj*nx] = uu[ii + jj*nx];
				vdudy_g[ii + jj*nx] = uu[ii + jj*nx];
				vdvdy_g[ii + jj*nx] = uu[ii + jj*nx];
				udvdx_g[ii + jj*nx] = uu[ii + jj*nx];


				H_g[ii + jj*nx] = uu[ii + jj*nx];
				Fx_g[ii + jj*nx] = uu[ii + jj*nx];
				Fy_g[ii + jj*nx] = uu[ii + jj*nx];
				urms_g[ii + jj*nx] = uu[ii + jj*nx];
				ust_g[ii + jj*nx] = uu[ii + jj*nx];
				D_g[ii + jj*nx] = uu[ii + jj*nx];


				for (int nt = 0; nt < ntheta; nt++)
				{
					ee_g[ii + jj*nx + nt*nx*ny] = ee[ii + jj*nx + nt*nx*ny];
					rr_g[ii + jj*nx + nt*nx*ny] = rr[ii + jj*nx + nt*nx*ny];
					cxsth_g[nt] = cxsth[nt];
					sxnth_g[nt] = sxnth[nt];
					theta_g[nt] = theta[nt];
				}
				thetamean_g[ii + jj*nx] = uu[ii + jj*nx];
				R_g[ii + jj*nx] = uu[ii + jj*nx];
				DR_g[ii + jj*nx] = uu[ii + jj*nx];
				c_g[ii + jj*nx] = uu[ii + jj*nx];
				cg_g[ii + jj*nx] = uu[ii + jj*nx];



				stdep_g[ii + jj*nx] = stdep[ii + jj*nx];
				Cc_g[ii + jj*nx] = uu[ii + jj*nx];
				kturb_g[ii + jj*nx] = uu[ii + jj*nx];
				rolthick_g[ii + jj*nx] = uu[ii + jj*nx];
				dzb_g[ii + jj*nx] = dzb[ii + jj*nx];



				//Averaged variables
				Hmean_g[ii + jj*nx] = uu[ii + jj*nx];
				uumean_g[ii + jj*nx] = uu[ii + jj*nx];
				vvmean_g[ii + jj*nx] = uu[ii + jj*nx];
				hhmean_g[ii + jj*nx] = uu[ii + jj*nx];
				zsmean_g[ii + jj*nx] = uu[ii + jj*nx];
				Cmean_g[ii + jj*nx] = uu[ii + jj*nx];
			}
		}
	}


	if (GPUDEVICE >= 0)
	{
		if (imodel == 1 || imodel > 2)
		{

			CUDA_CHECK(cudaMemcpy(qbndold_g, qbndold, 4 * ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
			CUDA_CHECK(cudaMemcpy(qbndnew_g, qbndnew, 4 * ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		}
		/*
		for (int ix=0; ix<nx; ix++)
		{
		for (int iy=0; iy<ny; iy++)
		{
		sigm[ix+iy*nx]=2*pi/Trep;
		for (int itheta=0; itheta<ntheta; itheta++)
		{
		sigt[ix+iy*nx+itheta*nx*ny] = sigm[ix+iy*nx];
		thet[ix+iy*nx+itheta*nx*ny] = theta[itheta];
		costhet[ix+iy*nx+itheta*nx*ny]= cos(theta[itheta]);
		sinthet[ix+iy*nx+itheta*nx*ny]= sin(theta[itheta]);
		}
		}
		}
		*/

		//Tsout= fopen (tsoutfile,"w");
		//fprintf(Tsout,"Totaltime\t hh\t zs\t uu\t vv\t H\n");


		
		dim3 blockDimLine(32, 1, 1);
		dim3 gridDimLine(ceil((nx*ny*1.0f) / blockDimLine.x), 1, 1);
		dim3 blockDim(16, 16, 1);
		dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);
		

		printf("gridDim=%i,%i,%i\n", gridDim.x, gridDim.y, gridDim.z);
		printf("gridDim=%i,%i,%i\n", gridDimLine.x, gridDimLine.y, gridDimLine.z);

		//Calculate bottomm friction based on initial hard layer file
		updatezom << <gridDim, blockDim, 0 >> >(nx, ny, cf, cf2, fw, fw2, stdep_g, cfm_g, fwm_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");`
		CUDA_CHECK(cudaThreadSynchronize());


		// Calculate initial maximum timestep


		FLOWDT << <gridDim, blockDim, 0 >> >(nx, ny, dx, 0.25f, dtflow_g, hh_g);
		CUDA_CHECK(cudaThreadSynchronize());


		minmaxKernel << <gridDimLine, blockDimLine, 0 >> >(nx*ny, arrmax_g, arrmin_g, dtflow_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		finalminmaxKernel << <1, blockDimLine, 0 >> >(arrmax_g, arrmin_g);
		CUDA_CHECK(cudaThreadSynchronize());

		//CUDA_CHECK(cudaMemcpy(arrmax, arrmax_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(arrmin, arrmin_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));

		//CUDA_CHECK(cudaMemcpy(hh, hh_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));

		//float hhmin=hh[0];

		//for (int ix = 0; ix < nx; ix++)
		//{
		//	for (int iy = 0; iy < ny; iy++)
		//	{
		//		hhmin = min(hhmin, hh[ix + iy*nx]);
		//	}
		//}



		dt = arrmin[0];

		printf("Initial timestep: dt=%f\n", arrmin[0]);


		
		if (imodel == 1 || imodel > 2)
		{
			set_bnd << <gridDim, blockDim, 0 >> >(nx, ny, Trep, ntheta, theta_g, sigm_g);

			//CUT_CHECK_ERROR("set_bnd() execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

			//also run dispersion relation, cg is needed in the first iteration of the bnd flow
			dispersion_init << <gridDim, blockDim, 0 >> >(nx, ny, twopi, g, aphi, bphi, sigm_g, hh_g, cg_g);
			//CUT_CHECK_ERROR("dispersion execution failed\n");
			CUDA_CHECK(cudaThreadSynchronize());

		}
		
	}
	else
	{


		//Calculate bottomm friction based on initial hard layer file
		updatezomCPU(nx, ny, cf, cf2, fw, fw2, stdep_g, cfm_g, fwm_g);


		if (imodel == 1 || imodel > 2)
		{
			set_bndCPU(nx, ny, Trep, ntheta, theta_g, sigm_g);



			//also run dispersion relation, cg is needed in the first iteration of the bnd flow
			dispersion_initCPU(nx, ny, twopi, g, aphi, bphi, sigm_g, hh_g, cg_g);


		}
	}
	// prepare output file
	printf("prepare output");
	creatncfile(tsoutfile, nx, ny, dx, 0.0f, imodel, zb, zs, uu, vv, H, H, thetamean, uu, uu, uu, uu, uu, uu, uu, hh, uu, uu, uu, uu, uu, uu);

	//create3dnc(nx,ny,ntheta,dx,0.0f,theta,ctheta);

	istepout = istepout + nstepout;
	printf("...done\n");


	printf("Starting Computation \n");

	// for hot start purposes
	//read2Dnc(nx,ny,"uufile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(uu_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"vvfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(vv_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"zsfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(zs_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read2Dnc(nx,ny,"hhfile.nc",uu);
	//CUDA_CHECK( cudaMemcpy(hh_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );


	// Run the model
	if (GPUDEVICE >= 0)
	{
		mainloopGPU();
	}
	else
	{
		mainloopCPU();
	}




	//close the bnd files and clean up a bit
	fclose(fsl);
	fclose(fwind);

	endcputime = clock();

	printf("Total runtime= %d  seconds\n", (endcputime - startcputime) / CLOCKS_PER_SEC);

	cudaThreadExit();

	//CUT_EXIT(argc, argv);










}

