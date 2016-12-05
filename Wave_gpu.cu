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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "XBeachGPU.h"
using DECNUM = float;

//global variables

DECNUM Trep, Trepold, Trepnew;
DECNUM * St, *Stnew, *Stold;
double * Stfile;
double * qfile;
double * Tpfile;
int nwbndstep = 0;
//int wavebndtype;
int nstep = 0;
//int breakmod = 1;

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
//DECNUM  nuh, nuhfac;
DECNUM *nuh_g;
DECNUM * viscu_g, *viscv_g;

DECNUM uumin = 0.0f;


//int nx, ny;
//DECNUM  dt, eps;
DECNUM *arrmin, *arrmax;
DECNUM *arrmin_g, *arrmax_g;
//DECNUM dx,grdalpha;
double totaltime;
int nwstp;//nb of hd step between wave step and next step for calculating waves and
DECNUM wdt;// wave model time step
DECNUM wavbndtime;
DECNUM slbndtime;
DECNUM windtime;

int SLstepinbnd, WNDstepinbnd;
//DECNUM Cd; //Wind drag
DECNUM fp, hm0gew, mainang, rt, scoeff, gam;
int nwavbnd, nwavfile;
DECNUM dtwavbnd;
int roller;
//DECNUM wci;
DECNUM *wci_g;
//DECNUM gammax, hwci, gammaa, n, alpha, beta, t1;
//DECNUM fw, fw2;
DECNUM *fwm_g;//Wave dissipation factor map

DECNUM phi = (1.0f + sqrt(5.0f)) / 2;
DECNUM aphi = 1 / (phi + 1);
DECNUM bphi = phi / (phi + 1);
DECNUM twopi = 8 * atan(1.0f);

DECNUM g = 9.81f;
DECNUM rho = 1025.0f;
DECNUM zo;
//DECNUM cf, cf2;//friction
DECNUM *cfm, *cfm_g; //friction map

//DECNUM lat; //lattitude 
//DECNUM fc; //coriolis

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

//DECNUM D50, D90, rhosed;

DECNUM * kturb_g, *rolthick_g, *dzsdt_g;
DECNUM * Ceq_g, *ceqsg_g, *ceqbg_g, *Tsg_g, *facero_g;
DECNUM * C, *Cc_g, *stdep, *stdep_g;
DECNUM * Sus_g, *Svs_g, *Sub_g, *Svb_g;
//DECNUM morfac, por;
DECNUM * ero_g, *depo_g;
DECNUM *dzb, *dzb_g, *ddzb_g;
DECNUM * Sout_g;
int * indSub_g, *indSvb_g;
//DECNUM sus = 1.0f;
//DECNUM bed = 1.0f;
//DECNUM facsk, facas;

DECNUM * zsbnd;
DECNUM rtsl;
DECNUM zsbndnew, zsbndold;


DECNUM windu, windv;

FILE * fwav;


FILE * Tsout;
//FILE * fXq,* fXE;

//char tsoutfile[256];
int iout, jout;
//int imodel; //1: Wave only; 2: Current Only; 3: Wave + current; 4: Wave + current + Sediment transport



int nstepout = 0;//nb of step between outputs. 0= no output
int istepout = 0;//

//double outputtimestep;
double nextoutputtime;



int nstepplot = 0;//nb of step between plots. 0= no plotting
int istepplot = 1;
int displayon = 0;
//double endtime;



//DECNUM wws;
//DECNUM drydzmax;
//DECNUM wetdzmax;
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




int startflowstep;
//int usesmago;

// Particle stuff
int npart;
//double sedstart; // Warm up time in seconds before starting sediment transport simulation 
int wxstep = 1;


//char wavebndfile[256];




//#include "read_input.cpp"



#include "Wave_kernel.cu"
#include "Flow_kernel.cu"
#include "Sediment_kernel.cu"
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
void mainloopGPU(XBGPUParam Param, std::vector<SLBnd> slbnd, std::vector<WindBnd> wndbnd)
{
	double dt = Param.dt;
	int nx, ny;
	nx = Param.nx;
	ny = Param.ny;
	
	
	//Warning hard wired stuff here
	
	std::vector<Pointout> zsout;
	
	std::vector< std::vector< Pointout > > zsAllout;
	
	Pointout stepread;
	




	FILE * fsSLTS;
	if (!Param.TSoutfile.empty())
	{
		for (int o=0; o < Param.TSoutfile.size(); o++)
		{
			//Overwrite existing files
			fsSLTS = fopen(Param.TSoutfile[o].c_str(), "w");
			fclose(fsSLTS);

			// Add empty row for each output point
			zsAllout.push_back(std::vector<Pointout>());
		}
	}
	
	//Map a link between the variable s a string and the variable pointer
	// This may need to be stored somewhere else??
	std::map<std::string, DECNUM *> OutputVarMapCPU;
	std::map<std::string, DECNUM *> OutputVarMapGPU;
	std::map<std::string, int> OutputVarMaplen;

	OutputVarMapCPU["zb"] = zb;
	OutputVarMapGPU["zb"] = zb_g;
	if (Param.morphology == 1)
	{
		OutputVarMaplen["zb"] = nx*ny;
	}
	else
	{
		OutputVarMaplen["zb"] = 0;
	}
	OutputVarMapCPU["uu"] = uu;
	OutputVarMapGPU["uu"] = uu_g;
	OutputVarMaplen["uu"] = nx*ny;

	OutputVarMapCPU["vv"] = vv;
	OutputVarMapGPU["vv"] = vv_g;
	OutputVarMaplen["vv"] = nx*ny;

	OutputVarMapCPU["zs"] = zs;
	OutputVarMapGPU["zs"] = zs_g;
	OutputVarMaplen["zs"] = nx*ny;

	OutputVarMapCPU["hh"] = hh;
	OutputVarMapGPU["hh"] = hh_g;
	OutputVarMaplen["hh"] = nx*ny;

	OutputVarMapCPU["H"] = H;
	OutputVarMapGPU["H"] = H_g;
	OutputVarMaplen["H"] = nx*ny;

	OutputVarMapCPU["thetamean"] = thetamean;
	OutputVarMapGPU["thetamean"] = thetamean_g;
	OutputVarMaplen["thetamean"] = nx*ny;

	OutputVarMapCPU["ee"] = ee;
	OutputVarMapGPU["ee"] = ee_g;
	OutputVarMaplen["ee"] = nx*ny*ntheta;

	OutputVarMapCPU["rr"] = rr;
	OutputVarMapGPU["rr"] = rr_g;
	OutputVarMaplen["rr"] = nx*ny*ntheta;


	OutputVarMapCPU["cfm"] = cfm;
	OutputVarMapGPU["cfm"] = cfm_g;
	OutputVarMaplen["cfm"] = nx*ny;

	OutputVarMapCPU["dzb"] = dzb;
	OutputVarMapGPU["dzb"] = dzb_g;
	OutputVarMaplen["dzb"] = nx*ny;

	OutputVarMapCPU["stdep"] = stdep;
	OutputVarMapGPU["stdep"] = stdep_g;
	OutputVarMaplen["stdep"] = nx*ny;

	OutputVarMapCPU["Fx"] = Fx;
	OutputVarMapGPU["Fx"] = Fx_g;
	OutputVarMaplen["Fx"] = nx*ny;

	OutputVarMapCPU["Fy"] = Fy;
	OutputVarMapGPU["Fy"] = Fy_g;
	OutputVarMaplen["Fy"] = nx*ny;

	OutputVarMapCPU["cgx"] = cgx;
	OutputVarMapGPU["cgx"] = cgx_g;
	OutputVarMaplen["cgx"] = nx*ny*ntheta;

	OutputVarMapCPU["cgy"] = cgy;
	OutputVarMapGPU["cgy"] = cgy_g;
	OutputVarMaplen["cgy"] = nx*ny*ntheta;

	OutputVarMapCPU["cx"] = cx;
	OutputVarMapGPU["cx"] = cx_g;
	OutputVarMaplen["cx"] = nx*ny*ntheta;

	OutputVarMapCPU["cy"] = cy;
	OutputVarMapGPU["cy"] = cy_g;
	OutputVarMaplen["cy"] = nx*ny*ntheta;

	OutputVarMapCPU["ctheta"] = ctheta;
	OutputVarMapGPU["ctheta"] = ctheta_g;
	OutputVarMaplen["ctheta"] = nx*ny*ntheta;


	OutputVarMapCPU["D"] = D;
	OutputVarMapGPU["D"] = D_g;
	OutputVarMaplen["D"] = nx*ny;

	OutputVarMapCPU["E"] = E;
	OutputVarMapGPU["E"] = E_g;
	OutputVarMaplen["E"] = nx*ny;


	OutputVarMapCPU["urms"] = urms;
	OutputVarMapGPU["urms"] = urms_g;
	OutputVarMaplen["urms"] = nx*ny;

	OutputVarMapCPU["ueu"] = ueu;
	OutputVarMapGPU["ueu"] = ueu_g;
	OutputVarMaplen["ueu"] = nx*ny;

	OutputVarMapCPU["vev"] = vev;
	OutputVarMapGPU["vev"] = vev_g;
	OutputVarMaplen["vev"] = nx*ny;


	OutputVarMapCPU["hhmean"] = hhmean;
	OutputVarMapGPU["hhmean"] = hhmean_g;
	OutputVarMaplen["hhmean"] = nx*ny;

	OutputVarMapCPU["uumean"] = uumean;
	OutputVarMapGPU["uumean"] = uumean_g;
	OutputVarMaplen["uumean"] = nx*ny;

	OutputVarMapCPU["vvmean"] = vvmean;
	OutputVarMapGPU["vvmean"] = vvmean_g;
	OutputVarMaplen["vvmean"] = nx*ny;

	OutputVarMapCPU["zsmean"] = zsmean;
	OutputVarMapGPU["zsmean"] = zsmean_g;
	OutputVarMaplen["zsmean"] = nx*ny;

	OutputVarMapCPU["Hmean"] = Hmean;
	OutputVarMapGPU["Hmean"] = Hmean_g;
	OutputVarMaplen["Hmean"] = nx*ny;

	OutputVarMapCPU["Cmean"] = Cmean;
	OutputVarMapGPU["Cmean"] = Cmean_g;
	OutputVarMaplen["Cmean"] = nx*ny;

	OutputVarMapCPU["sigm"] = sigm;
	OutputVarMapGPU["sigm"] = sigm_g;
	OutputVarMaplen["sigm"] = nx*ny;

	OutputVarMapCPU["k"] = k;
	OutputVarMapGPU["k"] = k_g;
	OutputVarMaplen["k"] = nx*ny;

	OutputVarMapCPU["c"] = c;
	OutputVarMapGPU["c"] = c_g;
	OutputVarMaplen["c"] = nx*ny;

	OutputVarMapCPU["kh"] = kh;
	OutputVarMapGPU["kh"] = kh_g;
	OutputVarMaplen["kh"] = nx*ny;

	OutputVarMapCPU["cg"] = cg;
	OutputVarMapGPU["cg"] = cg_g;
	OutputVarMaplen["cg"] = nx*ny;

	OutputVarMapCPU["sinh2kh"] = sinh2kh;
	OutputVarMapGPU["sinh2kh"] = sinh2kh_g;
	OutputVarMaplen["sinh2kh"] = nx*ny;

	OutputVarMapCPU["dhdx"] = dhdx;
	OutputVarMapGPU["dhdx"] = dhdx_g;
	OutputVarMaplen["dhdx"] = nx*ny;

	OutputVarMapCPU["dhdy"] = dhdy;
	OutputVarMapGPU["dhdy"] = dhdy_g;
	OutputVarMaplen["dhdy"] = nx*ny;

	OutputVarMapCPU["dudx"] = dudx;
	OutputVarMapGPU["dudx"] = dudx_g;
	OutputVarMaplen["dudx"] = nx*ny;

	OutputVarMapCPU["dudy"] = dudy;
	OutputVarMapGPU["dudy"] = dudy_g;
	OutputVarMaplen["dudy"] = nx*ny;

	OutputVarMapCPU["dvdx"] = dvdx;
	OutputVarMapGPU["dvdx"] = dvdx_g;
	OutputVarMaplen["dvdx"] = nx*ny;

	OutputVarMapCPU["dvdy"] = dvdy;
	OutputVarMapGPU["dvdy"] = dvdy_g;
	OutputVarMaplen["dvdy"] = nx*ny;

	OutputVarMapCPU["C"] = C;
	OutputVarMapGPU["C"] = Cc_g;
	OutputVarMaplen["C"] = nx*ny;

	OutputVarMapCPU["R"] = R;
	OutputVarMapGPU["R"] = R_g;
	OutputVarMaplen["R"] = nx*ny;

	OutputVarMapCPU["DR"] = DR;
	OutputVarMapGPU["DR"] = DR_g;
	OutputVarMaplen["DR"] = nx*ny;

	///////////////////////////////////////////////////////
	// Warning the variable bellow were never allocted on the CPU
	// This could lead to errors if asynchroneous meme copy is to be implemented for output

	OutputVarMapCPU["wci"] = dummy;
	OutputVarMapGPU["wci"] = wci_g;
	OutputVarMaplen["wci"] = nx*ny;

	OutputVarMapCPU["vmageu"] = dummy;
	OutputVarMapGPU["vmageu"] = vmageu_g;
	OutputVarMaplen["vmageu"] = nx*ny;

	OutputVarMapCPU["vmagev"] = dummy;
	OutputVarMapGPU["vmagev"] = vmagev_g;
	OutputVarMaplen["vmagev"] = nx*ny;

	OutputVarMapCPU["dzsdx"] = dummy;
	OutputVarMapGPU["dzsdx"] = dzsdx_g;
	OutputVarMaplen["dzsdx"] = nx*ny;

	OutputVarMapCPU["dzsdy"] = dummy;
	OutputVarMapGPU["dzsdy"] = dzsdy_g;
	OutputVarMaplen["dzsdy"] = nx*ny;
	
	OutputVarMapCPU["dzsdt"] = dummy;
	OutputVarMapGPU["dzsdt"] = dzsdt_g;
	OutputVarMaplen["dzsdt"] = nx*ny;

	OutputVarMapCPU["fwm"] = dummy;
	OutputVarMapGPU["fwm"] = fwm_g;
	OutputVarMaplen["fwm"] = nx*ny;

	OutputVarMapCPU["hu"] = dummy;
	OutputVarMapGPU["hu"] = hu_g;
	OutputVarMaplen["hu"] = nx*ny;

	OutputVarMapCPU["hv"] = dummy;
	OutputVarMapGPU["hv"] = hv_g;
	OutputVarMaplen["hv"] = nx*ny;

	OutputVarMapCPU["hum"] = dummy;
	OutputVarMapGPU["hum"] = hum_g;
	OutputVarMaplen["hum"] = nx*ny;

	OutputVarMapCPU["hvm"] = dummy;
	OutputVarMapGPU["hvm"] = hvm_g;
	OutputVarMaplen["hvm"] = nx*ny;

	OutputVarMapCPU["uv"] = dummy;
	OutputVarMapGPU["uv"] = uv_g;
	OutputVarMaplen["uv"] = nx*ny;

	OutputVarMapCPU["vu"] = dummy;
	OutputVarMapGPU["vu"] = vu_g;
	OutputVarMaplen["vu"] = nx*ny;

	//OutputVarMapCPU["wetu"] = dummy;
	//OutputVarMapGPU["wetu"] = wetu_g;
	//OutputVarMaplen["wetu"] = nx*ny;

	//OutputVarMapCPU["wetv"] = dummy;
	//OutputVarMapGPU["wetv"] = wetv_g;
	//OutputVarMaplen["wetv"] = nx*ny;

	OutputVarMapCPU["ududx"] = dummy;
	OutputVarMapGPU["ududx"] = ududx_g;
	OutputVarMaplen["ududx"] = nx*ny;

	OutputVarMapCPU["vdvdy"] = dummy;
	OutputVarMapGPU["vdvdy"] = vdvdy_g;
	OutputVarMaplen["vdvdy"] = nx*ny;

	OutputVarMapCPU["udvdx"] = dummy;
	OutputVarMapGPU["udvdx"] = udvdx_g;
	OutputVarMaplen["udvdx"] = nx*ny;

	OutputVarMapCPU["vdudy"] = dummy;
	OutputVarMapGPU["vdudy"] = ududx_g;
	OutputVarMaplen["vdudy"] = nx*ny;

	OutputVarMapCPU["ust"] = dummy;
	OutputVarMapGPU["ust"] = ust_g;
	OutputVarMaplen["ust"] = nx*ny;

	OutputVarMapCPU["kturb"] = dummy;
	OutputVarMapGPU["kturb"] = kturb_g;
	OutputVarMaplen["kturb"] = nx*ny;

	OutputVarMapCPU["rolthick"] = dummy;
	OutputVarMapGPU["rolthick"] = rolthick_g;
	OutputVarMaplen["rolthick"] = nx*ny;

	OutputVarMapCPU["ceqsg"] = dummy;
	OutputVarMapGPU["ceqsg"] = ceqsg_g;
	OutputVarMaplen["ceqsg"] = nx*ny;

	//< or <= ? crashes with <= if the boundary limit is == to endtime
	while (totaltime < Param.endtime)
	{

		//std::cout << totaltime << " "<< Param.endtime << std::endl;
		
		dim3 blockDim(16, 16, 1);// This means that the grid has to be a factor of 16 on both x and y
		dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);

		dim3 blockDimLine(32, 1, 1);
		dim3 gridDimLine(ceil((nx*ny*1.0f) / blockDimLine.x), 1, 1);

		nstep++;
		wdt = dt; // previous timestep

		//Calculate timestep
		FLOWDT << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.CFL, dtflow_g, hh_g);
		CUDA_CHECK(cudaDeviceSynchronize());


		minmaxKernel << <gridDimLine, blockDimLine, 0 >> >(nx*ny, arrmax_g, arrmin_g, dtflow_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		finalminmaxKernel << <1, blockDimLine, 0 >> >(arrmax_g, arrmin_g);
		CUDA_CHECK(cudaDeviceSynchronize());

		//CUDA_CHECK(cudaMemcpy(arrmax, arrmax_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(arrmin, arrmin_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));

		dt = arrmin[0] * 0.5; // Not sure why this is but it is in the original XBeach!!

		
		
		if ((Param.swave == 1 ) && totaltime>0.0)
		{
			double dtwave;
			// Make sure the CFL condition for flow do not violate CFL condition for Waves
			WAVEDT << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.CFL, dtheta, dtflow_g, ctheta_g);
			CUDA_CHECK(cudaDeviceSynchronize());


			minmaxKernel << <gridDimLine, blockDimLine, 0 >> >(nx*ny, arrmax_g, arrmin_g, dtflow_g);
			//CUT_CHECK_ERROR("UpdateZom execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			finalminmaxKernel << <1, blockDimLine, 0 >> >(arrmax_g, arrmin_g);
			CUDA_CHECK(cudaDeviceSynchronize());

			//CUDA_CHECK(cudaMemcpy(arrmax, arrmax_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(arrmin, arrmin_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));

			dtwave = arrmin[0]*0.5;

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

		double nextstop = min(nextoutputtime, Param.endtime);

		//need to limit timestep so that it matches the exact output time or endtime
		if (ceil((nextstop - totaltime) / dt)> 0.0)
		{
			dt = (nextstop - totaltime) / ceil((nextstop - totaltime) / dt);
		}


		//printf("Timestep: %f\n", dt);
		Param.dt = dt;

		totaltime = totaltime + dt;	//total run time acheived until now in s
		//std::cout << totaltime << std::endl;
		

		if (Param.swave == 1 )
		{
			wavebnd(Param); // Calculate the boundary condition for this step
		}

		if (Param.flow == 1)
		{
			flowbnd(Param, slbnd, wndbnd);// Calculate the flow boundary for this step
		}

		if (Param.swave == 1 )
		{

			wavestep(Param); // Calculate the wave action ballance for this step
		}



		if (Param.flow == 1)
		{
			flowstep(Param);// solve the shallow water and continuity for this step
		}
		if (Param.morphology == 1 && totaltime >= Param.sedstart)
		{
			//Sediment step
			sedimentstep(Param);//solve the sediment dispersion, and morphology
		}

		//add last value for avg calc
		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Hmean_g, H_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, uumean_g, uu_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, vvmean_g, vv_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, hhmean_g, hh_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, zsmean_g, zs_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		addavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Cmean_g, Cc_g);
		//CUT_CHECK_ERROR("Add avg execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());


		//////////////////////////////////////////
		//BIG
		//WARNING HERE -- NEED TO MAKE ASYNC
		/////////////////////////////////////////
		if (!Param.TSnodesout.empty())
		{
			CUDA_CHECK(cudaMemcpy(zs, zs_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(H, H_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
		}
		//after the transfert


		if (!Param.TSnodesout.empty())
		{
			for (int o=0; o < Param.TSnodesout.size(); o++)
			{
				stepread.time = totaltime;
				stepread.zs = zs[Param.TSnodesout[o].i + Param.TSnodesout[o].j*nx];
				stepread.H = H[Param.TSnodesout[o].i + Param.TSnodesout[o].j*nx];
				zsAllout[o].push_back(stepread);
			}
		}


		if (nextoutputtime-totaltime <= dt*0.001  && Param.outputtimestep > 0)
		{
			istepout = istepout + nstepout;

			//
			nextoutputtime = nextoutputtime + Param.outputtimestep;

			//Avg mean variables

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstep, Hmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstep, uumean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstep, vvmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstep, hhmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstep, zsmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			divavg_var << <gridDim, blockDim, 0 >> >(nx, ny, nstep, Cmean_g);
			//CUT_CHECK_ERROR("Div avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			//For each requested variables

			if (!Param.outvars.empty())
			{
				writenctimestep(Param, totaltime);
				
				for (int ivar = 0; ivar < Param.outvars.size(); ivar++)
				{
					if (OutputVarMaplen[Param.outvars[ivar]] > 0)
					{
						//Should be async
						CUDA_CHECK(cudaMemcpy(OutputVarMapCPU[Param.outvars[ivar]], OutputVarMapGPU[Param.outvars[ivar]], OutputVarMaplen[Param.outvars[ivar]] * sizeof(DECNUM), cudaMemcpyDeviceToHost));
						//Create definition for each variable and store it
						writencvarstep(Param, Param.outvars[ivar], OutputVarMapCPU[Param.outvars[ivar]]);
					}
				}
			}



			// Download mean vars
			//CUDA_CHECK(cudaMemcpy(Hmean, Hmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(uumean, uumean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(vvmean, vvmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(hhmean, hhmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(zsmean, zsmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(Cmean, Cmean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));


			//CUDA_CHECK(cudaMemcpy(H, H_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(uu, uu_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(vv, vv_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(zs, zs_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(Fx, Fx_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(Fy, Fy_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(thetamean, thetamean_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(D, D_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(urms, urms_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(ueu, ueu_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK(cudaMemcpy(vev, vev_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK( cudaMemcpy(C, ceqsg_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			//CUDA_CHECK(cudaMemcpy(C, hum_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//CUDA_CHECK( cudaMemcpy(C,k_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			//CUDA_CHECK( cudaMemcpy(ctheta,ee_g, nx*ny*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			//CUDA_CHECK(cudaMemcpy(hh, hh_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//if (Param.morphology == 1 )// If moprhology is on
			//{
			//	CUDA_CHECK(cudaMemcpy(zb, zb_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//	CUDA_CHECK(cudaMemcpy(dzb, dzb_g, nx*ny*sizeof(DECNUM), cudaMemcpyDeviceToHost));
			//}
			//CUDA_CHECK( cudaMemcpy(xxp, xxp_g, npart*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			//CUDA_CHECK( cudaMemcpy(yyp, yyp_g, npart*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
			printf("Writing output, totaltime:%f s, Mean dt=%f\n", totaltime, Param.outputtimestep/nstep);
			write_text_to_log_file("Writing outputs, totaltime: " + std::to_string(totaltime) + ", Mean dt= " + std::to_string(Param.outputtimestep / nstep));
			//writestep2nc(tsoutfile, nx, ny,/*npart,*/(float) totaltime, imodel,/*xxp,yyp,*/zb, zs, uu, vv, H, H, thetamean, D, urms, ueu, vev, C, dzb, Fx, Fy, hh, Hmean, uumean, vvmean, hhmean, zsmean, Cmean);
			//writestep2nc(Param,(float)totaltime, zb, zs, uu, vv, H, H, thetamean, D, urms, ueu, vev, C, dzb, Fx, Fy, hh, Hmean, uumean, vvmean, hhmean, zsmean, Cmean);
			
			// Write a sigle timeseries in the file
			if (!Param.TSoutfile.empty())
			{
				for (int o=0; o < Param.TSoutfile.size(); o++)
				{
					//Overwrite existing files
					fsSLTS = fopen(Param.TSoutfile[o].c_str(), "a");
					for (int n = 0; n < zsAllout[o].size(); n++)
					{
						fprintf(fsSLTS, "%f\t%.4f\t%.4f\n", zsAllout[o][n].time, zsAllout[o][n].zs, zsAllout[o][n].H);
					}
					fclose(fsSLTS);
					//reset zsout
					zsAllout[o].clear();
					//zsAllout.push_back(std::vector<SLBnd>());
				}
			}
			
			
			
			
			//write3dvarnc(nx,ny,ntheta,totaltime,ctheta);
			//outfile[],nx,ny,npart,totaltime,xxp,yyp,zs,uu, vv, H,Tp,Dp,      D,Urms,ueu,vev)
			//fprintf(Tsout,"%f\t%f\t%f\t%f\t%f\t%f\n",totaltime,hh[iout+jout*nx],zs[iout+jout*nx],uu[iout+jout*nx],vv[iout+jout*nx],H[iout+jout*nx]);


			//Clear avg vars
			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Hmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, uumean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, vvmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, hhmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, zsmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			resetavg_var << <gridDim, blockDim, 0 >> >(nx, ny, Cmean_g);
			//CUT_CHECK_ERROR("Reset avg execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			nstep = 0;
		}
	}
	

	
	
}


void mainloopCPU(XBGPUParam Param,std::vector<SLBnd> slbnd, std::vector<WindBnd> wndbnd)
{
	printf("Computing CPU mode\n");

	int nx, ny;
	nx = Param.nx;
	ny = Param.ny;

	while (totaltime <= Param.endtime)
	{

		nstep++;
		wdt = Param.dt; // Sometinmes in stationary wave run one can have a larger wdt (wave time step)
		totaltime = totaltime + Param.dt;	//total run time acheived until now in s



		if (Param.swave == 1 )
		{
			wavebnd(Param); // Calculate the boundary condition for this step
		}

		if (Param.flow == 1)
		{
			flowbnd(Param,slbnd,wndbnd);// Calculate the flow boundary for this step
		}
		if (Param.swave == 1)
		{
			wavestepCPU(Param); // Calculate the wave action ballance for this step
		}
		if (Param.flow == 1)
		{
			//flowstepCPU();// solve the shallow water and continuity for this step
		}
		if (Param.morphology == 1 && totaltime >= Param.sedstart)
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
			writestep2nc(Param, (float)totaltime, zb_g, zs_g, uu_g, vv_g, H_g, xadvec_g, thetamean_g, D_g, urms_g, ueu_g, vev_g, Cc_g, dzb_g, Fx_g, Fy_g, hh_g, Hmean_g, uumean_g, vvmean_g, hhmean_g, zsmean_g, Cmean_g);
			
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




void flowbnd(XBGPUParam Param,std::vector<SLBnd> slbnd, std::vector<WindBnd> wndbnd)
{
	
	double zsbndi;
	int stepinbnd;
	int nx, ny;
	
	nx = Param.nx;
	ny = Param.ny;
	//update sl bnd
	
	// find next timestep
	
	double difft = slbnd[SLstepinbnd].time - totaltime;
	
	if (difft < 0.0)
	{
		SLstepinbnd++;
	}

	zsbndi = interptime(slbnd[SLstepinbnd].wlev, slbnd[SLstepinbnd - 1].wlev, slbnd[SLstepinbnd].time - slbnd[SLstepinbnd - 1].time, totaltime - slbnd[SLstepinbnd - 1].time);


	//std::cout << "i= " << stepinbnd << "; " << zsbndi << "\n" << std::endl;



	//if (totaltime >= slbndtime)
	//{

	//	zsbndold = zsbndnew;
	//	rtsl = slbndtime;
	//	fscanf(fsl, "%f\t%f", &slbndtime, &zsbndnew);
		//slbndtime=+rtsl;
		//zsbnd=zsbndold+(t-slbndtime+rtsl)*(zsbndnew-zsbndold)/rtsl;
	//}





	if (Param.wavebndtype == 1)
	{
		for (int ni = 0; ni < ny; ni++)
		{
			zsbnd[ni] = zsbndi;//zsbndold + ((float) totaltime - rtsl)*(zsbndnew - zsbndold) / (slbndtime - rtsl);
		}
	}

	if (Param.wavebndtype == 2)
	{
		if (Param.GPUDEVICE >= 0)
		{
			dim3 blockDim(16, 16, 1);
			dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);
			// FLow abs_2d should be here not at the flow step		
			// Set weakly reflective offshore boundary
			ubnd1D << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.dt, Param.g, Param.rho, (float)totaltime, wavbndtime, dtwavbnd, zsbndi, Trep, qbndold_g, qbndnew_g, zs_g, uu_g, vv_g, vu_g, umeanbnd_g, vmeanbnd_g, zb_g, cg_g, hum_g, cfm_g, Fx_g, hh_g);
			//CUT_CHECK_ERROR("ubnd execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());
		}
		else
		{
			ubndCPU(nx, ny, Param.dx, Param.dt, Param.g, Param.rho, (float)totaltime, wavbndtime, dtwavbnd, zsbndi, Trep, qbndold_g, qbndnew_g, zs_g, uu_g, vv_g, vu_g, umeanbnd_g, vmeanbnd_g, zb_g, cg_g, hum_g, cfm_g, Fx_g, hh_g);

		}
	}

	
	difft = wndbnd[WNDstepinbnd].time - totaltime;
	if (difft < 0.0)
	{
		WNDstepinbnd++;
	}
	windu = interptime(wndbnd[WNDstepinbnd].U, wndbnd[WNDstepinbnd - 1].U, wndbnd[WNDstepinbnd].time - wndbnd[WNDstepinbnd - 1].time, totaltime - wndbnd[WNDstepinbnd - 1].time);
	windv = interptime(wndbnd[WNDstepinbnd].V, wndbnd[WNDstepinbnd - 1].V, wndbnd[WNDstepinbnd].time - wndbnd[WNDstepinbnd - 1].time, totaltime - wndbnd[WNDstepinbnd - 1].time);

}


void flowstep(XBGPUParam Param)
{
	int nx, ny;

	nx = Param.nx;
	ny = Param.ny;
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

	
	// Set weakly reflective offshore boundary ! MOVED TO FLOW BND SUBROUTINE!!
	//ubnd<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,g,rho,totaltime,wavbndtime,dtwavbnd,slbndtime,rtsl,zsbndold,zsbndnew,Trep,qbndold_g,qbndnew_g,zs_g,uu_g,vv_g,vu_g,umeanbnd_g,vmeanbnd_g,zb_g,cg_g,hum_g,cfm_g,Fx_g,hh_g);
	//CUT_CHECK_ERROR("ubnd execution failed\n");
	//CUDA_CHECK( cudaDeviceSynchronize() );


	//
	// Water level slopes
	//
	wlevslopes << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.eps, zs_g, dzsdx_g, dzsdy_g, hh_g);
	//CUT_CHECK_ERROR("wlevslopes execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());



	//
	// Water depth at u pts for momentum and continuity eq (hum hu)
	//

	udepthmomcont << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.eps, uumin, wetu_g, zs_g, uu_g, hh_g, hum_g, hu_g, zb_g);
	//CUT_CHECK_ERROR("udepthmomcont execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	//
	// Water depth at v pts for momentum and continuity eq (hvm hv)
	//

	vdepthmomcont << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.eps, uumin, wetv_g, zs_g, vv_g, hh_g, hvm_g, hv_g, zb_g);
	//CUT_CHECK_ERROR("vdepthmomcont execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());






	//
	// Advection in the x direction using 2n order finite difference
	//

	ududx_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, hu_g, hum_g, uu_g, ududx_g);
	//CUT_CHECK_ERROR("uadvec execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//vdudy
	vdudy_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, hv_g, hum_g, uu_g, vv_g, vdudy_g);
	//CUT_CHECK_ERROR("uadvec execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());




	//
	// Smagorinsky formulation or Normal eddy viscosity
	//
	CUDA_CHECK(cudaMalloc((void **)&nuh_g, nx*ny*sizeof(DECNUM)));
	smago << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, uu_g, vv_g, Param.nuh, nuh_g, Param.usesmago);
	//CUT_CHECK_ERROR("uadvec execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	//
	// increase eddy viscosity by wave induced breaking as in Reniers 2004 & Set viscu = 0.0 near water line
	//
	CUDA_CHECK(cudaMalloc((void **)&viscu_g, nx*ny*sizeof(DECNUM)));
	viscou << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.rho, Param.eps, Param.nuhfac, nuh_g, hh_g, hum_g, hvm_g, DR_g, uu_g, wetu_g, viscu_g);
	//CUT_CHECK_ERROR("visco execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Explicit Euler step momentum u-direction
	//

	eulerustep << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.dt, Param.g, Param.rho, cfm_g, Param.fc, windu, Param.Cd, uu_g, urms_g, ududx_g, vdudy_g, viscu_g, dzsdx_g, hu_g, hum_g, Fx_g, vu_g, ueu_g, vmageu_g, wetu_g);
	//CUT_CHECK_ERROR("eulerustep execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd << <gridDim, blockDim, 0 >> >(nx, ny, uu_g, vv_g, zs_g);
	//CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Advection in the y direction using 2n order finite difference
	//
	//vdvdy
	vdvdy_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, hv_g, hvm_g, vv_g, vdvdy_g);
	//CUT_CHECK_ERROR("vadvec for v execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());
	//udvdx

	udvdx_adv2 << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, hu_g, hvm_g, uu_g, vv_g, udvdx_g);
	//CUT_CHECK_ERROR("vadvec for v execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	//
	// increase eddy viscosity by wave induced breaking as in Reniers 2004 & Set viscv = 0.0 near water line
	//
	CUDA_CHECK(cudaMalloc((void **)&viscv_g, nx*ny*sizeof(DECNUM)));
	viscov << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.rho, Param.eps, Param.nuhfac, nuh_g, hh_g, hum_g, hvm_g, DR_g, vv_g, wetv_g, viscv_g);
	//CUT_CHECK_ERROR("visco v execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());
	CUDA_CHECK(cudaFree(nuh_g));


	viscovbnd << <gridDim, blockDim, 0 >> >(nx, ny, viscv_g);
	//CUT_CHECK_ERROR("visco v execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Explicit Euler step momentum v-direction
	//

	eulervstep << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.dt, Param.g, Param.rho, cfm_g, Param.fc, windv, Param.Cd, vv_g, urms_g, udvdx_g, vdvdy_g, viscv_g, dzsdy_g, hv_g, hvm_g, Fy_g, uv_g, vev_g, vmagev_g, wetv_g);
	//CUT_CHECK_ERROR("eulervstep execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd << <gridDim, blockDim, 0 >> >(nx, ny, uu_g, vv_g, zs_g);
	//CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());



	//
	//v velocities at u pts and u velocities at v pts
	//

	calcuvvu << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, uu_g, vv_g, vu_g, uv_g, ust_g, thetamean_g, ueu_g, vev_g, vmageu_g, vmagev_g, wetu_g, wetv_g);
	//CUT_CHECK_ERROR("calcuvvu execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	uvlatbnd << <gridDim, blockDim, 0 >> >(nx, ny, vu_g, uv_g, ueu_g, vev_g, vmageu_g, vmagev_g);
	//fix side bnd for vu
	//twodimbndnoix<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,vu_g);
	//CUT_CHECK_ERROR("wave force X bnd execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());









	//
	//Calculate hu
	//
	depthhu << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, uumin, Param.eps, hh_g, uu_g, hu_g, zs_g, zb_g);
	//CUT_CHECK_ERROR("depthhu execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	//
	//Calculate hv
	//
	depthhv << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, uumin, Param.eps, hh_g, vv_g, hv_g, zs_g, zb_g);
	//CUT_CHECK_ERROR("depthhv execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Update water level using continuity eq.
	//
	continuity << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, Param.dt, Param.eps, uu_g, hu_g, vv_g, hv_g, zs_g, hh_g, zb_g, dzsdt_g);
	//CUT_CHECK_ERROR("continuity execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd << <gridDim, blockDim, 0 >> >(nx, ny, uu_g, vv_g, zs_g);
	//CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	hsbnd << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, hh_g, zb_g, zs_g);
	//CUT_CHECK_ERROR("hh lateral bnd execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());



	CUDA_CHECK(cudaFree(viscu_g));
	CUDA_CHECK(cudaFree(viscv_g));


}

void sedimentstep(XBGPUParam Param)
{
	int nx, ny;

	nx = Param.nx;
	ny = Param.ny;
	double dx = Param.dx;
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
	longturb << <gridDim, blockDim, 0 >> >(nx, ny, dx, Param.rho, Param.g, Param.dt, Param.beta, c_g, kturb_g, rolthick_g, dzsdt_g, uu_g, vv_g, hu_g, hv_g, wetu_g, wetv_g, hh_g);
	//CUT_CHECK_ERROR("longturb execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	//
	// Calculate Equilibrium concentration Ceq
	//
	//CUDA_CHECK( cudaMalloc((void **)&ceqsg_g, nx*ny*sizeof(DECNUM )) );
	CUDA_CHECK(cudaMalloc((void **)&ceqbg_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Tsg_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&ua_g, nx*ny*sizeof(DECNUM)));
	//BEWARE BELOW SHOULD BE hh_old_g
	//Sbvr or Sednew
	Sbvr << <gridDim, blockDim, 0 >> >(nx, ny, Param.rho, Param.g, Param.eps, Trep, Param.D50, Param.D90, Param.rhosed, Param.wws, Param.nuhfac, ueu_g, vev_g, H_g, DR_g, R_g, c_g, hh_g, urms_g, ceqsg_g, ceqbg_g, Tsg_g, cfm_g, kturb_g);
	//CUT_CHECK_ERROR("CalcCeq execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	Rvr << <gridDim, blockDim, 0 >> >(nx, ny, Trep, Param.facsk, Param.facas, H_g, hh_g, urms_g, c_g, ua_g);
	//CUT_CHECK_ERROR("Rvr execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());

	// CUDA_CHECK( cudaMemcpy(uumean, ua_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );




	//
	// Limit erosion to available sediment on top of hard layer
	//
	CUDA_CHECK(cudaMalloc((void **)&facero_g, nx*ny*sizeof(DECNUM)));
	Erosus << <gridDim, blockDim, 0 >> >(nx, ny, Param.dt, Param.morfac, Param.por, hh_g, ceqsg_g, ceqbg_g, Tsg_g, facero_g, stdep_g);
	//CUT_CHECK_ERROR("Erosus execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());



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
	Susp << <gridDim, blockDim, 0 >> >(nx, ny, dx, Param.eps, Param.nuh, Param.nuhfac, Param.rho, Param.sus, Param.bed, ueu_g, vev_g, uu_g, uv_g, hu_g, vv_g, vu_g, hv_g, zb_g, hh_g, DR_g, Cc_g, ceqbg_g, Sus_g, Svs_g, Sub_g, Svb_g, thetamean_g, ua_g);
	//CUT_CHECK_ERROR("Susp execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	//Calculate suspended concentration
	//

	CUDA_CHECK(cudaMalloc((void **)&ero_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&depo_g, nx*ny*sizeof(DECNUM)));
	Conc << <gridDim, blockDim, 0 >> >(nx, ny, dx, Param.dt, Param.eps, hh_g, Cc_g, ceqsg_g, Tsg_g, facero_g, ero_g, depo_g, Sus_g, Svs_g);
	//CUT_CHECK_ERROR("Conc execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	//
	// Update global variables and fix bnds
	//
	CClatbnd << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, hh_g, Cc_g);
	//CUT_CHECK_ERROR("CClatbnd execution failed\n");
	CUDA_CHECK(cudaDeviceSynchronize());


	if (Param.morfac > 0.0f)// Only if morphology is need i.e. if also if morphac >0.0
	{
		//
		// Adjust sediment fluxes for rocklayer
		//

		CUDA_CHECK(cudaMalloc((void **)&Sout_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMalloc((void **)&indSub_g, nx*ny*sizeof(int)));
		CUDA_CHECK(cudaMalloc((void **)&indSvb_g, nx*ny*sizeof(int)));
		hardlayer << <gridDim, blockDim, 0 >> >(nx, ny, dx, Param.dt, Sub_g, Svb_g, Sout_g, indSub_g, indSvb_g);
		//CUT_CHECK_ERROR("hardlayer execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, Sout_g);
		//CUT_CHECK_ERROR("Sout twodimbnd execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());


		//
		// Bed update
		//
		bedupdate << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, dx, Param.dt, Param.morfac, Param.por, hh_g, ero_g, depo_g, Sub_g, Svb_g, Sout_g, indSub_g, indSvb_g, zb_g, dzb_g, stdep_g);
		//CUT_CHECK_ERROR("bedupdate execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());


		//
		// Update lateral bnd	
		//	
		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, zb_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());


		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, stdep_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());


		//
		// Avalanching
		//
		CUDA_CHECK(cudaMalloc((void **)&ddzb_g, nx*ny*sizeof(DECNUM)));
		CUDA_CHECK(cudaMemcpy(ddzb_g, zeros, nx*ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		avalanching << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, dx, Param.dt, Param.por, Param.drydzmax, Param.wetdzmax, Param.maxslpchg, hh_g, zb_g, ddzb_g, stdep_g);
		//CUT_CHECK_ERROR("avalanching execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		//
		// Update Zb for avalanching
		//

		updatezb << <gridDim, blockDim, 0 >> >(nx, ny, dx, Param.dt, zb_g, ddzb_g, dzb_g, zs_g, hh_g, stdep_g);
		//CUT_CHECK_ERROR("avalanching execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		//
		// Update lateral bnd	
		//	
		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, zb_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		//zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,dzb_g);
		////CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		//CUDA_CHECK( cudaDeviceSynchronize() );

		zblatbnd << <gridDim, blockDim, 0 >> >(nx, ny, stdep_g);
		//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		updatezom << <gridDim, blockDim, 0 >> >(nx, ny, Param.cfsand, Param.cfreef, Param.fwsand, Param.fwreef, stdep_g, cfm_g, fwm_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());
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
	totaltime = 0.0;
	nextoutputtime = 0.0;


	int nx, ny;
	float dx, grdalpha;
	double dt;
	
	


	// Reset the log file 
	FILE * flog;
	flog = fopen("XBG_log.txt", "w");
	fclose(flog);
		
	//Logfile header
	time_t rawtime, dstart;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%d-%m-%Y %H:%M:%S", timeinfo);
	std::string strtimenow(buffer);
	write_text_to_log_file("#################################");
	write_text_to_log_file("XBeach_GPU v0.0");
	write_text_to_log_file("#################################");
	write_text_to_log_file("model started at " + strtimenow);
	

	//////////////////////////////////////////////////////
	/////             Read Operational file          /////
	//////////////////////////////////////////////////////
	XBGPUParam XParam;

	std::vector<SLBnd> slbnd;
	std::vector<WindBnd> wndbnd;

	std::ifstream fs("XBG_param.txt");

	if (fs.fail()){
		std::cerr << "XBG_param.txt file could not be opened" << std::endl;
		write_text_to_log_file("ERROR: XBG_param.txt file could not be opened...use this log file to create a file named XBG_param.txt");
		SaveParamtolog(XParam);
		exit(1);
	}

	std::string line;
	while (std::getline(fs, line))
	{
		//std::cout << line << std::endl;

		//Get param or skip empty lines
		if (!line.empty() && line.substr(0, 1).compare("#") != 0)
		{
			XParam = readparamstr(line, XParam);
			//std::cout << line << std::endl;
		}

	}
	fs.close();

	// This is done after reading the boundaries to better constrain endtine when it is not specified
	//XParam = checkparamsanity(XParam);
	//std::cout << XParam.Bathymetryfile << std::endl;

	


	//filename = XParam.Bathymetryfile.c_str();

	wdt = 0.0;

	FILE * fid;
	FILE * fiz;

	std::string bathyext;
	//read bathy input data:
	if (!XParam.Bathymetryfile.empty())
	{
		printf("bathy: %s\n", XParam.Bathymetryfile.c_str());

		write_text_to_log_file("bathy: " + XParam.Bathymetryfile);

		std::vector<std::string> extvec = split(XParam.Bathymetryfile, '.');
		bathyext = extvec.back();
		write_text_to_log_file("bathy extension: " + bathyext);
		if (bathyext.compare("md") == 0)
		{
			write_text_to_log_file("Reading 'md' file");
			readbathyHead(XParam.Bathymetryfile, XParam.nx, XParam.ny, XParam.dx, XParam.grdalpha);
		}
		if (bathyext.compare("nc") == 0)
		{
			write_text_to_log_file("Reading 'nc' file");
			readgridncsize(XParam.Bathymetryfile, XParam.nx, XParam.ny, XParam.dx);
			write_text_to_log_file("For nc of bathy file please specify grdalpha in the XBG_param.txt (default 0)");
			
		}
		if (bathyext.compare("dep") == 0 || bathyext.compare("bot") == 0)
		{
			//XBeach style file
			write_text_to_log_file("Reading " + bathyext + " file");
			write_text_to_log_file("For this type of bathy file please specify nx, ny, dx and grdalpha in the XBG_param.txt");
		}
			
		XParam.grdalpha = XParam.grdalpha*pi / 180; // grid rotation

		//fid = fopen(XParam.Bathymetryfile.c_str(), "r");
		//fscanf(fid, "%u\t%u\t%lf\t%*f\t%lf", &XParam.nx, &XParam.ny, &XParam.dx, &XParam.grdalpha);
		printf("nx=%d\tny=%d\tdx=%f\talpha=%f\n", XParam.nx, XParam.ny, XParam.dx, XParam.grdalpha*180/pi);
		write_text_to_log_file("nx=" + std::to_string(XParam.nx) + " ny=" + std::to_string(XParam.ny) + " dx=" + std::to_string(XParam.dx) + " grdalpha=" + std::to_string(XParam.grdalpha*180.0/pi));
	}
	else
	{
		std::cerr << "Fatal error:No bathymetry file specified. Please specify using 'bathy = Filename.bot'" << std::endl;
		write_text_to_log_file("Fatal error : No bathymetry file specified.Please specify using 'bathy = Filename.bot'");
		exit(1);
	}
	

	//printf("output time step=%f\n", XParam.outputtimestep);

	//READ INITIAL ZS CONDITION
	//fiz=fopen("zsinit.md","r");
	//fscanf(fiz,"%u\t%u\t%f\t%*f\t%f",&nx,&ny,&dx,&grdalpha);

	//XParam.grdalpha = XParam.grdalpha*pi / 180; // grid rotation

	printf("Opening sea level bnd...");
	write_text_to_log_file("Opening sea level bnd...");

	if (!XParam.slbnd.empty())
	{
		slbnd = readWLfile(XParam.slbnd);
	}
	else
	{
		printf("No file was specified. Setting Offshore water level boundary to zero...");
		write_text_to_log_file("WARNING: No file was specified. Setting Offshore water level boundary to zero...");
		if (XParam.endtime < 0.000001)
		{
			printf("WARNING: endtime is also set to zero so the model will initialise but no calculation will be performed. either specifiy a sea level boundary or endtime");
			write_text_to_log_file("WARNING: endtime is also set to zero so the model will initialise but no calculation will be performed. either specifiy a sea level boundary or endtime");
			
		}
		SLBnd slbndline;

		slbndline.time = 0.0;
		slbndline.wlev = 0.0;
		slbnd.push_back(slbndline);

		slbndline.time = XParam.endtime;
		slbndline.wlev = 0.0;
		slbnd.push_back(slbndline);
	}
	
	
	//Note: the first rtsl should be 0 
	
	SLstepinbnd = 1;

	printf("done\n");
	write_text_to_log_file("done");

	// Read Wind forcing
	printf("Opening wind forcing...");
	write_text_to_log_file("Opening wind forcing...");
	if (!XParam.windfile.empty())
	{
		wndbnd = readWNDfile(XParam.windfile, XParam.grdalpha);
	}
	else
	{
		printf("No wind file was specified. Setting wind forcing to zero...");
		write_text_to_log_file("No wind file was specified. Setting wind forcing to zero...");
		WindBnd wndbndline;

		wndbndline.time = 0.0;
		wndbndline.spd = 0.0;
		wndbndline.dir = 0.0;
		wndbndline.theta = 0.0;
		wndbndline.U = 0.0;
		wndbndline.V = 0.0;
		wndbnd.push_back(wndbndline);

		wndbndline.time = max(XParam.endtime,slbnd.back().time);
		wndbndline.spd = 0.0;
		wndbndline.dir = 0.0;
		wndbndline.theta = 0.0;
		wndbndline.U = 0.0;
		wndbndline.V = 0.0;
		wndbnd.push_back(wndbndline);
		
	}
	WNDstepinbnd = 1;

	printf("done\n");
	write_text_to_log_file("done");
	XParam = checkparamsanity(XParam, slbnd,wndbnd);
	



	nx = XParam.nx;
	ny = XParam.ny;
	dx = XParam.dx;

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
	dummy = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	umeanbnd = (DECNUM *)malloc(ny*sizeof(DECNUM));
	
	
	

	// set initital condition and read bathy file
	printf("Read bathy data...");
	write_text_to_log_file("Read bathy data...");

	if (bathyext.compare("md") == 0)
	{
		readbathy(XParam.Bathymetryfile, zb);
	}
	if (bathyext.compare("nc") == 0)
	{
		readnczb(XParam.nx, XParam.ny, XParam.Bathymetryfile, zb);
	}
	if (bathyext.compare("bot") == 0 || bathyext.compare("dep") == 0)
	{
		readXBbathy(XParam.Bathymetryfile, XParam.nx, XParam.ny, zb);
	}

	

	// set initital condition and read bathy file
	printf("Set initial condition...");
	write_text_to_log_file("Set initial condition...");

	int jread;
	//int jreadzs;
	for (int fnod = ny; fnod >= 1; fnod--)
	{

		//fscanf(fid, "%u", &jread);
		//fscanf(fiz,"%u",&jreadzs);
		umeanbnd[(fnod - 1)] = 0.0f;
		for (int inod = 0; inod < nx; inod++)
		{
			//fscanf(fid, "%f", &zb[inod + (jread - 1)*nx]);
			uu[inod + (fnod - 1)*nx] = 0.0f;
			vv[inod + (fnod - 1)*nx] = 0.0f;
			dzb[inod + (fnod - 1)*nx] = 0.0f;
			cfm[inod + (fnod - 1)*nx] = XParam.cf;
			stdep[inod + (fnod - 1)*nx] = 0.0f;
			zeros[inod + (fnod - 1)*nx] = 0.0f;
			//fscanf(fiz,"%f",&zs[inod+(jreadzs-1)*nx]);

			//hh[inod+(jread-1)*nx]=max(zb[inod+(jread-1)*nx]+zs[inod+(jreadzs-1)*nx],eps);
			//zs[inod+(jread-1)*nx]=max(zs[inod+(jreadzs-1)*nx],-1*zb[inod+(jread-1)*nx]);

			zs[inod + (fnod - 1)*nx] = max(slbnd[0].wlev*1.0f, -1 * zb[inod + (fnod - 1)*nx]);
			hh[inod + (fnod - 1)*nx] = max(zb[inod + (fnod - 1)*nx] + slbnd[0].wlev, XParam.eps);



		}
	}

	//fclose(fid);
	printf("...done\n");
	write_text_to_log_file("...done");


	//// read Hard layer file

	if (!XParam.SedThkfile.empty())
	{
		printf("Hard layer file found\n");
		write_text_to_log_file("Hard layer file found");
		int STnx, STny;
		double STdx, STgrdalpha;

		readbathyHead(XParam.SedThkfile, STnx, STny, STdx, STgrdalpha);
		readbathy(XParam.SedThkfile, stdep);
		if (STnx != nx || STny != ny)
		{
			printf("Error Sediment thickness file (Hard layer file) dimension mismatch. Model will run with constant sediment thickness.\n");
			write_text_to_log_file("ERROR: Sediment thickness file (Hard layer file) dimension mismatch. Model will run with constant sediment thickness.");
		}

		
		
		
	}
	else
	{
		printf("No hard layer file found, Model will run with constant sediment thickness\n");
		write_text_to_log_file("No hard layer file found, Model will run with constant sediment thickness");

		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++)
			{
				stdep[i + j*nx] = 5.0f;
			}
		}
	}

	
	
	//calculate coriolis force
	XParam.lat = XParam.lat*pi / 180.0f;
	DECNUM wearth = pi*(1.0f / 24.0f) / 1800.0f;
	XParam.fc = 2.0f*wearth*sin(XParam.lat);
	
	roller = 1; // option to turn off/on roller model (0/1)

	
	double t1 = -(pi) / 2.0;
	
	//Allocate More array on CPU

	Fx = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	Fy = (DECNUM *)malloc(nx*ny*sizeof(DECNUM));
	zsbnd = (DECNUM *)malloc(ny*sizeof(DECNUM));

	

	cgx = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cgy = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cx = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cy = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	ctheta = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	


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


	if (XParam.GPUDEVICE >= 0)
	{
		// Init GPU

		CUDA_CHECK(cudaSetDevice(XParam.GPUDEVICE));


		if (XParam.swave == 1)
		{
			waveinitGPU(XParam);
		}

		//CUT_DEVICE_INIT(argc, argv);
		printf("Allocating GPU memory\n");
		write_text_to_log_file("Allocating GPU memory");



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

		if (XParam.swave == 1 )
		{
			waveinitGPU(XParam);
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

	if (XParam.GPUDEVICE >= 0)
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


	if (XParam.GPUDEVICE >= 0)
	{
		if (XParam.swave == 1 )
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
		write_text_to_log_file("gridDim= " + std::to_string(gridDim.x) + "," + std::to_string(gridDim.y) + "," + std::to_string(gridDim.z));
		//printf("gridDim=%i,%i,%i\n", gridDimLine.x, gridDimLine.y, gridDimLine.z);

		//Calculate bottomm friction based on initial hard layer file
		updatezom << <gridDim, blockDim, 0 >> >(nx, ny, XParam.cfsand, XParam.cfreef, XParam.fwsand, XParam.fwreef, stdep_g, cfm_g, fwm_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");`
		CUDA_CHECK(cudaDeviceSynchronize());


		// Calculate initial maximum timestep


		FLOWDT << <gridDim, blockDim, 0 >> >(nx, ny, dx, 0.5*XParam.CFL, dtflow_g, hh_g);
		CUDA_CHECK(cudaDeviceSynchronize());


		minmaxKernel << <gridDimLine, blockDimLine, 0 >> >(nx*ny, arrmax_g, arrmin_g, dtflow_g);
		//CUT_CHECK_ERROR("UpdateZom execution failed\n");
		CUDA_CHECK(cudaDeviceSynchronize());

		finalminmaxKernel << <1, blockDimLine, 0 >> >(arrmax_g, arrmin_g);
		CUDA_CHECK(cudaDeviceSynchronize());

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



		dt = arrmin[0]*0.5;

		printf("Initial timestep: dt=%f\n", dt);
		write_text_to_log_file("Initial timestep: dt=" + std::to_string(dt));

		double tiny = 0.00000001;

		if (dt < tiny)
		{
			std::cerr << " Error: timestep too small;" << std::endl;
			write_text_to_log_file("ERROR: timestep too small");
			exit(EXIT_FAILURE);
		}

		XParam.dt = dt;
		
		if (XParam.swave == 1)
		{
			set_bnd << <gridDim, blockDim, 0 >> >(nx, ny, Trep, ntheta, theta_g, sigm_g);

			//CUT_CHECK_ERROR("set_bnd() execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

			//also run dispersion relation, cg is needed in the first iteration of the bnd flow
			dispersion_init << <gridDim, blockDim, 0 >> >(nx, ny, twopi, XParam.g, aphi, bphi, sigm_g, hh_g, cg_g);
			//CUT_CHECK_ERROR("dispersion execution failed\n");
			CUDA_CHECK(cudaDeviceSynchronize());

		}
		
	}
	else
	{


		//Calculate bottomm friction based on initial hard layer file
		updatezomCPU(nx, ny, XParam.cfsand, XParam.cfreef, XParam.fwsand, XParam.fwreef, stdep_g, cfm_g, fwm_g);


		if (XParam.swave == 1 )
		{
			set_bndCPU(nx, ny, Trep, ntheta, theta_g, sigm_g);



			//also run dispersion relation, cg is needed in the first iteration of the bnd flow
			dispersion_initCPU(nx, ny, twopi, g, aphi, bphi, sigm_g, hh_g, cg_g);


		}
	}

	//Map a link between the variable s a string and the variable pointer
	//{ "hh", "uu", "vv", "zs", "zb", "cfm", "dzb", "stdep", "Fx", "Fy", "cgx", "cgy", "cx", "cy", "ctheta", "usd", "D", "E", "H", "urms", "ueu", "vev", "thetamean", "Hmean", "uumean", "vvmean", "hhmean", "zsmean", "Cmean", "sigm", "k", "c", "kh", "cg", "sinh2kh", "dhdx", "dhdy", "dudx", "dudy", "dvdx", "dvdy", "C", "R", "DR" };
	std::map<std::string, DECNUM *> OutputVarMapCPU;
	std::map<std::string, int> OutputVarMapndim;

	
	OutputVarMapCPU["hh"] = hh;
	OutputVarMapndim["hh"] = 3;

	OutputVarMapCPU["uu"] = uu;
	OutputVarMapndim["uu"] = 3;

	OutputVarMapCPU["vv"] = vv;
	OutputVarMapndim["vv"] = 3;

	OutputVarMapCPU["zs"] = zs;
	OutputVarMapndim["zs"] = 3;
	
	OutputVarMapCPU["zb"] = zb;
	if (XParam.morphology == 1)
	{
		OutputVarMapndim["zb"] = 3;
	}
	else
	{
		OutputVarMapndim["zb"] = 2;
	}
	
	OutputVarMapCPU["cfm"] = cfm;
	OutputVarMapndim["cfm"] = 3;

	OutputVarMapCPU["dzb"] = dzb;
	OutputVarMapndim["dzb"] = 3;

	OutputVarMapCPU["stdep"] = stdep;
	OutputVarMapndim["stdep"] = 3;

	OutputVarMapCPU["Fx"] = Fx;
	OutputVarMapndim["Fx"] = 3;

	OutputVarMapCPU["Fy"] = Fy;
	OutputVarMapndim["Fy"] = 3;

	OutputVarMapCPU["cgx"] = cgx;
	OutputVarMapndim["cgx"] = 4;

	OutputVarMapCPU["cgy"] = cgy;
	OutputVarMapndim["cgy"] = 4;

	OutputVarMapCPU["cx"] = cx;
	OutputVarMapndim["cx"] = 4;

	OutputVarMapCPU["cy"] = cy;
	OutputVarMapndim["cy"] = 4;

	OutputVarMapCPU["ctheta"] = ctheta;
	OutputVarMapndim["ctheta"] = 4;

	OutputVarMapCPU["ee"] = ee;
	OutputVarMapndim["ee"] = 4;

	OutputVarMapCPU["rr"] = rr;
	OutputVarMapndim["rr"] = 4;

	OutputVarMapCPU["D"] = D;
	OutputVarMapndim["D"] = 3;

	OutputVarMapCPU["E"] = E;
	OutputVarMapndim["E"] = 3;

	OutputVarMapCPU["H"] = H;
	OutputVarMapndim["H"] = 3;

	OutputVarMapCPU["thetamean"] = thetamean;
	OutputVarMapndim["thetamean"] = 3;

	OutputVarMapCPU["urms"] = urms;
	OutputVarMapndim["urms"] = 3;

	OutputVarMapCPU["ueu"] = ueu;
	OutputVarMapndim["ueu"] = 3;

	OutputVarMapCPU["vev"] = vev;
	OutputVarMapndim["vev"] = 3;


	OutputVarMapCPU["hhmean"] = hhmean;
	OutputVarMapndim["hhmean"] = 3;

	OutputVarMapCPU["uumean"] = uumean;
	OutputVarMapndim["uumean"] = 3;

	OutputVarMapCPU["vvmean"] = vvmean;
	OutputVarMapndim["vvmean"] = 3;

	OutputVarMapCPU["zsmean"] = zsmean;
	OutputVarMapndim["zsmean"] = 3;

	OutputVarMapCPU["Hmean"] = Hmean;
	OutputVarMapndim["Hmean"] = 3;

	OutputVarMapCPU["Cmean"] = Cmean;
	OutputVarMapndim["Cmean"] = 3;

	OutputVarMapCPU["sigm"] = sigm;
	OutputVarMapndim["sigm"] = 3;

	OutputVarMapCPU["k"] = k;
	OutputVarMapndim["k"] = 3;

	OutputVarMapCPU["c"] = c;
	OutputVarMapndim["c"] = 3;

	OutputVarMapCPU["kh"] = kh;
	OutputVarMapndim["kh"] = 3;

	OutputVarMapCPU["cg"] = cg;
	OutputVarMapndim["cg"] = 3;

	OutputVarMapCPU["sinh2kh"] = sinh2kh;
	OutputVarMapndim["sinh2kh"] = 3;

	OutputVarMapCPU["dhdx"] = dhdx;
	OutputVarMapndim["dhdx"] = 3;

	OutputVarMapCPU["dhdy"] = dhdy;
	OutputVarMapndim["dhdy"] = 3;

	OutputVarMapCPU["dudx"] = dudx;
	OutputVarMapndim["dudx"] = 3;

	OutputVarMapCPU["dudy"] = dudy;
	OutputVarMapndim["dudy"] = 3;

	OutputVarMapCPU["dvdx"] = dvdx;
	OutputVarMapndim["dvdx"] = 3;

	OutputVarMapCPU["dvdy"] = dvdy;
	OutputVarMapndim["dvdy"] = 3;

	OutputVarMapCPU["C"] = C;
	OutputVarMapndim["C"] = 3;

	OutputVarMapCPU["R"] = R;
	OutputVarMapndim["R"] = 3;

	OutputVarMapCPU["DR"] = DR;
	OutputVarMapndim["DR"] = 3;

	OutputVarMapCPU["wci"] = zeros;
	OutputVarMapndim["wci"] = 3;

	OutputVarMapCPU["vmageu"] = zeros;
	OutputVarMapndim["vmageu"] = 3;

	OutputVarMapCPU["vmagev"] = zeros;
	OutputVarMapndim["vmagev"] = 3;

	OutputVarMapCPU["dzsdx"] = zeros;
	OutputVarMapndim["dzsdx"] = 3;

	OutputVarMapCPU["dzsdy"] = zeros;
	OutputVarMapndim["dzsdy"] = 3;

	OutputVarMapCPU["dzsdt"] = zeros;
	OutputVarMapndim["dzsdt"] = 3;

	OutputVarMapCPU["fwm"] = zeros;
	OutputVarMapndim["fwm"] = 3;

	OutputVarMapCPU["hu"] = hh;
	OutputVarMapndim["hu"] = 3;

	OutputVarMapCPU["hv"] = hh;
	OutputVarMapndim["hv"] = 3;

	OutputVarMapCPU["hum"] = hh;
	OutputVarMapndim["hum"] = 3;

	OutputVarMapCPU["hvm"] = hh;
	OutputVarMapndim["hvm"] = 3;

	OutputVarMapCPU["uv"] = zeros;
	OutputVarMapndim["uv"] = 3;

	OutputVarMapCPU["vu"] = zeros;
	OutputVarMapndim["vu"] = 3;

	OutputVarMapCPU["ududx"] = zeros;
	OutputVarMapndim["ududx"] = 3;

	OutputVarMapCPU["vdvdy"] = zeros;
	OutputVarMapndim["vdvdy"] = 3;

	OutputVarMapCPU["vdudy"] = zeros;
	OutputVarMapndim["vdudy"] = 3;

	OutputVarMapCPU["udvdx"] = zeros;
	OutputVarMapndim["udvdx"] = 3;

	OutputVarMapCPU["ust"] = zeros;
	OutputVarMapndim["ust"] = 3;

	OutputVarMapCPU["stdep"] = stdep;
	OutputVarMapndim["stdep"] = 3;

	OutputVarMapCPU["kturb"] = zeros;
	OutputVarMapndim["kturb"] = 3;

	OutputVarMapCPU["rolthick"] = zeros;
	OutputVarMapndim["rolthick"] = 3;

	OutputVarMapCPU["ceqsg"] = zeros;
	OutputVarMapndim["ceqsg"] = 3;

	// prepare output file
	printf("prepare output");
	write_text_to_log_file("prepare output");

	// Proof of concept for map
	//creatncfile(tsoutfile, nx, ny, dx, 0.0f, imodel, zb, zs, uu, vv, H, H, thetamean, uu, uu, uu, uu, uu, uu, uu, hh, uu, uu, uu, uu, uu, uu);
	//creatncfile(XParam, 0.0f, OutputVarMapCPU["zb"], OutputVarMapCPU["zs"], OutputVarMapCPU["uu"], OutputVarMapCPU["vv"], H, H, thetamean, uu, uu, uu, uu, uu, uu, uu, hh, uu, uu, uu, uu, uu, uu);

	if (!XParam.outvars.empty())
	{
		//create nc file with no variables
		creatncfileUD(XParam, 0.0, ntheta, dtheta, thetamin, thetamax);
		for (int ivar = 0; ivar < XParam.outvars.size(); ivar++)
		{
			//Create definition for each variable and store it
			defncvar(XParam, XParam.outvars[ivar], OutputVarMapndim[XParam.outvars[ivar]], OutputVarMapCPU[XParam.outvars[ivar]]);
		}
	}



	//create3dnc(nx,ny,ntheta,dx,0.0f,theta,ctheta);

	istepout = istepout + nstepout; // Depreciated ?
	nextoutputtime = nextoutputtime + XParam.outputtimestep;

	printf("...done\n");
	write_text_to_log_file("...done");


	write_text_to_log_file("#################################");
	write_text_to_log_file("XBeach_GPU Parameters:");

	SaveParamtolog(XParam);
	write_text_to_log_file("#################################");
	write_text_to_log_file("Starting Computation ");
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
	if (XParam.GPUDEVICE >= 0)
	{
		mainloopGPU(XParam, slbnd,wndbnd);
	}
	else
	{
		mainloopCPU(XParam, slbnd,wndbnd);
	}




	//close the bnd files and clean up a bit
	//fclose(fsl);
	//fclose(fwind);

	endcputime = clock();
	printf("End Computation");
	printf("Total runtime= %d  seconds\n", (endcputime - startcputime) / CLOCKS_PER_SEC);
	write_text_to_log_file("End Computation");
	write_text_to_log_file("#################################");
	write_text_to_log_file("Total runtime= " + std::to_string((endcputime - startcputime) / CLOCKS_PER_SEC) + " seconds");
	cudaDeviceReset();

	//CUT_EXIT(argc, argv);










}

