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

//#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cutil.h>
#include <string>
#include <netcdf.h>
#include <time.h>



using namespace std;


#include "Wave_kernel.cu"
#include "Flow_kernel.cu"
#include "Sediment_kernel.cu"
#include "read_input.cpp"

#define pi 3.14159265


// additional functions
void makjonswap(float hm0gew,float fp,float mainang,float rt,float scoeff,float gam,float * theta,int ntheta,float& TTrep,float * &Stt);
extern "C" void creatncfile(char outfile[], int nx,int ny,/*int npart,*/float dx,float totaltime,int imodel,/*float * xxp,float * yyp,*/float *zb,float *zs,float * uu, float * vv, float * H,float * Tp,float * Dp,float * D,float * Urms,float * ueu,float * vev,float * C,float *Fx,float *Fy,float *hh,float *Hmean,float *uumean,float *vvmean,float *hhmean,float *zsmean,float *Cmean);
extern "C" void writestep2nc(char outfile[], int nx,int ny,/*int npart,*/float totaltime,int imodel/*,float *xxp,float *yyp*/,float *zb,float *zs,float * uu, float * vv, float * H, float * Tp, float *Dp,float *D,float *Urms,float *ueu,float *vev,float *C,float *dzb,float *Fx,float *Fy,float *hh,float *Hmean,float *uumean,float *vvmean,float *hhmean,float *zsmean,float *Cmean);

extern "C" void create3dnc(int nx,int ny,int nt,float dx,float totaltime,float *theta,float * var);
extern "C" void write3dvarnc(int nx,int ny,int nt,float totaltime,float * var);

extern "C" void read3Dnc(int nx, int ny,int ntheta,char ncfile[],float * &ee);
extern "C" void read2Dnc(int nx, int ny,char ncfile[],float * &hh);

extern "C" void readXbbndhead(char * wavebndfile,float &thetamin,float &thetamax,float &dtheta,float &dtwavbnd,int &nwavbnd,int &nwavfile);
extern "C" void readXbbndstep(int nx, int ny,int ntheta,char * wavebndfile,int step,float &Trep,double *&qfile,double *&Stfile );
extern "C" void readStatbnd(int nx, int ny,int ntheta,float rho,float g,char * wavebndfile,double *&Tpfile,double *&Stfile );
extern "C" void readbndhead(char * wavebndfile,float &thetamin,float &thetamax,float &dtheta,float &dtwavbnd,int &nwavbnd);

//global variables

float Trep,Trepold,Trepnew;
float * St,* Stnew,* Stold;
double * Stfile;
double * qfile;
double * Tpfile;
int nwbndstep=0;
int wavebndtype;
int nstep=0;
int breakmod=1;

float * hh;//=ones(20,40).*10; //Depth+SL
float * zb,* qbndold,* qbndnew;
float * qbndold_g,*qbndnew_g;
float * umeanbnd_g;
float * vmeanbnd_g;
float * umeanbnd;
float * hh_g,*uu_g,*vv_g,*zs_g,* zb_g,*hhold_g;
float * ueu_g, * vev_g;
float * vmageu_g,* vmagev_g;
float * uu;
float  *vv;
float *zs;
float *dummy;

float *xadvec_g,*yadvec_g,*thetaadvec_g;

float *hum_g,* hu_g;
float *hvm_g,* hv_g;
int *wetu_g,*wetv_g;
float *ududx_g,* vdudy_g;
float *udvdx_g, *vdvdy_g;
float *vu_g,*uv_g;
float  nuh,nuhfac;
float *nuh_g;
float * viscu_g,* viscv_g;

float uumin=0.0f;


int nx,ny;
float dx,dt,eps;
float grdalpha;
float totaltime=0.0f;
int nstpw,nwstp;//nb of hd step between wave step and next step for calculating waves and
float wdt;// wave model time step
float wavbndtime;
float slbndtime;
float windtime;
float Cd; //Wind drag
float fp,hm0gew,mainang,rt,scoeff,gam;
int nwavbnd,nwavfile;
float dtwavbnd;
int roller;
float wci;
float *wci_g;
float gammax,hwci,gammaa,n,alpha,beta,t1;
float fw,fw2;
float *fwm_g;//Wave dissipation factor map

float phi    = (1.0f + sqrt(5.0f))/2;
float aphi   = 1/(phi+1);
float bphi   = phi/(phi+1);
float twopi  = 8*atan(1.0f);

float g=9.81f;
float rho=1025.0f;
float zo;
float cf,cf2;//friction
float *cfm, *cfm_g; //friction map

float lat; //lattitude 
float fc; //coriolis

int ntheta;
float thetamin,thetamax;
float *theta;//=(0:100)*((pi)/100)+t1;see below
float *theta_g;

float * cgx, * cgy, *cx, *cy, * ctheta,* cxsth,* sxnth;//
float * cgx_g, * cgy_g, *cx_g, *cy_g, * ctheta_g,* cxsth_g,* sxnth_g,*eect_g;//

int var2plot=2;// 1: wave height 2: eta 3: u 4: v
int colorindx;

float dang;
float dtheta;

float * ee;//=zeros(nx+1,ny+1,ntheta);
float * dd;
float * wete;

float * ee_g, *St_g;

float * tm_g;

float * rr;//;rr=zeros(nx+1,ny+1,ntheta);

float * drr;//=zeros(nx+1,ny+1,ntheta);
float * usd;//=zeros(nx+1,ny+1);
float * D;//=zeros(nx+1,ny+1);
float * D_g;
float * E, * H;
float * E_g, * H_g;
float * drr_g, * rr_g;
float * DR_g, * R_g;

float * DR, * R;

float * Sxx_g,* Syy_g,* Sxy_g,* Fx_g,* Fy_g;
float /** Sxx,* Syy,* Sxy,*/* Fx,* Fy;
float * thetamean;
float * thetamean_g;

float * urms_g;
float * ust_g;
float omega;// = 2*pi/Trep;

float D50, D90,rhosed;

float * kturb_g,* rolthick_g,*dzsdt_g;
float * Ceq_g,* ceqsg_g,* ceqbg_g, * Tsg_g, * facero_g;
float * C,* Cc_g,* stdep, * stdep_g;
float * Sus_g, * Svs_g, * Sub_g,* Svb_g;
float morfac,por;
float * ero_g, *depo_g;
float *dzb,*dzb_g,* ddzb_g;
float * Sout_g;
int * indSub_g,* indSvb_g;
float sus=1.0f;
float bed=1.0f;
float facsk, facas;

float * zsbnd;
float rtsl;
float zsbndnew, zsbndold;

float windth,windthold,windv,windvold,windvnew,windthnew,rtwind;
FILE * fsl;
FILE * fwav;
FILE * fwind;

FILE * Tsout;
//FILE * fXq,* fXE;

char tsoutfile[256];
int iout,jout;
int imodel; //1: Wave only; 2: Current Only; 3: Wave + current; 4: Wave + current + Sediment transport
int nstepout=0;//nb of step between outputs. 0= no output
int istepout=0;//
int nstepplot=0;//nb of step between plots. 0= no plotting
int istepplot=1;
int displayon=0;
int endstep;



float wws;
float drydzmax;
float wetdzmax;
float maxslpchg;



float Hplotmax;

float * sigm, * thet, * costhet, * sinthet;
float * sigm_g, * thet_g, * costhet_g, * sinthet_g;
float * ueu,*vev,*urms;
float * ua_g;

float *k,*c,*kh,*cg,*sinh2kh;
float *k_g, * c_g, * kh_g, * cg_g, * sinh2kh_g;

float * dhdx, * dhdy, * dudx, * dudy, * dvdx, * dvdy;
float * dhdx_g, * dhdy_g, * dudx_g, * dudy_g, * dvdx_g, * dvdy_g;
float *dzsdx_g,* dzsdy_g;
float *zeros;

float * Hmean_g,* uumean_g,*vvmean_g,*hhmean_g,*zsmean_g,*Cmean_g;
float *Hmean,*uumean,*vvmean,*hhmean,*zsmean,*Cmean;

int GPUDEVICE;

int startflowstep;
int usesmago;

// Particle stuff
int npart;
int sedstart;
int wxstep=1;


char wavebndfile[256];

extern "C"
void wavebnd(void);
void flowbnd(void);
void wavestep(void);
void flowstep(void);
void sedimentstep(void);


// Main loop that actually runs the model
void mainloop(void)
{

while (nstep<=endstep)
{
	
	nstep++; 
	wdt=dt; // Sometinmes in stationary wave run one can have a larger wdt (wave time step)
	totaltime=nstep*dt;	//total run time acheived until now in s

	dim3 blockDim(16, 16, 1);// This means that the grid has to be a factor of 16 on both x and y
	dim3 gridDim(nx / blockDim.x, ny / blockDim.y, 1);
	
	if(imodel==1 || imodel>2)
	{
		wavebnd(); // Calculate the boundary condition for this step
		wavestep(); // Calculate the wave action ballance for this step
	}


	
	if(imodel>=2)
	{
		flowbnd();// Calculate the flow boundary for this step
		flowstep();// solve the shallow water and continuity for this step
	}
	if(imodel>=4 && nstep>=sedstart)
	{
	//Sediment step
	sedimentstep();//solve the sediment dispersion, and morphology
	}

	//add last value for avg calc
	addavg_var<<<gridDim, blockDim, 0>>>(nx,ny,Hmean_g,H_g);
	CUT_CHECK_ERROR("Add avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	addavg_var<<<gridDim, blockDim, 0>>>(nx,ny,uumean_g,uu_g);
	CUT_CHECK_ERROR("Add avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	addavg_var<<<gridDim, blockDim, 0>>>(nx,ny,vvmean_g,vv_g);
	CUT_CHECK_ERROR("Add avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	addavg_var<<<gridDim, blockDim, 0>>>(nx,ny,hhmean_g,hh_g);
	CUT_CHECK_ERROR("Add avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	addavg_var<<<gridDim, blockDim, 0>>>(nx,ny,zsmean_g,zs_g);
	CUT_CHECK_ERROR("Add avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	addavg_var<<<gridDim, blockDim, 0>>>(nx,ny,Cmean_g,Cc_g);
	CUT_CHECK_ERROR("Add avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	
	if (nstep==istepout && nstepout>0)
	{
	istepout=istepout+nstepout;

	//Avg mean variables

	divavg_var<<<gridDim, blockDim, 0>>>(nx,ny,nstepout,Hmean_g);
	CUT_CHECK_ERROR("Div avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	divavg_var<<<gridDim, blockDim, 0>>>(nx,ny,nstepout,uumean_g);
	CUT_CHECK_ERROR("Div avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	divavg_var<<<gridDim, blockDim, 0>>>(nx,ny,nstepout,vvmean_g);
	CUT_CHECK_ERROR("Div avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	divavg_var<<<gridDim, blockDim, 0>>>(nx,ny,nstepout,hhmean_g);
	CUT_CHECK_ERROR("Div avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	divavg_var<<<gridDim, blockDim, 0>>>(nx,ny,nstepout,zsmean_g);
	CUT_CHECK_ERROR("Div avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	divavg_var<<<gridDim, blockDim, 0>>>(nx,ny,nstepout,Cmean_g);
	CUT_CHECK_ERROR("Div avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	// Download mean vars
	CUDA_SAFE_CALL( cudaMemcpy(Hmean, Hmean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(uumean, uumean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(vvmean, vvmean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(hhmean, hhmean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(zsmean, zsmean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(Cmean, Cmean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	

	CUDA_SAFE_CALL( cudaMemcpy(H, H_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(uu, uu_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(vv, vv_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(zs, zs_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(Fx, Fx_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(Fy, Fy_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(thetamean, thetamean_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(D, D_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(urms, urms_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(ueu, ueu_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(vev, vev_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	//CUDA_SAFE_CALL( cudaMemcpy(C, ceqsg_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(C, Cc_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	//CUDA_SAFE_CALL( cudaMemcpy(C,k_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	//CUDA_SAFE_CALL( cudaMemcpy(ctheta,ee_g, nx*ny*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(hh, hh_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	if(imodel==4)// If moprhology is on
	{
		CUDA_SAFE_CALL( cudaMemcpy(zb, zb_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
		CUDA_SAFE_CALL( cudaMemcpy(dzb, dzb_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );
	}
	//CUDA_SAFE_CALL( cudaMemcpy(xxp, xxp_g, npart*sizeof(float ), cudaMemcpyDeviceToHost) );
	//CUDA_SAFE_CALL( cudaMemcpy(yyp, yyp_g, npart*sizeof(float ), cudaMemcpyDeviceToHost) );
	printf("Writing output, totaltime:%f s\n",totaltime);
	writestep2nc(tsoutfile,nx,ny,/*npart,*/totaltime,imodel,/*xxp,yyp,*/zb,zs,uu, vv, H,H,thetamean,D,urms,ueu,vev,C,dzb,Fx,Fy,hh,Hmean,uumean,vvmean,hhmean,zsmean,Cmean);

	//write3dvarnc(nx,ny,ntheta,totaltime,ctheta);
		   //outfile[],nx,ny,npart,totaltime,xxp,yyp,zs,uu, vv, H,Tp,Dp,      D,Urms,ueu,vev)
	//fprintf(Tsout,"%f\t%f\t%f\t%f\t%f\t%f\n",totaltime,hh[iout+jout*nx],zs[iout+jout*nx],uu[iout+jout*nx],vv[iout+jout*nx],H[iout+jout*nx]);


	//Clear avg vars
	resetavg_var<<<gridDim, blockDim, 0>>>(nx,ny,Hmean_g);
	CUT_CHECK_ERROR("Reset avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	resetavg_var<<<gridDim, blockDim, 0>>>(nx,ny,uumean_g);
	CUT_CHECK_ERROR("Reset avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	resetavg_var<<<gridDim, blockDim, 0>>>(nx,ny,vvmean_g);
	CUT_CHECK_ERROR("Reset avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	resetavg_var<<<gridDim, blockDim, 0>>>(nx,ny,hhmean_g);
	CUT_CHECK_ERROR("Reset avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	resetavg_var<<<gridDim, blockDim, 0>>>(nx,ny,zsmean_g);
	CUT_CHECK_ERROR("Reset avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	resetavg_var<<<gridDim, blockDim, 0>>>(nx,ny,Cmean_g);
	CUT_CHECK_ERROR("Reset avg execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	}
}
}


void waveinit(void)
{
	// Initialize wave model
	
	//Wave input bnd
	printf("Opening wave bnd\n");

	if(wavebndtype==1)
	{
		
		readbndhead(wavebndfile,thetamin,thetamax,dtheta,dtwavbnd,nwavbnd);
		

	}
	else
	{
		readXbbndhead(wavebndfile,thetamin,thetamax,dtheta,dtwavbnd,nwavbnd,nwavfile);
		
	}
	



	//printf("Hs=%f\tTp=%f\ta=%f\tscoef=%f\tgam=%f\trt=%f\n",hm0gew,fp,mainang,scoeff,gam,rt);
	
	fp=1/fp;
	
	





	//hm0gew=1.0f;//Significant wave height (m)
	//fp=1.0f/12.0f; //Wave peak frequency (Hz)
	//mainang=0; //wave mean direction (angle of incidence º)
	//rt=1000; //Boundary duration
	//scoeff=100;// spread coef n.u.
	//gam=3.3f;//: peak enhancement factor, optional parameter (DEFAULT 3.3)
	
	thetamin=thetamin*pi/180;
	thetamax=thetamax*pi/180;
	dtheta=dtheta*pi/180;
	
	ntheta=round((thetamax-thetamin)/dtheta);
	printf("ntheta=%d\tdtheta=%f\n",ntheta,dtheta);
	printf("nwavbnd=%d\n",nwavbnd);

	theta=(float *)malloc(ntheta*sizeof(float));
	
	Stfile=(double *)malloc(ntheta*ny*nwavbnd*sizeof(double));
	qfile=(double *)malloc(3*ny*nwavbnd*sizeof(double));
	Tpfile=(double *)malloc(nwavbnd*sizeof(double));
	//dummy=(double *)malloc(1000*sizeof(double));

	qbndnew=(float *)malloc(3*ny*sizeof(float));
	qbndold=(float *)malloc(3*ny*sizeof(float));
	St= (float *)malloc(ntheta*ny*sizeof(float));
	Stold= (float *)malloc(ntheta*ny*sizeof(float));
	Stnew= (float *)malloc(ntheta*ny*sizeof(float));
	cxsth=(float *)malloc(ntheta*sizeof(float));
	sxnth=(float *)malloc(ntheta*sizeof(float));
	
	for(int i=0; i<ntheta; i++)
	{
		theta[i]=i*(dtheta)+thetamin+0.5f*dtheta;
		cxsth[i]=cosf(theta[i]);
		sxnth[i]=sinf(theta[i]);

		//printf("theta=%f\tcxsth=%f\tsxnth=%f\n",theta[i],cxsth[i],sxnth[i]);
	}

	dang=theta[1]-theta[0];
	//dtheta=dang;
	

	ee=(float *)malloc(nx*ny*ntheta*sizeof(float));
	dd=(float *)malloc(nx*ny*ntheta*sizeof(float));
	wete=(float *)malloc(nx*ny*ntheta*sizeof(float));

	rr=(float *)malloc(nx*ny*ntheta*sizeof(float));
	
	//drr=(float *)malloc(nx*ny*ntheta*sizeof(float));


	printf("Reading bnd data\n");
	if(wavebndtype==1)
	{
		readStatbnd(nx,ny,ntheta,rho,g,wavebndfile,Tpfile,Stfile );
		Trepold=Tpfile[0];
		Trepnew=Tpfile[1];
		rt=dtwavbnd;

	}
	
	if(wavebndtype==2)
	{
		readXbbndstep(nx,ny,ntheta,wavebndfile,1,Trepold,qfile,Stfile );
	

	for (int ni=0; ni<ny; ni++)
		{
		for (int itheta=0; itheta<ntheta; itheta++)
		{
			Stold[ni+itheta*ny]=Stfile[ni+itheta*ny+nwbndstep*ny*ntheta];
			Stnew[ni+itheta*ny]=Stfile[ni+itheta*ny+(nwbndstep+1)*ny*ntheta];
			
			
		}
		for (int xi=0; xi<3; xi++)
		{
			qbndold[ni+xi*ny]=qfile[ni+xi*ny+nwbndstep*ny*3];
			qbndnew[ni+xi*ny]=qfile[ni+xi*ny+(nwbndstep+1)*ny*3];
		}
		}
		//CUDA_SAFE_CALL( cudaMemcpy(qbndold_g,qbndold, 3*ny*sizeof(float ), cudaMemcpyHostToDevice) );
		//CUDA_SAFE_CALL( cudaMemcpy(qbndnew_g,qbndnew, 3*ny*sizeof(float ), cudaMemcpyHostToDevice) );
		//printf("qfile[0]=%f\n",qfile[0]);
	}
	else
	{
		

		for (int ni=0; ni<ny; ni++)
		{
			for (int itheta=0; itheta<ntheta; itheta++)
			{
				Stold[ni+itheta*ny]=Stfile[itheta];
				Stnew[ni+itheta*ny]=Stfile[itheta+ntheta];
			}
	
		}
		
	}



	

	//fscanf(fwav,"%f\t%f\t%f\t%f\t%f\t%f",&hm0gew,&fp,&mainang,&scoeff,&gam,&rt);
	//mainang=(1.5*pi-grdalpha)-mainang*pi/180;
	//fp=1/fp;
	//printf("init rt=%f\n",rt);
	

	//makjonswap(hm0gew,fp,mainang,rt,scoeff,gam,theta,ntheta,Trepnew, Stnew);
	
	
	
	
	Trep=Trepold;
	for (int i=0; i<ntheta; i++)                             //! Fill St
	{
		//St[i]=Stold[i];
		//printf("St[%d]=%f\n",i,St[i]);
		for (int ii=0; ii<ny; ii++)
		{
			St[ii+i*ny]=Stold[ii+i*ny];
		}
	}
	
	







	//printf("hh=%f\n",hh[0]);





for (int ii=0; ii<nx; ii++)
{

	for (int jj=0; jj<ny; jj++)
	{

		for (int nt=0; nt<ntheta; nt++)
		{
			if(ii==0)
			{
				ee[0+jj*nx+nt*nx*ny]=St[jj+nt*ny];// not on gpu since it is a bank conflicting problem
			}
			else{
				ee[ii+jj*nx+nt*nx*ny]=0.0f;
			}
			rr[ii+jj*nx+nt*nx*ny]=0.0f;
			
		}
	}
}
	

	
	
	
	
}



void wavebnd(void)
{
    if (totaltime > dtwavbnd*(nwavbnd-1)*wxstep)
	{
		if (wavebndtype==2)
		{
			readXbbndstep(nx,ny,ntheta,wavebndfile,wxstep,Trep,qfile,Stfile);
			
		}
		nwbndstep=0;
		
		
		
		for (int ni=0; ni<ny; ni++)
			{
			
				for (int itheta=0; itheta<ntheta; itheta++)
				{
					Stold[ni+itheta*ny]=Stfile[ni+itheta*ny+nwbndstep*ny*ntheta];
					Stnew[ni+itheta*ny]=Stfile[ni+itheta*ny+(nwbndstep+1)*ny*ntheta];
				}
				if (wavebndtype==2)
				{
				for (int xi=0; xi<3; xi++)
				{
					qbndold[ni+xi*ny]=qfile[ni+xi*ny+nwbndstep*ny*3];
					qbndnew[ni+xi*ny]=qfile[ni+xi*ny+(nwbndstep+1)*ny*3];
				}
				}
			}
				
		
		wxstep=wxstep+1;

	}
	


	

	
	//if ((nstep==1 || nstep==nwstp) && (imodel==1 || imodel>=3))
	//{
		//update wave bnd
		
		if (totaltime>=wavbndtime /*&& wavebndtype==2*/)
		{
			for (int i=0; i<ntheta; i++)                             //! Fill Stold
			{
				for (int ni=0; ni<ny; ni++)
				{
					Stold[ni+i*ny]=Stfile[ni+i*ny+nwbndstep*ntheta*ny];
					
				}
			}
			//fscanf(fwav,"%f\t%f\t%f\t%f\t%f\t%f",&hm0gew,&fp,&mainang,&scoeff,&gam,&rt);
			//mainang=(1.5*pi-grdalpha)-mainang*pi/180;
			//printf("rt=%f\n",rt);
	
			//fp=1/fp;

			//fscanf(fwav,"%f\t%f",&rt,&Trepnew);
			nwbndstep=nwbndstep+1;
			for (int i=0; i<ntheta; i++)                             //! Fill St
			{
				for (int ni=0; ni<ny; ni++)
				{
					
					
						Stnew[ni+i*ny]=Stfile[ni+i*ny+nwbndstep*ntheta*ny];
					
					if(wavebndtype==1)
					{
						
						Trep=Tpfile[nwbndstep];
					}
				}
			}
			if (wavebndtype==2)
			{
				for (int ni=0; ni<ny; ni++)
					{
					for (int xi=0; xi<3; xi++)
					{
						qbndold[ni+xi*ny]=qbndnew[ni+ny*xi];
						qbndnew[ni+xi*ny]=qfile[ni+xi*ny+nwbndstep*ny*3];
					}	
					}
			    CUDA_SAFE_CALL( cudaMemcpy(qbndold_g,qbndold, 3*ny*sizeof(float ), cudaMemcpyHostToDevice) );
				CUDA_SAFE_CALL( cudaMemcpy(qbndnew_g,qbndnew, 3*ny*sizeof(float ), cudaMemcpyHostToDevice) );
				//printf("qbndold[300]=%f\n",qbndold[300]);
				//printf("qbndnew[300]=%f\n",qbndnew[300]);
			}
			//printf("Stfile[0]=%f\n",Stfile[0]);
			

			//makjonswap(hm0gew,fp,mainang,rt,scoeff,gam,theta,ntheta,Trepnew, Stnew);
		wavbndtime=wavbndtime+dtwavbnd;
		}
		
		for (int i=0; i<ntheta; i++)                             //! Fill St
		{
			for (int ni=0; ni<ny; ni++)
			{
				St[ni+i*ny]=Stold[ni+i*ny]+(totaltime-wavbndtime+dtwavbnd)*(Stnew[ni+i*ny]-Stold[ni+i*ny])/dtwavbnd;
			}
			//printf("St[%d]=%f\n",i,St[i*ny]);
		}
		//printf("Wave timestep:%f\n",wdt);
		//Wave model step
		//wavestep();
		nwstp=nstep+nstpw;
		wdt=dt;
	//}

}

void flowbnd(void)
{
	//update sl bnd

	if (totaltime>=slbndtime)
	{
		zsbndold=zsbndnew;
		rtsl=slbndtime;
		fscanf(fsl,"%f\t%f",&slbndtime,&zsbndnew);		
		//slbndtime=+rtsl;
		//zsbnd=zsbndold+(t-slbndtime+rtsl)*(zsbndnew-zsbndold)/rtsl;
	}





	if (wavebndtype==1)
	{
		for (int ni=0; ni<ny; ni++)
		{
			zsbnd[ni]=zsbndold+(totaltime-rtsl)*(zsbndnew-zsbndold)/(slbndtime-rtsl);
		}
	}
	
			

	if (totaltime>=windtime)
	{
		windthold=windthnew;
		windvold=windvnew;
		rtwind=windtime;
		fscanf(fwind,"%f\t%f\t%f",&windtime,&windvnew,&windthnew);
		//windtime=windtime+rtwind;
		//printf("windthold=%f\n",windthold);
		//printf("windthnew=%f\n",windthnew);
	}
		windth=windthold+(totaltime-rtwind)*(windthnew-windthold)/(windtime-rtwind);
		windv=windvold+(totaltime-rtwind)*(windvnew-windvold)/(windtime-rtwind);
		//printf("windv=%f\n",windv);

		windth=(1.5*pi-grdalpha)-windth*pi/180;
		//printf("windv=%f\twindth=%f\n",windv,windth);


	
	
	
}


void wavestep(void)
{

//Subroutine runs the wave model

	dim3 blockDim(16, 16, 1);
	dim3 gridDim(nx / blockDim.x, ny / blockDim.y, 1);

	dim3 blockDim4(4, 4, 1);
	dim3 gridDim4(nx / blockDim4.x, ny / blockDim4.y, 1);
	
	
	CUDA_SAFE_CALL( cudaMemcpy(St_g, St, ny*ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	//offshorebndWav(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g)
	offshorebndWav<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g);
	CUT_CHECK_ERROR("Offshore Wave bnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
//Sanity check
	sanity<<<gridDim, blockDim, 0>>>(nx, ny,eps,hh_g,sigm_g,ntheta,ee_g);
	
	CUT_CHECK_ERROR("sanity execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


	//CUDA_SAFE_CALL( cudaMalloc((void **)&cg_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&cx_g, nx*ny*ntheta*sizeof(float )) );
//	CUDA_SAFE_CALL( cudaMalloc((void **)&c_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&cy_g, nx*ny*ntheta*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&k_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&kh_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&sinh2kh_g, nx*ny*sizeof(float )) );

	
	
//dispersion
	dispersion<<<gridDim, blockDim, 0>>>(nx,ny,twopi,g,aphi,bphi,sigm_g,hh_g,k_g,c_g,kh_g,sinh2kh_g,cg_g);
	CUT_CHECK_ERROR("dispersion execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    
    //CUDA_SAFE_CALL( cudaMemcpy(C,kh_g,  ny*nx*sizeof(float ), cudaMemcpyDeviceToHost) );

	CUDA_SAFE_CALL( cudaMalloc((void **)&dhdx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dhdy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dudx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dudy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dvdx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dvdy_g, nx*ny*sizeof(float )) );
	
	
// Wave current interaction	(i.e remove wci in shallow water)
	calcwci<<<gridDim, blockDim, 0>>>(nx,ny,wci,hwci,hh_g,wci_g);
	CUT_CHECK_ERROR("calcwci execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


// // Slopes of water depth and velocities
    slopes<<<gridDim, blockDim, 0>>>(nx,ny,dx,hh_g,uu_g,vv_g,dhdx_g,dhdy_g,dudx_g,dudy_g,dvdx_g,dvdy_g);//
	CUT_CHECK_ERROR("slopes execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//CUDA_SAFE_CALL( cudaMalloc((void **)&cgx_g, nx*ny*ntheta*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&cgy_g, nx*ny*ntheta*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&ctheta_g, nx*ny*ntheta*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMemcpy(C,kh_g,  ny*nx*sizeof(float ), cudaMemcpyDeviceToHost) );
//Propagation speed in theta space
	propagtheta<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,wci_g,ctheta_g,/*c_g,cx_g,cy_g,*/cxsth_g,sxnth_g,/*uu_g,vv_g,*/dhdx_g,dhdy_g,dudx_g,dudy_g,dvdx_g,dvdy_g,sigm_g,kh_g);//
	CUT_CHECK_ERROR("propagtheta execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    
    
	//////////
	//CUDA_SAFE_CALL( cudaMemcpy(ctheta,ctheta_g,  ny*nx*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );
	//////////
	
	CUDA_SAFE_CALL( cudaFree(dhdx_g));
	CUDA_SAFE_CALL( cudaFree(dhdy_g) );
	CUDA_SAFE_CALL( cudaFree(dudx_g) );
	CUDA_SAFE_CALL( cudaFree(dudy_g) );
	CUDA_SAFE_CALL( cudaFree(dvdx_g) );
	CUDA_SAFE_CALL( cudaFree(dvdy_g) );


	//


	//read3Dnc(nx,ny,ntheta,"eeX.nc",ee);
	//CUDA_SAFE_CALL( cudaMemcpy(ee_g, ee, nx*ny*ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	
	


//
// transform to wave action
//
	action<<<gridDim, blockDim, 0>>>(ntheta,nx,ny,ee_g,sigm_g);
	CUT_CHECK_ERROR("action execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	


//
// Upwind Euler timestep propagation
//
	CUDA_SAFE_CALL( cudaMalloc((void **)&xadvec_g, nx*ny*ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&yadvec_g, nx*ny*ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&thetaadvec_g, nx*ny*ntheta*sizeof(float )) );

	xadvecupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,dt,wci_g,ee_g,cg_g,cxsth_g,uu_g,xadvec_g);
	CUT_CHECK_ERROR("eulerupwind xadvec execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
   

	yadvecupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,dt,wci_g,ee_g,cg_g,sxnth_g,vv_g,yadvec_g);
	CUT_CHECK_ERROR("eulerupwind yadvec execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

		
	//CUDA_SAFE_CALL( cudaMalloc((void **)&eect_g, nx*ny*ntheta*sizeof(float )) );
	
	//eectheta<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,ee_g,ctheta_g,eect_g);
	//CUT_CHECK_ERROR("eulerupwind eectheta execution failed\n");
    //CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    //thetaadvecuw<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,eect_g,thetaadvec_g);
    //CUT_CHECK_ERROR("eulerupwind thetaadvecuw execution failed\n");
    //CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
	thetaadvecupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,wdt,wci,ee_g,ctheta_g,thetaadvec_g);
	CUT_CHECK_ERROR("eulerupwind thetaadvec execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
     //CUDA_SAFE_CALL( cudaMemcpy(ctheta,yadvec_g,  ny*nx*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );
    
    //CUDA_SAFE_CALL( cudaMemcpy(ctheta,thetaadvec_g,  ny*nx*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );
    
    
    //read3Dnc(nx,ny,ntheta,"xadvecX.nc",ee);
	//CUDA_SAFE_CALL( cudaMemcpy(xadvec_g, ee, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read3Dnc(nx,ny,ntheta,"yadvecX.nc",ee);
	//CUDA_SAFE_CALL( cudaMemcpy(yadvec_g, ee, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read3Dnc(nx,ny,ntheta,"thetaadvecX.nc",ee);
	//CUDA_SAFE_CALL( cudaMemcpy(thetaadvec_g, ee, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
    
    

	eulerupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,dt,wci,ee_g,xadvec_g,yadvec_g,thetaadvec_g);
	CUT_CHECK_ERROR("eulerupwind  execution failed\n");
        CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	
	



	//CUDA_SAFE_CALL( cudaFree(cgx_g));
	//CUDA_SAFE_CALL( cudaFree(cgy_g));
	//Fix lateraL BND
	rollerlatbnd<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,eps,hh_g,ee_g);
	CUT_CHECK_ERROR("energy latbnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );



//
// transform back to wave energy
//
	energy<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,ee_g,sigm_g);
	CUT_CHECK_ERROR("energy execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
    //CUDA_SAFE_CALL( cudaMemcpy(ctheta,ee_g,  ny*nx*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );

	//CUDA_SAFE_CALL( cudaMalloc((void **)&H_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&E_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&D_g, nx*ny*sizeof(float )) );
	
	//CUDA_SAFE_CALL( cudaMemcpy(ctheta,ee_g,  ny*nx*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );

//
// Energy integrated over wave directions,Hrms
//
	energint<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,rho,g,gammax,E_g,H_g,hh_g,ee_g);
	CUT_CHECK_ERROR("energint execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	
	

//
// calculate change in intrinsic frequency // removed because it is super slow and doesn't do much
//
// tm is thetamean and it is calculated in the mean dir scheme
//	CUDA_SAFE_CALL( cudaMalloc((void **)&tm_g, nx*ny*sizeof(float )) );
//	calctm<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,tm_g,theta_g,ee_g);
//	CUT_CHECK_ERROR("energint execution failed\n");
//    CUDA_SAFE_CALL( cudaThreadSynchronize() );
/*
//Change of intrinsec frequency
*/





// 
//  Total dissipation from breaking  and bottom friction
//

	if (breakmod==1)
	{
		roelvink<<<gridDim, blockDim, 0>>>(nx,ny,rho,g,gammaa,alpha,n,Trep,fwm_g,cfm_g,hh_g,H_g,E_g,D_g,k_g);
		CUT_CHECK_ERROR("roelvink execution failed\n");
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
	}
	else
	{
		baldock<<<gridDim, blockDim, 0>>> (nx,ny,rho,g,gammaa,alpha,n,Trep,fwm_g,cfm_g,hh_g,H_g,E_g,D_g,k_g);//Baldock more appropriate for pseudo stationary cases
		CUT_CHECK_ERROR("baldoc execution failed\n");
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
	}
//
//  Calculate roller energy balance
//
//CUDA_SAFE_CALL( cudaMemcpy(hhmean,E_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );

if (roller==1)
{
	xadvecupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,dt,wci_g,rr_g,c_g,cxsth_g,uu_g,xadvec_g);
	CUT_CHECK_ERROR("eulerupwind xadvec execution failed\n");
        CUDA_SAFE_CALL( cudaThreadSynchronize() );

	yadvecupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,dt,wci_g,rr_g,c_g,sxnth_g,vv_g,yadvec_g);
	CUT_CHECK_ERROR("eulerupwind yadvec execution failed\n");
        CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//eectheta<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,rr_g,ctheta_g,eect_g);
	//CUT_CHECK_ERROR("eulerupwind eectheta execution failed\n");
    //CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    //thetaadvecuw<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,eect_g,thetaadvec_g);
    //CUT_CHECK_ERROR("eulerupwind thetaadvecuw execution failed\n");
    //CUDA_SAFE_CALL( cudaThreadSynchronize() );	

	thetaadvecupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,wdt,wci,rr_g,ctheta_g,thetaadvec_g);
	CUT_CHECK_ERROR("eulerupwind thetaadvec execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	eulerupwind<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,dx,dt,wci,rr_g,xadvec_g,yadvec_g,thetaadvec_g);
	CUT_CHECK_ERROR("eulerupwind  execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

//
//  Adjust lateral bnds
//

	rollerlatbnd<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,eps,hh_g,rr_g);
	CUT_CHECK_ERROR("rollerlatbnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
}
	//CUDA_SAFE_CALL( cudaFree(eect_g));
	CUDA_SAFE_CALL( cudaFree(xadvec_g));
	CUDA_SAFE_CALL( cudaFree(yadvec_g));
	CUDA_SAFE_CALL( cudaFree(thetaadvec_g));

//read2Dnc(nx,ny,"D.nc",uu);
//CUDA_SAFE_CALL( cudaMemcpy(D_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );


// 
//  Distribution of dissipation over directions and frequencies
//                               
	dissipation<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,eps,dt,g,beta,wci_g,hh_g,ee_g,D_g,E_g,rr_g,c_g,cxsth_g,sxnth_g,uu_g,vv_g,DR_g,R_g);
	CUT_CHECK_ERROR("dissipation execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );



//
//Fix lateraL BND
//
	rollerlatbnd<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,eps,hh_g,ee_g);
	CUT_CHECK_ERROR("energy latbnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	

	
// 
//  Compute mean wave direction
// 

	meandir<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,rho,g,dtheta,ee_g, theta_g,thetamean_g,E_g,H_g);
	CUT_CHECK_ERROR("meandir execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );



//
//	Constant warm start // WARNING ONLY TO BE USED FOR DEBUGGING
//
//	read3Dnc(nx,ny,ntheta,"eeX.nc",ee);
//	CUDA_SAFE_CALL( cudaMemcpy(ee_g, ee, nx*ny*ntheta*sizeof(float ), cudaMemcpyHostToDevice) );

	
	
	

//
// Radiation stresses and forcing terms
//
	
	CUDA_SAFE_CALL( cudaMalloc((void **)&Sxx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Sxy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Syy_g, nx*ny*sizeof(float )) );

	radstress<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dx,dtheta,ee_g,rr_g,cxsth_g,sxnth_g,cg_g,c_g,Sxx_g,Sxy_g,Syy_g);

	CUT_CHECK_ERROR("radstress execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

//	
// Wave forces
//
	wavforce<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dx,dtheta,Sxx_g,Sxy_g,Syy_g,Fx_g,Fy_g,hh_g);
	CUT_CHECK_ERROR("wavforce execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	twodimbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,Fx_g);
	CUT_CHECK_ERROR("wave force X bnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	twodimbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,Fy_g);
	CUT_CHECK_ERROR("wave force Y bnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//CUDA_SAFE_CALL( cudaMemcpy(ctheta,ctheta_g,  ny*nx*ntheta*sizeof(float ), cudaMemcpyDeviceToHost) );
	

//
// CAlculate stokes velocity and breaker delay //Breaker delay removed because it is slow and kinda useless
//
	breakerdelay<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,g,rho,Trep,eps,urms_g,ust_g,H_g,E_g,c_g,k_g,hh_g,R_g);
	CUT_CHECK_ERROR("breakerdelay execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	twodimbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,urms_g);
	CUT_CHECK_ERROR("wave force Y bnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	twodimbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,ust_g);
	CUT_CHECK_ERROR("wave force Y bnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


	CUDA_SAFE_CALL( cudaFree(Sxy_g) );
	CUDA_SAFE_CALL( cudaFree(Sxx_g) );
	CUDA_SAFE_CALL( cudaFree(Syy_g) );
	//CUDA_SAFE_CALL( cudaFree(cg_g));
	//CUDA_SAFE_CALL( cudaFree(c_g));
	CUDA_SAFE_CALL( cudaFree(tm_g));




//
// Adjust Offshore Bnd
//

//CUDA_SAFE_CALL( cudaMemcpy(St_g, St, ny*ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
//offshorebndWav(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g)
//offshorebndWav<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g);
//CUT_CHECK_ERROR("Offshore Wave bnd execution failed\n");
//CUDA_SAFE_CALL( cudaThreadSynchronize() );



CUDA_SAFE_CALL( cudaFree(E_g));
//CUDA_SAFE_CALL( cudaFree(H_g));
//CUDA_SAFE_CALL( cudaFree(D_g));

//CUDA_SAFE_CALL( cudaFree(k_g));
CUDA_SAFE_CALL( cudaFree(kh_g));
CUDA_SAFE_CALL( cudaFree(sinh2kh_g));




	



}
void flowstep(void)
{
// Flow model timestep
	dim3 blockDim(16, 16, 1);
	dim3 gridDim(nx / blockDim.x, ny / blockDim.y, 1);

	dim3 blockDim4(4, 4, 1);
	dim3 gridDim4(nx / blockDim4.x, ny / blockDim4.y, 1);
	
	
	//////////////////////////////////////////
	// BELOW IS FOR DEBUGGING ONLY
	/////////////////////////////////////////
	//read2Dnc(nx,ny,"Hfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(H_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"Fxfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(Fx_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"Fyfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(Fy_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"urmsfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(urms_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"ustfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(ust_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"thetameanfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(thetamean_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"uufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(uu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vvfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vv_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"zsfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(zs_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hhfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hh_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"ueufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(ueu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vevfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vev_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	
	
	//read2Dnc(nx,ny,"hufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"humfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hum_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hvfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hv_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hvmfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hvm_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vmageufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vmageu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vmagevfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vmagev_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );


	// Set weakly reflective offshore boundary
	ubnd<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,g,rho,totaltime,wavbndtime,dtwavbnd,slbndtime,rtsl,zsbndold,zsbndnew,Trep,qbndold_g,qbndnew_g,zs_g,uu_g,vv_g,vu_g,umeanbnd_g,vmeanbnd_g,zb_g,cg_g,hum_g,cfm_g,Fx_g,hh_g);
	CUT_CHECK_ERROR("ubnd execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


	//
	// Water level slopes
	//
	wlevslopes<<<gridDim, blockDim, 0>>>(nx,ny,dx,eps,zs_g,dzsdx_g,dzsdy_g,hh_g);
	CUT_CHECK_ERROR("wlevslopes execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	
	 
	//
	// Water depth at u pts for momentum and continuity eq (hum hu)
	//
	
	udepthmomcont<<<gridDim, blockDim, 0>>>(nx,ny,dx,eps,uumin,wetu_g,zs_g,uu_g,hh_g,hum_g,hu_g,zb_g);
	CUT_CHECK_ERROR("udepthmomcont execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//
	// Water depth at v pts for momentum and continuity eq (hvm hv)
	//
	
	vdepthmomcont<<<gridDim, blockDim, 0>>>(nx,ny,dx,eps,uumin,wetv_g,zs_g,vv_g,hh_g,hvm_g,hv_g,zb_g);
	CUT_CHECK_ERROR("vdepthmomcont execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	
	
	
	
	
	//
	// Advection in the x direction using 2n order finite difference
	//
	
	ududx_adv2<<<gridDim, blockDim, 0>>>(nx,ny,dx,hu_g,hum_g,uu_g,ududx_g);
	CUT_CHECK_ERROR("uadvec execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    

	//vdudy
	vdudy_adv2<<<gridDim, blockDim, 0>>>(nx,ny,dx,hv_g,hum_g,uu_g,vv_g,vdudy_g);
	CUT_CHECK_ERROR("uadvec execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    


	//
	// Smagorinsky formulation or Normal eddy viscosity
	//
	CUDA_SAFE_CALL( cudaMalloc((void **)&nuh_g, nx*ny*sizeof(float )) );
	smago<<<gridDim, blockDim, 0>>>(nx,ny,dx, uu_g,vv_g,nuh, nuh_g,usesmago);
	CUT_CHECK_ERROR("uadvec execution failed\n");
    	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//
	// increase eddy viscosity by wave induced breaking as in Reniers 2004 & Set viscu = 0.0 near water line
	//
	CUDA_SAFE_CALL( cudaMalloc((void **)&viscu_g, nx*ny*sizeof(float )) );
	viscou<<<gridDim, blockDim, 0>>>(nx,ny,dx,rho,eps,nuhfac,nuh_g,hh_g,hum_g,hvm_g,DR_g,uu_g,wetu_g,viscu_g);
	CUT_CHECK_ERROR("visco execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


    //
    // Explicit Euler step momentum u-direction
    //
      
	eulerustep<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,g,rho,cfm_g,fc,windth,windv,Cd,uu_g,urms_g,ududx_g,vdudy_g,viscu_g,dzsdx_g,hu_g,hum_g,Fx_g,vu_g,ueu_g,vmageu_g,wetu_g);
	CUT_CHECK_ERROR("eulerustep execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    
    
    

	//
	// Advection in the y direction using 2n order finite difference
	//
	//vdvdy
	vdvdy_adv2<<<gridDim, blockDim, 0>>>(nx,ny,dx,hv_g,hvm_g,vv_g,vdvdy_g);
	CUT_CHECK_ERROR("vadvec for v execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	//udvdx
	
	udvdx_adv2<<<gridDim, blockDim, 0>>>(nx,ny,dx,hu_g,hvm_g,uu_g,vv_g,udvdx_g);
	CUT_CHECK_ERROR("vadvec for v execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//
	// increase eddy viscosity by wave induced breaking as in Reniers 2004 & Set viscv = 0.0 near water line
	//
	CUDA_SAFE_CALL( cudaMalloc((void **)&viscv_g, nx*ny*sizeof(float )) );
	viscov<<<gridDim, blockDim, 0>>>(nx,ny,dx,rho,eps,nuhfac,nuh_g,hh_g,hum_g,hvm_g,DR_g,vv_g,wetv_g,viscv_g);
	CUT_CHECK_ERROR("visco v execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	CUDA_SAFE_CALL( cudaFree(nuh_g));
	
	
	viscovbnd<<<gridDim, blockDim, 0>>>(nx,ny,viscv_g );
	CUT_CHECK_ERROR("visco v execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	

	//
    // Explicit Euler step momentum v-direction
	//
	
	eulervstep<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,g,rho,cfm_g,fc,windth,windv,Cd,vv_g,urms_g,udvdx_g,vdvdy_g,viscv_g,dzsdy_g,hv_g,hvm_g,Fy_g,uv_g,vev_g,vmagev_g,wetv_g);
	CUT_CHECK_ERROR("eulervstep execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
	//
	// Adjust lateral bnds
	//
	uuvvzslatbnd<<<gridDim, blockDim, 0>>>(nx,ny,uu_g,vv_g,zs_g);
	CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    
	
	//
	//v velocities at u pts and u velocities at v pts
	//
	
	calcuvvu<<<gridDim, blockDim, 0>>>(nx,ny,dx,uu_g,vv_g,vu_g,uv_g,ust_g,thetamean_g,ueu_g,vev_g,vmageu_g,vmagev_g,wetu_g,wetv_g);
	CUT_CHECK_ERROR("calcuvvu execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


	//
	//Calculate hu
	//
	depthhu<<<gridDim, blockDim, 0>>>(nx,ny,dx,uumin,eps,hh_g,uu_g,hu_g,zs_g,zb_g);
	CUT_CHECK_ERROR("depthhu execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//
	//Calculate hv
	//
	depthhv<<<gridDim, blockDim, 0>>>(nx,ny,dx,uumin,eps,hh_g,vv_g,hv_g,zs_g,zb_g);
	CUT_CHECK_ERROR("depthhv execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );


	//
	// Update water level using continuity eq.
	//
	continuity<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,eps,uu_g,hu_g,vv_g,hv_g,zs_g,hh_g,zb_g,dzsdt_g);
	CUT_CHECK_ERROR("continuity execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	
		//
	// Adjust lateral bnds
	//
	uuvvzslatbnd<<<gridDim, blockDim, 0>>>(nx,ny,uu_g,vv_g,zs_g);
	CUT_CHECK_ERROR("uu vv zs lateral bnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );


	hsbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,zb_g,zs_g);
	CUT_CHECK_ERROR("hh lateral bnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	


	CUDA_SAFE_CALL( cudaFree(viscu_g));
	CUDA_SAFE_CALL( cudaFree(viscv_g));


}

void sedimentstep(void)
{
// suspended sediment timestep
	dim3 blockDim(16,16, 1);
	dim3 gridDim(nx / blockDim.x, ny / blockDim.y, 1);

	dim3 blockDim4(4, 4, 1);
	dim3 gridDim4(nx / blockDim4.x, ny / blockDim4.y, 1);
	
	/////////////////////////////////////////////////////
	// BELOW IS FOR DEBUGGING ONLY
	/////////////////////////////////////////////////////
	
	//read2Dnc(nx,ny,"uufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(uu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vvfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vv_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hvfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hv_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	
	//read2Dnc(nx,ny,"ueufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(ueu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	
	//read2Dnc(nx,ny,"vevfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vev_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"Hfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(H_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hhfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hh_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"urmsfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(urms_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"uvfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(uv_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	
//
// Compute long wave turbulence due to breaking
//
	longturb<<<gridDim, blockDim, 0>>>(nx,ny,dx,rho,g,dt,beta,c_g,kturb_g,rolthick_g,dzsdt_g,uu_g,vv_g,hu_g,hv_g,wetu_g,wetv_g,hh_g);
	CUT_CHECK_ERROR("longturb execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
//
// Calculate Equilibrium concentration Ceq
//
	//CUDA_SAFE_CALL( cudaMalloc((void **)&ceqsg_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&ceqbg_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Tsg_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&ua_g, nx*ny*sizeof(float )) );
	//BEWARE BELOW SHOULD BE hh_old_g
	//Sbvr or Sednew
	Sbvr<<<gridDim, blockDim, 0>>>(nx,ny,rho,g,eps,Trep,D50,D90,rhosed,wws,nuhfac,ueu_g,vev_g,H_g,DR_g,R_g,c_g,hh_g,urms_g,ceqsg_g,ceqbg_g,Tsg_g,cfm_g,kturb_g);
	CUT_CHECK_ERROR("CalcCeq execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );

	Rvr<<<gridDim, blockDim, 0>>>(nx,ny,Trep,facsk,facas,H_g, hh_g,urms_g, c_g, ua_g);
	CUT_CHECK_ERROR("Rvr execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );	
    
   // CUDA_SAFE_CALL( cudaMemcpy(uumean, ua_g, nx*ny*sizeof(float ), cudaMemcpyDeviceToHost) );

	
	
	
//
// Limit erosion to available sediment on top of hard layer
//
	CUDA_SAFE_CALL( cudaMalloc((void **)&facero_g, nx*ny*sizeof(float )) );
	Erosus<<<gridDim, blockDim, 0>>>(nx,ny,dt,morfac,por,hh_g,ceqsg_g,ceqbg_g,Tsg_g,facero_g,stdep_g);
	CUT_CHECK_ERROR("Erosus execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );



//////////////////
// LOOP for sediment fractions should come here
/////////////////////

//
// suspended load 
// 
	CUDA_SAFE_CALL( cudaMalloc((void **)&Sus_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Svs_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Sub_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Svb_g, nx*ny*sizeof(float )) );
	Susp<<<gridDim, blockDim, 0>>>(nx,ny,dx,eps,nuh,nuhfac,rho,sus,bed,ueu_g,vev_g,uu_g,uv_g,hu_g,vv_g,vu_g,hv_g,zb_g,hh_g,DR_g, Cc_g,ceqbg_g,Sus_g, Svs_g,Sub_g,Svb_g,thetamean_g,ua_g);
	CUT_CHECK_ERROR("Susp execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
 
//
//Calculate suspended concentration
//

	CUDA_SAFE_CALL( cudaMalloc((void **)&ero_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&depo_g, nx*ny*sizeof(float )) );
	Conc<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,eps,hh_g,Cc_g,ceqsg_g,Tsg_g,facero_g,ero_g,depo_g,Sus_g,Svs_g);
	CUT_CHECK_ERROR("Conc execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    

//
// Update global variables and fix bnds
//
	CClatbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,Cc_g);
	CUT_CHECK_ERROR("CClatbnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	
if (morfac>0.0f)// Only if morphology is need i.e. if imodel=4 and morphac >0.0
	{	
//
// Adjust sediment fluxes for rocklayer
//

	CUDA_SAFE_CALL( cudaMalloc((void **)&Sout_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&indSub_g, nx*ny*sizeof(int )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&indSvb_g, nx*ny*sizeof(int )) );
	hardlayer<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,Sub_g,Svb_g,Sout_g, indSub_g,indSvb_g);
	CUT_CHECK_ERROR("hardlayer execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,Sout_g);
	CUT_CHECK_ERROR("Sout twodimbnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
		
	
//
// Bed update
//
	bedupdate<<<gridDim, blockDim, 0>>>(nx,ny,eps,dx,dt,morfac,por,hh_g,ero_g,depo_g,Sub_g,Svb_g,Sout_g,indSub_g,indSvb_g,zb_g,dzb_g,stdep_g);
	CUT_CHECK_ERROR("bedupdate execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	
//
// Update lateral bnd	
//	
	zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,zb_g);
	CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	
	zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,stdep_g);
	CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	
//
// Avalanching
//
	CUDA_SAFE_CALL( cudaMalloc((void **)&ddzb_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMemcpy(ddzb_g,zeros, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	avalanching<<<gridDim, blockDim, 0>>>(nx,ny,eps,dx,dt,por,drydzmax,wetdzmax,maxslpchg,hh_g,zb_g,ddzb_g,stdep_g);
	CUT_CHECK_ERROR("avalanching execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );

//
// Update Zb for avalanching
//

	updatezb<<<gridDim, blockDim, 0>>>(nx,ny,dx,dt,zb_g,ddzb_g,dzb_g,zs_g,hh_g,stdep_g);
	CUT_CHECK_ERROR("avalanching execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
//
// Update lateral bnd	
//	
	zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,zb_g);
	CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,dzb_g);
	//CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
	//CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	zblatbnd<<<gridDim, blockDim, 0>>>(nx,ny,stdep_g);
	CUT_CHECK_ERROR("Zb twodimbnd execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	updatezom<<<gridDim, blockDim, 0>>>(nx,ny,cf,cf2,fw,fw2,stdep_g,cfm_g,fwm_g);
	CUT_CHECK_ERROR("UpdateZom execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
}


	
	
///////////////////
// END LOOP sediment fraction

	CUDA_SAFE_CALL( cudaFree(ddzb_g));
	CUDA_SAFE_CALL( cudaFree(Sout_g));
	CUDA_SAFE_CALL( cudaFree(indSub_g));
	CUDA_SAFE_CALL( cudaFree(indSvb_g));
	CUDA_SAFE_CALL( cudaFree(ua_g));
	
	
	CUDA_SAFE_CALL( cudaFree(ero_g));
	CUDA_SAFE_CALL( cudaFree(depo_g));
	CUDA_SAFE_CALL( cudaFree(Sus_g));
	CUDA_SAFE_CALL( cudaFree(Svs_g));
	CUDA_SAFE_CALL( cudaFree(Sub_g));
	CUDA_SAFE_CALL( cudaFree(Svb_g));
	CUDA_SAFE_CALL( cudaFree(facero_g));
	//CUDA_SAFE_CALL( cudaFree(ceqsg_g));
	CUDA_SAFE_CALL( cudaFree(ceqbg_g));
	CUDA_SAFE_CALL( cudaFree(Tsg_g));



}



int main(int argc, char **argv)
{
	
	// Start timer to keep track of time 
	clock_t starttime, endtime;
	
	
	starttime=clock();

	//////////////////////////////////////////////////////
    /////             Read Operational file          /////
    //////////////////////////////////////////////////////
 
 
    char opfile[]="opfile.dat"; // Compulsory input file
 
    char filename[256];
	
	char slbnd[256];
	char windfile[256];
	char zofile[256];
	char HLfile[256];
	
	
      
     
        FILE * fop;
        fop=fopen(opfile,"r");
        fscanf(fop,"%*s");//Dummy string
        fscanf(fop,"%s\t%*s",&filename);// Bathy file name needs to be md format
		fscanf(fop,"%d\t%*s",&imodel);// Type of model: 1: wave only; 2: currents only 3: waves+currents 4:waves+currents+sediment(+ morphology if morfac>0)
		fscanf(fop,"%d\t%*s",&GPUDEVICE);// What GPU device to use 
		fscanf(fop,"%f\t%*s",&dt);// Model time step in s. //This should be calculated by the model
		fscanf(fop,"%d\t%*s",&nstpw);// Number of flow and sediment step between wave step needs to be 1 for unsteady runs 
		fscanf(fop,"%f\t%*s",&eps);//drying height in m
		fscanf(fop,"%f,%f\t%*s",&cf,&cf2);// bottom friction for flow model cf is for sand and cf2 fro reef area (Reef and sand discrimination is done based on structure file if none is present cf2 should not be used )
		//fscanf(fop,"%f\t%*s",&cf);
		//fscanf(fop,"%s\t%*s",&zofile);
		fscanf(fop,"%f\t%*s",&nuh); // Viscosity coeff or samgo coeff depennding on usesmago
		fscanf(fop,"%f\t%*s",&nuhfac);//nuhfac=1.0f;//0.001f; //viscosity coefficient for roller induced turbulent horizontal viscosity// it should be small contrary to what XBeach recommend as default
		fscanf(fop,"%d\t%*s",&usesmago);// Uses smagorynsky formulation to calculate viscosity 0: No 1: Yes
		fscanf(fop,"%f\t%*s",&lat);// Latitude of the grid use negative for south hemisphere (this implies the grid is small on earth scale)
		fscanf(fop,"%f\t%*s",&Cd);// Wind drag coeff
		fscanf(fop,"%f\t%*s",&wci); // Wave current interaction switch (can also be used as a number between 0 and 1 to reduce the interaction if unstable)
		fscanf(fop,"%f\t%*s",&hwci);// hwci=0.010f;//min depth for wci
		fscanf(fop,"%d\t%*s",&breakmod); // Wave dissipation model 1: roelvink 2: Baldock. use 1 for unsteady runs (i.e. with wave group) and use 2 for steady runs
		fscanf(fop,"%f\t%*s",&gammaa);// Wave breaking gamma param 
		fscanf(fop,"%f\t%*s",&n);// exponential; in Roelving breaking model
		fscanf(fop,"%f\t%*s",&alpha);// calibration for wave dissipation (should be 1)
		fscanf(fop,"%f\t%*s",&gammax);//gammax=2.0f; //maximum ratio Hrms/hh
		fscanf(fop,"%f\t%*s",&beta);// Roller slope dissipation param
		fscanf(fop,"%f,%f\t%*s",&fw,&fw2);// Wave bottom dissipation parameters fw is for sand fw2 is for reefs. see cf comments
		//fscanf(fop,"%f\t%*s",&fw);
		fscanf(fop,"%f,%f\t%*s",&D50,&D90);// sand grain size in m
		//fscanf(fop,"%f\t%*s",&D50);
		//fscanf(fop,"%f\t%*s",&D90);
		fscanf(fop,"%f\t%*s",&rhosed);// sand density
		fscanf(fop,"%f\t%*s",&wws);// sand fall velocity (should be calculated)
		fscanf(fop,"%f,%f\t%*s",&drydzmax,&wetdzmax);// max slope in avalannching model
		//fscanf(fop,"%f\t%*s",&drydzmax);
		//fscanf(fop,"%f\t%*s",&wetdzmax);
		fscanf(fop,"%f\t%*s",&maxslpchg);// max change within a step to avoid avalanching tsunami
		fscanf(fop,"%f\t%*s",&por);// sand porosity (should not be constant)
		fscanf(fop,"%f\t%*s",&morfac);// morphological factor 0 no changes in morphology 1 normal changes in morpho >1 accelerated morphological changes (beware this doesn't accelerate the bnd you have to do this manually)
		fscanf(fop,"%f,%f\t%*s",&sus,&bed);// calibration coeff for suspended load and bed load
		//fscanf(fop,"%f\t%*s",&sus);
		//fscanf(fop,"%f\t%*s",&bed);
		fscanf(fop,"%f,%f\t%*s",&facsk,&facas);// calibration factor for wave skewness and Asymetry
		//fscanf(fop,"%f\t%*s",&facsk);
		//fscanf(fop,"%f\t%*s",&facas);
		fscanf(fop,"%s\t%*s",&HLfile);// Structure file write down "none" if none present
		fscanf(fop,"%d\t%*s",&wavebndtype); // 1 is quasistationary wave spectrum; 2 is for infrgravity and long bound waves Xbeach type
		fscanf(fop,"%s\t%*s",&wavebndfile);// wave bnd file see wiki for details
		fscanf(fop,"%s\t%*s",&slbnd); // tide/surge bnd file
		fscanf(fop,"%s\t%*s",&windfile); // Wind forcing file
		//fscanf(fop,"%d\t%*s",&npart);
		fscanf(fop,"%d\t%*s",&sedstart);// which step to start sediment transport and morpho
		//fscanf(fop,"%f\t%*s",&Hplotmax);
		//fscanf(fop,"%d\t%*s",&nstepplot);
		fscanf(fop,"%d\t%*s",&nstepout); // output step
		fscanf(fop,"%d\t%*s",&endstep);// end step
		//fscanf(fop,"%d\t%d\t%*s",&iout,&jout);
		fscanf(fop,"%s\t%*s",&tsoutfile);// output file
		fclose(fop);
       
        printf("bathy file: %s\n",filename);
		printf("Imodel: %d\n",imodel);
		//printf("nstepplot: %d\n",nstepplot);
		printf("smago?: %d\n",usesmago);
		printf("bed: %f\n",bed);
		
		printf("facsk=%f facas=%f\n",facsk,facsk);


        
	wdt=0.0;

	FILE * fid;
	
 
    //read input data:
	printf("bathy: %s\n",filename);
 
    	fid=fopen(filename,"r");
	fscanf(fid,"%u\t%u\t%f\t%*f\t%f",&nx,&ny,&dx,&grdalpha);
	printf("nx=%d\tny=%d\tdx=%f\talpha=%f\n",nx,ny,dx,grdalpha);

	grdalpha=grdalpha*pi/180; // grid rotation

	if(imodel>=2)
	{
	printf("Opening sea level bnd...");
	fsl=fopen(slbnd,"r");

	fscanf(fsl,"%f\t%f",&rtsl,&zsbndold);
	
	//Note: the first rtsl should be 0 
	fscanf(fsl,"%f\t%f",&slbndtime,&zsbndnew);
	printf("done\n");

	//zsbnd sea level in bnd file
	//rtsl bnd time
	//slbndtime=rtsl;
	}
	else
	{
		zsbndold=0;
		zsbndnew=0;
	}
	
	// Allocate CPU memory
	
	hh=(float *)malloc(nx*ny*sizeof(float));
	uu=(float *)malloc(nx*ny*sizeof(float));
	vv=(float *)malloc(nx*ny*sizeof(float));
	zs=(float *)malloc(nx*ny*sizeof(float));
	zb=(float *)malloc(nx*ny*sizeof(float));
	cfm=(float *)malloc(nx*ny*sizeof(float));
	dzb=(float *)malloc(nx*ny*sizeof(float));
	stdep=(float *)malloc(nx*ny*sizeof(float));
	zeros=(float *)malloc(nx*ny*sizeof(float));
	umeanbnd=(float *)malloc(ny*sizeof(float));

	

	// set initital condition and read bathy file
	printf("Set initial condition\n");
	
	int jread;
    for (int fnod=ny; fnod>=1;fnod--)
    {
		
		fscanf(fid,"%u",&jread);
		umeanbnd[(jread-1)]=0.0f;
        for(int inod=0; inod<nx; inod++)
        {
			fscanf(fid,"%f",&zb[inod+(jread-1)*nx]);
			uu[inod+(jread-1)*nx]=0.0f;
			vv[inod+(jread-1)*nx]=0.0f;
			dzb[inod+(jread-1)*nx]=0.0f;
			cfm[inod+(jread-1)*nx]=cf;
			stdep[inod+(jread-1)*nx]=0.0f;
			zeros[inod+(jread-1)*nx]=0.0f;
			zs[inod+(jread-1)*nx]=max(zsbndold,-1*zb[inod+(jread-1)*nx]);
			hh[inod+(jread-1)*nx]=max(zb[inod+(jread-1)*nx]+zsbndold,eps);
			
			

		}
	}

	fclose(fid);
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
	

	int testfrictfile=strcmp(HLfile,nofrictionfile);
	if (testfrictfile!=0)
	{
	printf("Hard layer file found\n");	
	fid=fopen(HLfile,"r");
	int jx,jy;
	fscanf(fid,"%u\t%u\t%*f\t%*f\t%*f",&jx,&jy,&dx,&grdalpha);
	
	if (jx!=nx || jy!=ny)
	{
		printf("Error Hard layer file dimension mismatch. Model will run with constant friction.\n");
	}
	}
	else
	{
		printf("No hard layer file found\n");
	}

	for (int fnod=ny; fnod>=1;fnod--)
	{
	 	if(testfrictfile!=0)
	 	{fscanf(fid,"%u",&jread);}
		for(int inod=0; inod<nx; inod++)
        {
			if(testfrictfile!=0)
			{
				fscanf(fid,"%f",&stdep[inod+(jread-1)*nx]);
			}
			else
			{
				stdep[inod+(fnod-1)*nx]=5.0f;
				
			}

		}
	}
	if(testfrictfile!=0)
	{
	fclose(fid);
	}
	

	if(imodel==1 || imodel>=3)
	{
		waveinit();
	}
	
	
	// Read Wind forcing
	printf("Opening wind bnd\n");
	fwind=fopen(windfile,"r");
	fscanf(fwind,"%f\t%f\t%f",&rtwind,&windvold,&windthold);
	fscanf(fwind,"%f\t%f\t%f",&windtime,&windvnew,&windthnew);
	
	
	
	

	windv=windvold;
	windth=(1.5f*pi-grdalpha)-windthold*pi/180.0f;


	
	//zo=0.1;//roughness length
	//cf=zo;//zo;
	//lat=-35.0;
	//calculate coriolis force
	lat = lat*pi/180.0f;
	float wearth = pi*(1/24)/1800.0f;
	fc = 2.0f*wearth*sin(lat);
	

	//gammax=2.0f; //maximum ratio Hrms/hh
	//wci=1.0f; //switch for wave/current interaction.
	//hwci=0.010f;//min depth for wci
	//gammaa = 0.55; //breaker parameter in Baldock or Roelvink formulation
	//n=10.0;// power in roelvink dissipation model
	//alpha = 1.;//! wave dissipation coefficient
	//beta=0.15f;//! breaker slope coefficient in roller model
	roller= 1; // option to turn off/on roller model (0/1)

	//nuh=0.05; // Eddy viscosity [m2/s]
	//nuhfac=1.0f;//0.001f; //viscosity coefficient for roller induced turbulent horizontal viscosity
	
	t1=-(pi)/2;
	//thetamin=-60;
	//thetamax=60;
	
	//dtheta=10;
	
	//Allocate More array on CPU

	Fx=(float *)malloc(nx*ny*sizeof(float));
	Fy=(float *)malloc(nx*ny*sizeof(float));	
	zsbnd=(float *)malloc(ny*sizeof(float));

	cgx=(float *)malloc(nx*ny*ntheta*sizeof(float));
	cgy=(float *)malloc(nx*ny*ntheta*sizeof(float));
	cx=(float *)malloc(nx*ny*ntheta*sizeof(float));
	cy=(float *)malloc(nx*ny*ntheta*sizeof(float));
	ctheta=(float *)malloc(nx*ny*ntheta*sizeof(float));
	



	

	
	//Sxx=(float *)malloc(nx*ny*sizeof(float));
	//Syy=(float *)malloc(nx*ny*sizeof(float));
	//Sxy=(float *)malloc(nx*ny*sizeof(float));
	usd=(float *)malloc(nx*ny*sizeof(float));
	D=(float *)malloc(nx*ny*sizeof(float));
	E=(float *)malloc(nx*ny*sizeof(float));
	H=(float *)malloc(nx*ny*sizeof(float));
	urms=(float *)malloc(nx*ny*sizeof(float));
	ueu=(float *)malloc(nx*ny*sizeof(float));
	vev=(float *)malloc(nx*ny*sizeof(float));
	thetamean=(float *)malloc(nx*ny*sizeof(float));


	Hmean=(float *)malloc(nx*ny*sizeof(float));
	uumean=(float *)malloc(nx*ny*sizeof(float));
	vvmean=(float *)malloc(nx*ny*sizeof(float));
	hhmean=(float *)malloc(nx*ny*sizeof(float));
	zsmean=(float *)malloc(nx*ny*sizeof(float));
	Cmean=(float *)malloc(nx*ny*sizeof(float));
	
	
	
// Set initial water level on offshore bnd
	
for (int jj=0; jj<ny; jj++)
{
		zs[0+jj*nx]=zsbndold;
		hh[0+jj*nx]=zb[0+jj*nx]+zs[0+jj*nx];
	

}
	// Allocate more CPU memory
	

	omega = 2*pi/Trep;

	sigm= (float *)malloc(nx*ny*sizeof(float));
	//sigt= (float *)malloc(nx*ny*ntheta*sizeof(float));
	thet= (float *)malloc(nx*ny*ntheta*sizeof(float));
	//costhet=(float *)malloc(nx*ny*ntheta*sizeof(float));
	//sinthet=(float *)malloc(nx*ny*ntheta*sizeof(float));



	

	k= (float *)malloc(nx*ny*sizeof(float));
	c= (float *)malloc(nx*ny*sizeof(float));
	kh= (float *)malloc(nx*ny*sizeof(float));
	cg= (float *)malloc(nx*ny*sizeof(float));
	sinh2kh= (float *)malloc(nx*ny*sizeof(float));

	dhdx=(float *)malloc(nx*ny*sizeof(float));
	dhdy=(float *)malloc(nx*ny*sizeof(float));
	dudx=(float *)malloc(nx*ny*sizeof(float));
	dudy=(float *)malloc(nx*ny*sizeof(float));
	dvdx=(float *)malloc(nx*ny*sizeof(float));
	dvdy=(float *)malloc(nx*ny*sizeof(float));
	
	C=(float *)malloc(nx*ny*sizeof(float));
	R=(float *)malloc(nx*ny*sizeof(float));
	DR=(float *)malloc(nx*ny*sizeof(float));
	
	// Init GPU

	CUDA_SAFE_CALL(cudaSetDevice(GPUDEVICE));
	
	//CUT_DEVICE_INIT(argc, argv);
	printf("Allocating GPU memory\n");


//printf("Loading GPU twisters configurations(Random number generator)...\n");
//   
//      //const char *raw_path = cutFindFilePath("MersenneTwister.raw", argv[0]);
//      //const char *dat_path = cutFindFilePath("MersenneTwister.dat", argv[0]);
//    	printf("loadMTGPU\n");
//    	loadMTGPU("MersenneTwister.dat");
//     	 printf("seedMTGPU\n");
//   	 seedMTGPU(SEED);


	//Allocate GPU memory
	CUDA_SAFE_CALL( cudaMalloc((void **)&hh_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&uu_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vv_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&wci_g, nx*ny*sizeof(float )) );
	
	CUDA_SAFE_CALL( cudaMalloc((void **)&ueu_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vev_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vmageu_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vmagev_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&zs_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&zb_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&Ceq_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dzsdx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dzsdy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&cfm_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&fwm_g, nx*ny*sizeof(float )) );
		CUDA_SAFE_CALL( cudaMalloc((void **)&dzsdt_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&hu_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&hv_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&hum_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&hvm_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vu_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&uv_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&wetu_g, nx*ny*sizeof(int )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&wetv_g, nx*ny*sizeof(int )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&ududx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vdudy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vdvdy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&udvdx_g, nx*ny*sizeof(float )) );
	
	

	CUDA_SAFE_CALL( cudaMalloc((void **)&D_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&urms_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&ust_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Fx_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Fy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&k_g, nx*ny*sizeof(float )) );
	

	CUDA_SAFE_CALL( cudaMalloc((void **)&ee_g, nx*ny*ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&rr_g, nx*ny*ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&St_g, ny*ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&sigm_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&DR_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&R_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&H_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&qbndold_g, 3*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&qbndnew_g, 3*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&umeanbnd_g, ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vmeanbnd_g, ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&sigt_g, nx*ny*ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&c_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&cg_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&theta_g, ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&cxsth_g, ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&sxnth_g, ntheta*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&thetamean_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&Sxx_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&Syy_g, nx*ny*sizeof(float )) );
	//CUDA_SAFE_CALL( cudaMalloc((void **)&Sxy_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&ctheta_g, nx*ny*ntheta*sizeof(float )) );



	CUDA_SAFE_CALL( cudaMalloc((void **)&stdep_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Cc_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&dzb_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&kturb_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&rolthick_g, nx*ny*sizeof(float )) );





	// Averaged variables
	CUDA_SAFE_CALL( cudaMalloc((void **)&Hmean_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&uumean_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&vvmean_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&hhmean_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&zsmean_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&Cmean_g, nx*ny*sizeof(float )) );
	CUDA_SAFE_CALL( cudaMalloc((void **)&ceqsg_g, nx*ny*sizeof(float )) );





////////////////////////////////////////////////////////////////////////////////////////////
// Copy CPU array to the GPU                                                       /////////
////////////////////////////////////////////////////////////////////////////////////////////



	CUDA_SAFE_CALL( cudaMemcpy(hh_g, hh, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(hu_g, hh, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(hv_g, hh, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(hum_g, hh, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(hvm_g, hh, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(zb_g, zb, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(zs_g, zs, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(uu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vv_g, vv, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	//CUDA_SAFE_CALL( cudaMemcpy(zom_g, zom, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vu_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(uv_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vev_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(ueu_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vmageu_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vmagev_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(ududx_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vdudy_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vdvdy_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(udvdx_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );


	CUDA_SAFE_CALL( cudaMemcpy(H_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(Fx_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(Fy_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(urms_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(ust_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(D_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );



	CUDA_SAFE_CALL( cudaMemcpy(ee_g, ee, nx*ny*ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(rr_g, rr, nx*ny*ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(cxsth_g, cxsth, ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(sxnth_g, sxnth, ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(theta_g, theta, ntheta*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(thetamean_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(R_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(DR_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(c_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(cg_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	
	
	CUDA_SAFE_CALL( cudaMemcpy(stdep_g,stdep, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(Cc_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(kturb_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(rolthick_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(dzb_g,dzb, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	CUDA_SAFE_CALL( cudaMemcpy(umeanbnd_g,umeanbnd, ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vmeanbnd_g,umeanbnd, ny*sizeof(float ), cudaMemcpyHostToDevice) );

	//Averaged variables
	CUDA_SAFE_CALL( cudaMemcpy(Hmean_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(uumean_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(vvmean_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(hhmean_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(zsmean_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(Cmean_g,uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	

	/*CUDA_SAFE_CALL( cudaMemcpy(xxp_g, xxp, npart*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(yyp_g, yyp, npart*sizeof(float ), cudaMemcpyHostToDevice) );

	*/
if (imodel==1 || imodel>2)
{
	CUDA_SAFE_CALL( cudaMemcpy(qbndold_g,qbndold, 3*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(qbndnew_g,qbndnew, 3*ny*sizeof(float ), cudaMemcpyHostToDevice) );
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


	dim3 blockDim(16, 16, 1);
	dim3 gridDim(nx / blockDim.x, ny / blockDim.y, 1);
	
	
	
	
	//Calculate bottomm friction based on initial hard layer file
	updatezom<<<gridDim, blockDim, 0>>>(nx,ny,cf,cf2,fw,fw2,stdep_g,cfm_g,fwm_g);
	CUT_CHECK_ERROR("UpdateZom execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
if (imodel==1 || imodel>2)
{	
	set_bnd<<<gridDim, blockDim, 0>>>(nx,ny,Trep,ntheta,theta_g,sigm_g);

	CUT_CHECK_ERROR("set_bnd() execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
}	
	// prepare output file
	printf("prepare output");
	creatncfile(tsoutfile, nx,ny,dx,0.0f,imodel,zb,zs,uu,vv,H,H,thetamean,uu,uu,uu,uu,uu,uu,uu,hh,uu,uu,uu,uu,uu,uu);

	//create3dnc(nx,ny,ntheta,dx,0.0f,theta,ctheta);

	istepout=istepout+nstepout;
	printf("...done\n");
	

	printf("Starting Computation \n");
	
	// for hot start purposes
	//read2Dnc(nx,ny,"uufile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(uu_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"vvfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(vv_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"zsfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(zs_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	//read2Dnc(nx,ny,"hhfile.nc",uu);
	//CUDA_SAFE_CALL( cudaMemcpy(hh_g, uu, nx*ny*sizeof(float ), cudaMemcpyHostToDevice) );
	
	
	// Run the model
	mainloop();


	
	

	//close the bnd files and clean up a bit
	fclose(fsl);
	fclose(fwind);
	
	endtime=clock();
	
	printf("Total runtime= %d  seconds\n",(endtime-starttime)/CLOCKS_PER_SEC);

	cudaThreadExit();

	//CUT_EXIT(argc, argv);

	








}

