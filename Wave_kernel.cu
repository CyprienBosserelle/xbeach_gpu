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

__device__ int mminus(int ix,int nx)
{
	int xminus;
	if (ix==0)
	{
		xminus=0;
	}	
	else
	{
		xminus=ix-1;	
	}
	return(xminus);
}
__device__ int pplus(int ix, int nx)
{
	int xplus;
	if (ix==nx-1)
	{
		xplus=nx-1;
	}	
	else
	{
		xplus=ix+1;	
	}
	return(xplus);

}

__device__ int mminus2(int ix,int nx)
{
	int xminus;
	if (ix<=1)
	{
		xminus=0;
	}	
	else
	{
		xminus=ix-2;	
	}
	return(xminus);
}
__device__ int pplus2(int ix, int nx)
{
	int xplus;
	if (ix>=nx-2)
	{
		xplus=nx-1;
	}	
	else
	{
		xplus=ix+2;	
	}
	return(xplus);

}

__device__ int sign(float x)
{
	return((x>0.0f) - (x<0.0f));
}


__global__ void addavg_var(int nx, int ny,float * Varmean,float * Var)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	__shared__ float mvari[16][16];
	__shared__ float vari[16][16];

	mvari[tx][ty]=Varmean[i];
	vari[tx][ty]=Var[i];

	Varmean[i]=mvari[tx][ty]+vari[tx][ty];


}

__global__ void divavg_var(int nx, int ny,float ntdiv,float * Varmean)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	__shared__ float mvari[16][16];
	mvari[tx][ty]=Varmean[i];
	Varmean[i]=mvari[tx][ty]/ntdiv;
	

}

__global__ void resetavg_var(int nx, int ny,float * Varmean)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	Varmean[i]=0.0f;
}

__global__ void offshorebndWav(int nx, int ny,int ntheta,float totaltime,float Trep,float *St,float *sigm, float *ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	float taper=min(totaltime/100.0f,1.0f);
	
	sigm[i]=2*pi/Trep;
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		ee[0+iy*nx+itheta*ny*nx]=St[iy+itheta*ny]*taper;
		
	}

}

__global__ void set_bnd(int nx, int ny,float Trep,int ntheta,float * theta,float *sigm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;

	
	
	
	//sigmm=2*pi/Trep;


	sigm[i]=2.0f*pi/Trep;
	
	
			//for (int itheta=0; itheta<ntheta; itheta++)
			//{
			//	sigt[i+itheta*nx*ny] = sigmm;
				//thet[i+itheta*nx*ny] = theta[itheta];
				
			//}
			
}

__global__ void inituv(int nx, int ny,float *uu,float *vv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int yminus=mminus(iy,ny);
	

	float dummy;
	uu[i]=0.0f;
	vv[i]=0.0f;

	//if (iy>0)
	//{
	dummy=vv[ix+(yminus)*nx];
	vv[ix+(iy)*nx]=dummy;
	//}
	
			
}

__global__ void sanity(int nx, int ny,float eps,float * hh,float * sigm, int ntheta,float * ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	

	//__shared__ float hhi[16][16];
	//__shared__ float ssigm[16][16];
	//ssigm[tx][ty]=0.0f;
	//hhi[tx][ty]=hh[i];
	//hh[i]=max(hhi[tx][ty],eps);
	

    for (int itheta=0; itheta<ntheta; itheta++)
	{
		//ssigm[tx][ty]=ssigm[tx][ty]+sigt[i+itheta*nx*ny]/ntheta;
		ee[i+itheta*nx*ny]=max(ee[i+itheta*nx*ny],0.0f);
    }
	//sigm[i]=max(ssigm[tx][ty],0.0001f);
        


}

__global__ void dispersion(int nx,int ny,float twopi,float g,float aphi,float bphi,float * sigm,float * hh,float * k,float * c,float * kh,float * sinh2kh,float * cg)
{
	float L0, L1,L2, errdisp;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	

	__shared__ float sigmi[16][16];
	__shared__ float hhi[16][16];
	sigmi[tx][ty]=sigm[i];
	hhi[tx][ty]=hh[i];
  
        
        L0=twopi*g/(sigmi[tx][ty]*sigmi[tx][ty]);
        L1=L0;
        
        errdisp=1000.0f;
        //while (errdisp > 0.0001f)
        for (int k=1; k<200; k++)
        {
           L2        = L0*tanh(2*pi*hhi[tx][ty]/L1);        
            errdisp       = abs(L2 - L1);
            L1 = (L1*aphi + L2*bphi);//          ! Golden ratio
            if (errdisp <= 0.0001f)
            {
				break;
			}
			if(k==199)
			{
				L1=L0*powf(tanh(powf(sigmi[tx][ty]*sigmi[tx][ty]*hhi[tx][ty]/g,3.0f/4.0f)),2.0f/3.0f);
				break;
			}
        }
        

		//L1=L0*powf(tanh(powf(sigmi[tx][ty]*sigmi[tx][ty]*hhi[tx][ty]/g,3/4)),2/3);
        
		float kk=2*pi/L1;
		k[i]  = kk;
		float cc=sigmi[tx][ty]/kk;
        c[i]  = cc;
		float kkhh=min(kk*hhi[tx][ty],10.0f);
        kh[i]   = kkhh;
		float s2kh=sinhf(2.0f*kkhh);
        sinh2kh[i]=s2kh;
		cg[i] = cc*(0.5f+kkhh/s2kh);
 
}

__device__ float slopes2Dx(int nx,float dx,int i,int ix, int iy,float * hh)
{


	float dhdx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	



	dhdx=(hh[xplus+iy*nx]-hh[xminus+iy*nx])/((xplus-xminus)*dx);
	return(dhdx);
}
__device__ float slopes2Dy(int nx,int ny,float dx,int i,int ix, int iy,float * hh)
{
	float dhdy;
	
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);



	dhdy=(hh[ix+yplus*nx]-hh[ix+yminus*nx])/((yplus-yminus)*dx);
	return(dhdy);
}

__global__ void slopes(int nx,int ny,float dx,float * hh,float * uu,float * vv,float * dhdx,float * dhdy,float * dudx,float * dudy,float * dvdx,float * dvdy)//
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);


			
            dhdx[i]=slopes2Dx(nx,dx,i,ix,iy,hh);
            dudx[i]=slopes2Dx(nx,dx,i,ix,iy,uu);
            dvdx[i]=slopes2Dx(nx,dx,i,ix,iy,vv);

			dhdy[i]=slopes2Dy(nx,ny,dx,i,ix,iy,hh);
       	    dudy[i]=slopes2Dy(nx,ny,dx,i,ix,iy,uu);
			dvdy[i]=slopes2Dy(nx,ny,dx,i,ix,iy,vv);
      
        
	
}

				
__global__ void propagtheta(int nx,int ny,int ntheta,float * wci,float *ctheta,/*float *c,float * cx,float *cy,*/float *cxsth,float *sxnth,/*float *uu,float *vv,*/float *dhdx,float *dhdy,float *dudx,float *dudy,float *dvdx,float *dvdy,float *sigm,float *kh)//
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	float ctht;

	__shared__ float sigmiosinh2kh[16][16];
	 __shared__ float dhdxi[16][16];
	 __shared__ float dhdyi[16][16];
	 __shared__ float dudxi[16][16];
	 __shared__ float dudyi[16][16];
	 __shared__ float dvdxi[16][16];
	 __shared__ float dvdyi[16][16];

	//float ci=c[i];
	//float uui=uu[i]*wci;
	//float vvi=vv[i]*wci;
    sigmiosinh2kh[tx][ty]=sigm[i]/sinhf(2.0f*kh[i]);
	 
	dhdxi[tx][ty]=dhdx[i];
	dhdyi[tx][ty]=dhdy[i];
	dudxi[tx][ty]=dudx[i];
	dudyi[tx][ty]=dudy[i];
	dvdxi[tx][ty]=dvdx[i];
	dvdyi[tx][ty]=dvdy[i];
	 __syncthreads();

	for (int itheta=0; itheta<ntheta; itheta++)
		{

			//cx[i+itheta*nx*ny] =ci*cxsth[itheta]+uui;
			//cy[i+itheta*nx*ny] =ci*sxnth[itheta]+vvi;
			ctht= (sigmiosinh2kh[tx][ty])*(dhdxi[tx][ty]*sxnth[itheta]-dhdyi[tx][ty]*cxsth[itheta]) + wci[i]*(cxsth[itheta]*(sxnth[itheta]*dudxi[tx][ty] - cxsth[itheta]*dudyi[tx][ty]) + sxnth[itheta]*(sxnth[itheta]*dvdxi[tx][ty] - cxsth[itheta]*dvdyi[tx][ty]));
			ctheta[i+itheta*nx*ny]=min(max(ctht,-0.25*sigm[i]),0.25*sigm[i]);
		}

}

__global__ void action(int ntheta,int nx,int ny,float * ee,float * sigm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	__shared__ float ssigm[16][16];
	ssigm[tx][ty]=sigm[i];
	
	
	for (int itheta=0; itheta<ntheta; itheta++)
	{	
		ee[i+itheta*nx*ny] = ee[i+itheta*nx*ny]/ssigm[tx][ty]; 
	}

}
		

__global__ void xadvecupwind(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float * wci,float *ee,float *cg,float *cxsth,float *uu,float * xadvec)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	float dxplus_i = 1.0f/dx;
	float dxcent_i = 1.0f/(2*dx);
	float xxadvec;
	float costhet;
	float arrinx, arrminx, arrmaxx;
	float cgx,cgxmin;
	__shared__ float ccg[16][16];
	__shared__ float ccgxmin[16][16];
	__shared__ float ccgxmax[16][16];
	__shared__ float uui[16][16];

	__shared__ float uuixmin[16][16];
	__shared__ float uuixmax[16][16];
	//__shared__ float eet[16][16];
	
	unsigned int xminus=mminus(ix,nx);//max(ix-1,0);
	unsigned int xplus=pplus(ix,nx);//min(ix+1,nx-1);
	unsigned int yminus=mminus(iy,ny);//max(iy-1,0);
	unsigned int yplus=pplus(iy,ny);//min(iy+1,ny-1);

	
	ccg[tx][ty]=cg[i];
	ccgxmin[tx][ty]=cg[xminus+iy*nx];
	ccgxmax[tx][ty]=cg[xplus+iy*nx];
	
	
	uui[tx][ty]=uu[i]*wci[i];
	uuixmin[tx][ty]=uu[xminus+iy*nx]*wci[i];
	uuixmax[tx][ty]=uu[xplus+iy*nx]*wci[i];
	
	__syncthreads();




	

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		costhet=cxsth[itheta];
		
		cgx=0.5f*(ccg[tx][ty]*costhet+uui[tx][ty]+ccgxmax[tx][ty]*costhet+uuixmax[tx][ty]);
		cgxmin=0.5f*(ccg[tx][ty]*costhet+uui[tx][ty]+ccgxmin[tx][ty]*costhet+uuixmin[tx][ty]);
		xxadvec=0;
		//eet[tx][ty]=ee[i+itheta*nx*ny];
		
		arrinx=	ee[i+itheta*nx*ny]*max(cgx,0.0f)+ee[xplus+iy*nx+itheta*nx*ny]*min(cgx,0.0f);
		arrminx=ee[xminus+iy*nx+itheta*nx*ny]*max(cgxmin,0.0f)+ee[i+itheta*nx*ny]*min(cgxmin,0.0f);
		
		/*if (cgx>0)
		{			
			arrinx=ee[i+itheta*nx*ny]*cgx;	
		}
		else
		{
			arrinx=ee[xplus+iy*nx+itheta*nx*ny]*cgx;	
		}
		if (cgxmin>0)
		{
			arrminx=ee[xminus+iy*nx+itheta*nx*ny]*cgxmin;
		}
		else
		{
			arrminx=ee[i+itheta*nx*ny]*cgxmin;
		}*/

		xxadvec=(arrinx-arrminx)*dxplus_i;
		xadvec[i+itheta*nx*ny]=xxadvec;
	}
	
}




__global__ void xadvecupwind2(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float * wci,float *ee,float *cg,float *cxsth,float *uu,float * xadvec)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	float dxplus_i = 1.0f/dx;
	float dxcent_i = 1.0f/(2*dx);
	float xxadvec;
	float costhet;
	float arrinx, arrminx, arrmaxx;
	float cgx,cgxmin;
	__shared__ float ccg[16][16];
	__shared__ float ccgxmin[16][16];
	__shared__ float ccgxmax[16][16];
	__shared__ float uui[16][16];

	__shared__ float uuixmin[16][16];
	__shared__ float uuixmax[16][16];
	
	unsigned int xminus=mminus(ix,nx);//max(ix-1,0);
	unsigned int xplus=pplus(ix,nx);//min(ix+1,nx-1);
	unsigned int xplus2=pplus2(ix,nx);//min(ix+1,nx-1);
	unsigned int xminus2=mminus2(ix,nx);//max(ix-1,0);
	unsigned int yminus=mminus(iy,ny);//max(iy-1,0);
	unsigned int yplus=pplus(iy,ny);//min(iy+1,ny-1);
	
	
	ccg[tx][ty]=cg[i];
	ccgxmin[tx][ty]=cg[xminus+iy*nx];
	ccgxmax[tx][ty]=cg[xplus+iy*nx];
	
	
	uui[tx][ty]=uu[i]*wci[i];
	uuixmin[tx][ty]=uu[xminus+iy*nx]*wci[i];
	uuixmax[tx][ty]=uu[xplus+iy*nx]*wci[i];
	
	__syncthreads();




	

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		costhet=cxsth[itheta];
		
		cgx=0.5f*(ccg[tx][ty]*costhet+uui[tx][ty]+ccgxmax[tx][ty]*costhet+uuixmax[tx][ty]);
		cgxmin=0.5f*(ccg[tx][ty]*costhet+uui[tx][ty]+ccgxmin[tx][ty]*costhet+uuixmin[tx][ty]);
		xxadvec=0;
			
		if (cgx>0.0f)
		{			
			arrinx=(1.5*ee[i+itheta*nx*ny]-0.5*ee[xminus+iy*nx+itheta*nx*ny]);
			if (arrinx<0.0f)
			{
				arrinx=ee[i+itheta*nx*ny];
			}
			arrinx=arrinx*cgx;	
		}
		else
		{
			arrinx=1.5*ee[xplus+iy*nx+itheta*nx*ny]-0.5*ee[xplus2+iy*nx+itheta*nx*ny];
			if (arrinx<0.0f)
			{
				arrinx=ee[xplus+iy*nx+itheta*nx*ny];
			}
			arrinx=arrinx*cgx;	
		}
		if (cgxmin>0.0f)
		{
			arrminx=1.5*ee[xminus+iy*nx+itheta*nx*ny]-0.5*ee[xminus2+iy*nx+itheta*nx*ny];
			if (arrminx<0.0f)
			{
				arrminx=ee[xminus+iy*nx+itheta*nx*ny];
			}
			arrminx=arrminx*cgxmin;
		}
		else
		{
			arrminx=1.5*ee[i+itheta*nx*ny]-0.5*ee[xplus+iy*nx+itheta*nx*ny];
			if (arrminx<0.0f)
			{
				arrminx=ee[i+itheta*nx*ny];
			}
			arrminx=arrminx*cgxmin;
		}

		xxadvec=(arrinx-arrminx)*dxplus_i;
		xadvec[i+itheta*nx*ny]=xxadvec;
	}
	
	
}


		
__global__ void yadvecupwind(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float * wci,float *ee,float *cg,float *sxnth,float *vv,float * yadvec){

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	float dxplus_i = 1.0f/dx;
	float dxcent_i = 1.0f/(2.0f*dx);
	float yyadvec;
	float sinthet;
	float  arriny, arrminy, arrmaxy;
	float cgy,cgymin;
	__shared__ float ccg[16][16];
	__shared__ float ccgymin[16][16];
	__shared__ float ccgymax[16][16];
	
	__shared__ float vvi[16][16];
	__shared__ float vviymin[16][16];
	__shared__ float vviymax[16][16];

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);


		
	ccg[tx][ty]=cg[i];
	ccgymin[tx][ty]=cg[ix+(yminus)*nx];
	ccgymax[tx][ty]=cg[ix+(yplus)*nx];
		
	vvi[tx][ty]=wci[i]*vv[i];
	vviymin[tx][ty]=wci[i]*vv[ix+(yminus)*nx];
	vviymax[tx][ty]=wci[i]*vv[ix+(yplus)*nx];
	__syncthreads();

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		
		sinthet=sxnth[itheta];
		yyadvec=0;
		cgy=0.5f*(ccg[tx][ty]*sinthet+vvi[tx][ty]+ccgymax[tx][ty]*sinthet+vviymax[tx][ty]);
		cgymin=0.5f*(ccg[tx][ty]*sinthet+vvi[tx][ty]+ccgymin[tx][ty]*sinthet+vviymin[tx][ty]);
		
		arriny=ee[i+itheta*nx*ny]*max(cgy,0.0f)+ee[ix+yplus*nx+itheta*nx*ny]*min(cgy,0.0f);
		arrminy=ee[ix+yminus*nx+itheta*nx*ny]*max(cgymin,0.0f)+ee[i+itheta*nx*ny]*min(cgymin,0.0f);
		/*		
		if (cgy>0)
		{
			arriny=ee[i+itheta*nx*ny]*cgy;
		}
		else
		{
			arriny=ee[ix+yplus*nx+itheta*nx*ny]*cgy;
		}
		if (cgymin>0)
		{
			arrminy=ee[ix+yminus*nx+itheta*nx*ny]*cgymin;
		}
		else
		{
			arrminy=ee[i+itheta*nx*ny]*cgymin;
		}*/
		yyadvec=(arriny-arrminy)*dxplus_i;		
		yadvec[i+itheta*nx*ny]=yyadvec;
	}
	
	

}



__global__ void yadvecupwind2(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float * wci,float *ee,float *cg,float *sxnth,float *vv,float * yadvec)
{

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	float dxplus_i = 1.0f/dx;
	float dxcent_i = 1.0f/(2.0f*dx);
	float yyadvec;
	float sinthet;
	float  arriny, arrminy, arrmaxy;
	float cgy,cgymin;
	__shared__ float ccg[16][16];
	__shared__ float ccgymin[16][16];
	__shared__ float ccgymax[16][16];
	
	__shared__ float vvi[16][16];
	__shared__ float vviymin[16][16];
	__shared__ float vviymax[16][16];

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	unsigned int yminus2=mminus2(iy,ny);
	unsigned int yplus2=pplus2(iy,ny);


		
	ccg[tx][ty]=cg[i];
	ccgymin[tx][ty]=cg[ix+(yminus)*nx];
	ccgymax[tx][ty]=cg[ix+(yplus)*nx];
		
	vvi[tx][ty]=wci[i]*vv[i];
	vviymin[tx][ty]=wci[i]*vv[ix+(yminus)*nx];
	vviymax[tx][ty]=wci[i]*vv[ix+(yplus)*nx];
	__syncthreads();

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		
		sinthet=sxnth[itheta];
		yyadvec=0;
		cgy=0.5f*(ccg[tx][ty]*sinthet+vvi[tx][ty]+ccgymax[tx][ty]*sinthet+vviymax[tx][ty]);
		cgymin=0.5f*(ccg[tx][ty]*sinthet+vvi[tx][ty]+ccgymin[tx][ty]*sinthet+vviymin[tx][ty]);
		
				
		if (cgy>0.0f)
		{
			arriny=1.5*ee[i+itheta*nx*ny]-0.5*ee[ix+yminus*nx+itheta*nx*ny];
			if (arriny<0.0f)
			{
				arriny=ee[i+itheta*nx*ny];
			}
			arriny=arriny*cgy;
		}
		else
		{
			arriny=1.5*ee[ix+yplus*nx+itheta*nx*ny]-0.5*ee[ix+yplus2*nx+itheta*nx*ny];
			if (arriny<0.0f)
			{
				arriny=ee[ix+yplus*nx+itheta*nx*ny];
			}
			arriny=arriny*cgy;
		}
		if (cgymin>0.0f)
		{
			arrminy=1.5*ee[ix+yminus*nx+itheta*nx*ny]-0.5*ee[ix+yminus2*nx+itheta*nx*ny];
			if (arrminy<0.0f)
			{
				arrminy=ee[ix+yminus*nx+itheta*nx*ny];
			}
			arrminy=arrminy*cgymin;
		}
		else
		{
			arrminy=1.5*ee[i+itheta*nx*ny]-0.5*ee[ix+yplus*nx+itheta*nx*ny];
			if(arrminy<0.0f)
			{
				arrminy=ee[i+itheta*nx*ny];
			}
			arrminy=arrminy*cgymin;
		}
		yyadvec=(arriny-arrminy)*dxplus_i;		
		yadvec[i+itheta*nx*ny]=yyadvec;
	}
	

}

__global__ void eectheta(int nx,int ny,int ntheta,float *ee,float *ctheta,float *eect)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		eect[i+itheta*nx*ny]=ee[i+itheta*nx*ny]*ctheta[i+itheta*nx*ny];
	}
		
	
}

__global__ void thetaadvecuw(int nx,int ny,int ntheta,float dtheta,float *eect,float * thetaadvec)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	float arrint,arrmint,arrmaxt;
	
	float tthetaadvec;
	
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		unsigned int Tminus=mminus(itheta,ntheta);
		unsigned int Tplus=pplus(itheta,ntheta);
		tthetaadvec=0;
		
		arrint=eect[i+itheta*nx*ny];
		arrmint=eect[i+Tminus*nx*ny]*(itheta>1);
		arrmaxt=eect[i+Tplus*nx*ny]*(itheta<ntheta-2);
		
		tthetaadvec=((arrint-arrmint)*(arrint>0)+(arrmaxt-arrint)*(arrint<0))/dtheta;
		
		
		thetaadvec[i+itheta*nx*ny]=tthetaadvec;
	}
}
		
__global__ void thetaadvecupwind(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float wci,float *ee,float *ctheta,float * thetaadvec){
	
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	float dxplus_i = 1.0f/dx;
	float dxcent_i = 1.0f/(2.0f*dx);
	float tthetaadvec;
	float costhet,sinthet;
	float arrint,arrmint,arrmaxt;
	
	
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	

	//tthetaadvec=((arrint-arrmint)*max(sign(arrint),0.0f)+(arrint-arrmaxt)*min(sign(arrint),0.0f))/dtheta;
	
	
	
	__syncthreads();




	

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		unsigned int Tminus=mminus(itheta,ntheta);
		unsigned int Tplus=pplus(itheta,ntheta);
		tthetaadvec=0;
		
		arrint=ee[i+itheta*nx*ny]*ctheta[i+itheta*nx*ny];
			
		if (arrint>0)
		{
			if (itheta==0)
			{
				arrmint=0.0f;
			}
			else
			{
				arrmint=ee[i+(Tminus)*nx*ny]*ctheta[i+(Tminus)*nx*ny];
			}
			tthetaadvec=(arrint-arrmint)/dtheta;
		}
		else
		{
			if (arrint<0)
			{
				if (itheta==ntheta-1)
				{
				arrmaxt=0.0f;
				}
				else
				{

				arrmaxt=ee[i+(Tplus)*nx*ny]*ctheta[i+(Tplus)*nx*ny];
				}
				tthetaadvec=(arrmaxt-arrint)/dtheta;
			}
			else
			{
				arrmint=ee[i+(Tminus)*nx*ny]*ctheta[i+(Tminus)*nx*ny];
				arrmaxt=ee[i+(Tplus)*nx*ny]*ctheta[i+(Tplus)*nx*ny];
				tthetaadvec=(arrmaxt-arrmint)/(2*dtheta);
			}
		}
        
		thetaadvec[i+itheta*nx*ny]=tthetaadvec;
	}
	
		
}

__global__ void thetaadvecupwind2(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float wci,float *ee,float *ctheta,float * thetaadvec)
{
	
		unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	
	float dxplus_i = 1.0f/dx;
	float dxcent_i = 1.0f/(2.0f*dx);
	float tthetaadvec;
	float costhet,sinthet;
	float arrint,arrmint,arrmaxt;
	
	
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	

	//tthetaadvec=((arrint-arrmint)*max(sign(arrint),0.0f)+(arrint-arrmaxt)*min(sign(arrint),0.0f))/dtheta;
	
	
	
	__syncthreads();




	

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		unsigned int Tminus=mminus(itheta,ntheta);
		unsigned int Tplus=pplus(itheta,ntheta);
		tthetaadvec=0;
		
		arrint=ee[i+itheta*nx*ny]*ctheta[i+itheta*nx*ny];
			
		if (arrint>0)
		{
			if (itheta==0)
			{
				arrmint=0.0f;
			}
			else
			{
				arrmint=ee[i+(Tminus)*nx*ny]*ctheta[i+(Tminus)*nx*ny];
			}
			tthetaadvec=(arrint-arrmint)/dtheta;
		}
		else
		{
			if (arrint<0)
			{
				if (itheta==ntheta-1)
				{
				arrmaxt=0.0f;
				}
				else
				{

				arrmaxt=ee[i+(Tplus)*nx*ny]*ctheta[i+(Tplus)*nx*ny];
				}
				tthetaadvec=(arrmaxt-arrint)/dtheta;
			}
			else
			{
				arrmint=ee[i+(Tminus)*nx*ny]*ctheta[i+(Tminus)*nx*ny];
				arrmaxt=ee[i+(Tplus)*nx*ny]*ctheta[i+(Tplus)*nx*ny];
				tthetaadvec=(arrmaxt-arrmint)/(2*dtheta);
			}
		}
        
		thetaadvec[i+itheta*nx*ny]=tthetaadvec;
	}
	
		
}


		
__global__ void eulerupwind(int nx,int ny,int ntheta,float dtheta,float dx,float dt,float wci,float *ee,float *xadvec,float *yadvec,float * thetaadvec){



	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		ee[i+itheta*nx*ny]=ee[i+itheta*nx*ny]-dt*(xadvec[i+itheta*nx*ny]+yadvec[i+itheta*nx*ny]+thetaadvec[i+itheta*nx*ny]);
	}
	

}

__global__ void energy(int nx,int ny,int ntheta,float * ee,float * sigm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	__shared__ float eetmp[16][16];
	__shared__ float ssigm[16][16];
	ssigm[tx][ty]=sigm[i];

	for (int itheta=0; itheta<ntheta; itheta++)
		{
			eetmp[tx][ty] = ee[i+itheta*nx*ny]*ssigm[tx][ty];
			
			ee[i+itheta*nx*ny]=max(eetmp[tx][ty],0.0f);
			
		}
}

__global__ void energint(int nx,int ny,int ntheta,float dtheta,float rho,float g,float gammax,float * E,float * H,float * hh,float * ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;
	__shared__ float Etmp[16][16];
	__shared__ float Htmp[16][16];
	__shared__ float hhtmp[16][16];
	
	Etmp[tx][ty]=0.0f;
	hhtmp[tx][ty]=hh[i];
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		Etmp[tx][ty]=Etmp[tx][ty]+ee[i+itheta*nx*ny];
	}
	Etmp[tx][ty]=Etmp[tx][ty]*dtheta;
	Htmp[tx][ty]=sqrtf(Etmp[tx][ty]/(rho*g/8.0f));//Hrms
	float idiv=max(1.0f,powf(Htmp[tx][ty]/(gammax*hhtmp[tx][ty]),2.0f));
	
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		ee[i+itheta*nx*ny]=ee[i+itheta*nx*ny]/idiv;
	}
	Htmp[tx][ty]=min(Htmp[tx][ty],gammax*hhtmp[tx][ty]);
	E[i]=(rho*g*Htmp[tx][ty]*Htmp[tx][ty])/8.0f;
	H[i]=Htmp[tx][ty];

	
}

__global__ void calctm(int nx, int ny,int ntheta, float * tm,float * theta,float * ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	//float sume=0;
	//float sumthet=0;

	__shared__ float sume[16][16];
	__shared__ float sumthet[16][16];

	sume[tx][ty]=0.0f;
	sumthet[tx][ty]=0.0f;




	for (int itheta=0; itheta<ntheta; itheta++)
	{
		sume[tx][ty]=sume[tx][ty]+ee[i+itheta*nx*ny];
		sumthet[tx][ty]=sumthet[tx][ty]+ee[i+itheta*nx*ny]*theta[itheta];
	}
	
	tm[i] = (sumthet[tx][ty]/ntheta)/(max(sume[tx][ty],0.00001)/ntheta);
}

__global__ void roelvink(int nx, int ny,float rho,float g,float gamma,float alpha,float n,float Trep,float * fwm,float * cfm,float *hh,float *H,float *E,float *D, float *k)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	float fac,hroelvink,Qb;
	float htmp,etmp;
	//int bk=3; //1 or 3
	
    fac=8.0f/rho/g;
    hroelvink=hh[i];
	etmp=E[i];

    htmp=sqrtf(fac*etmp);


    Qb=1-exp(-1*pow((htmp/gamma/hroelvink),n));

	Qb=min(Qb,1.0f);


 	D[i]=Qb*2*alpha/Trep*etmp;
	
	float fw;
	float sqrt2=1.4142136;
	float urms=pi*htmp/(Trep*sqrt2*sinh(k[i]*hroelvink));
	float uorb=pi*htmp/(Trep*sinhf(min(k[i]*hroelvink,10.0f)));
			/* float Ab=urms*Trep/(2*pi);
			 float kb=30*cfm[i];
			 float Abkb=Ab/kb;

			 if (Abkb<=0.2)
			 {
				 fw=0.5f;
			 }
			 else
			 {
				 if(Abkb>100)
				 {
					 fw=exp(-7.30+5.61*powf(Abkb,-0.109f));
				 }
				 else
				 {
					 fw=exp(-8.82+7.02*powf(Abkb,-0.078f));
				 }
			 }*/
			 fw=fwm[i];
			 D[i]=D[i]+0.66666666f/pi*rho*fw*uorb*uorb*uorb;

}

__global__ void baldock(int nx, int ny,float rho,float g,float gamma,float alpha,float n,float Trep,float *fwm,float * cfm,float *hh,float *H,float *E,float *D, float *k)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	float fac,hbald,Qb;
	float Htmp,etmp;
	float kh,Hb,R,ki;	
//int bk=3; //1 or 3
	
    //fac=8.0f/rho/g;
    hbald=hh[i];
	//etmp=E[i];
	ki=k[i];

    Htmp=H[i];
    
    kh=ki*hbald;

   Hb=(0.88f/ki)*tanhf(gamma*kh/0.88f);

	R=Hb/max(Htmp,0.00001f);
    Qb=expf(-1.0f*R*R);
    
    D[i]=0.25f*alpha*Qb*rho*(1.0f/Trep)*g*(Hb*Hb+Htmp*Htmp);

//ADD dissipation due to bottom friction
	
			// float fw;
			 //float sqrt2=1.4142136f;
			 //float urms=pi*Htmp/(Trep*sqrt2*sinh(min(max(k[i],0.01f)*hbald,10.0f)));
			 float uorb=pi*Htmp/(Trep*sinhf(min(max(k[i],0.01f)*hbald,10.0f)));
			 //float urms=uorb/sqrt2;
			/* float Ab=urms*Trep/(2.0f*pi);
			 float kb=30.0f*zo[i];
			 float Abkb=Ab/kb;

			 if (Abkb<=0.2)
			 {
				 fw=0.5;
			 }
			 else
			 {
				 if(Abkb>100)
				 {
					 fw=exp(-7.30+5.61*pow(Abkb,-0.109f));
				 }
				 else
				 {
					 fw=exp(-8.82+7.02*pow(Abkb,-0.078f));
				 }
			 }
			 //fw=0.01f;
			*/
			D[i]=D[i]+0.66666666f/pi*rho*fwm[i]*uorb*uorb*uorb;
			
	
}



__global__ void rollerlatbnd(int nx,int ny,int ntheta,float eps,float *hh,float *rr)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	
	__shared__ float hhi[16][16];
	__shared__ float rri[16][16];
	__shared__ float rrt[16][16];
	__shared__ float rrb[16][16];
	__shared__ float rrr[16][16];


	hhi[tx][ty]=hh[i];
	float wet=0.0f;

	if (hhi[tx][ty]>eps)
	{
		wet=1.0f;
	}
	
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		//rri[tx][ty]=rr[i+itheta*nx*ny];
		//rrt[tx][ty]=rr[ix+yplus*nx+itheta*nx*ny];
		//rrb[tx][ty]=rr[ix+yminus*nx+itheta*nx*ny];
		//rrr[tx][ty]=rr[xplus+iy*nx+itheta*nx*ny];

		//rr[i+itheta*nx*ny]=rri[tx][ty]*wet;
		if (iy==0)
		{
			rr[i+itheta*nx*ny]=rr[ix+yplus*nx+itheta*nx*ny]*wet;
		}
		if (iy==ny-1)
		{
			rr[i+itheta*nx*ny]=rr[ix+yminus*nx+itheta*nx*ny]*wet;
		}
		//if (ix==0)
		//{
		//	rr[i+itheta*nx*ny]=rri[tx][ty]*wet;
		//}
	}
			



}

__global__ void twodimbnd(int nx,int ny,float eps,float * hh,float * F)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	
	__shared__ float hhi[16][16];
	__shared__ float Fi[16][16];
	__shared__ float Ft[16][16];
	__shared__ float Fb[16][16];
	__shared__ float Fr[16][16];

	
	


	hhi[tx][ty]=hh[i];
	float wet=0.0f;

	if (hhi[tx][ty]>eps)
	{
		wet=1.0f;
	}
	

		Fi[tx][ty]=F[i];
		Ft[tx][ty]=F[ix+yplus*nx];
		Fb[tx][ty]=F[ix+yminus*nx];
		Fr[tx][ty]=F[xplus+iy*nx];

		//F[i]=Fi[tx][ty]*wet;
		if (iy==0)
		{
			F[i]=Ft[tx][ty]*wet;
		}
		if (iy==ny-1)
		{
			F[i]=Fb[tx][ty]*wet;
		}
		if (ix==0)
		{
			F[i]=Fr[tx][ty]*wet;
		}

}

__global__ void dissipation(int nx,int ny,int ntheta,float dtheta,float eps,float dt,float g,float beta,float * wci,float *hh,float *ee,float *D,float *E,float *rr,float *c,float *cxsth,float *sxnth,float * uu,float * vv,float *DR,float *R)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	float dd;
	float RR=0.0f;
	float drr=0.0f;
	float cx,cy;
	//float wci=0.0f;
	float rrd;
	float eet,rrt;

	__shared__ float cc[16][16];
	__shared__ float uui[16][16];
	__shared__ float vvi[16][16];

	cc[tx][ty]=c[i];
	uui[tx][ty]=uu[i]*wci[i];
	vvi[tx][ty]=vv[i]*wci[i];
	__syncthreads();

	for (int itheta=0; itheta<ntheta; itheta++)
	{
			cx=cc[tx][ty]*cxsth[itheta]+uui[tx][ty];
			cy=cc[tx][ty]*sxnth[itheta]+vvi[tx][ty];

			eet=ee[i+itheta*nx*ny];
			dd=eet*D[i]/max(E[i],0.0000001f);

			
			
			
			if(hh[i]>eps)
			{
				rrt=rr[i+itheta*nx*ny];
				rrd=2.0f*g*beta*rrt/sqrtf(cx*cx+cy*cy);
				eet=max(eet-dt*dd,0.0f);
				rrt=max(rrt+dt*(dd-rrd),0.0f);
				drr=drr+rrd;   
			}
			else
			{
				eet=0.0f;
				rrt=0.0f;
				drr=0.0f;

			}
			RR=RR+rrt;
			ee[i+itheta*nx*ny]=eet;
			rr[i+itheta*nx*ny]=rrt;
	}
	R[i]=RR*dtheta;
	DR[i]=drr*dtheta;
}

__global__ void meandir(int nx,int ny,int ntheta,float rho,float g,float dtheta,float * ee,float * thet, float * thetamean, float * E,float * H)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	
	float sumethet=0;
	
	float sume=0;
	for (int itheta=0; itheta<ntheta; itheta++)
	{
		sume=sume+ee[i+itheta*nx*ny];
		sumethet=sumethet+ee[i+itheta*nx*ny]*thet[itheta];
	}
	sume=max(sume,0.00001f);
	__syncthreads;
	thetamean[i]=(sumethet/ntheta)/(sume/ntheta);
	E[i]=sume*dtheta;

	H[i]=sqrtf(sume*dtheta/(rho*g/8));//sqrt(E[i]/(1/8*rho*g));
}

__global__ void radstress(int nx,int ny, int ntheta,float dx,float dtheta,float * ee,float *rr,float * cxsth,float * sxnth,float * cg,float * c,float * Sxx,float * Sxy,float * Syy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;
	


	 __shared__ float n[16][16];
	 

	n[tx][ty]=cg[i]/c[i];
	__syncthreads;
	float sume=0.0f;
	float tmpxx=0.0f;
	float tmpxy=0.0f;
	float tmpyy=0.0f;
	float rolxx=0.0f;
	float rolxy=0.0f;
	float rolyy=0.0f;

	for (int itheta=0; itheta<ntheta; itheta++)
	{
		sume=sume+ee[i+itheta*nx*ny];
		tmpxx=tmpxx+((1.0f+cxsth[itheta]*cxsth[itheta])*ee[i+itheta*nx*ny]);
		tmpyy=tmpyy+((1.0f+sxnth[itheta]*sxnth[itheta])*ee[i+itheta*nx*ny]);
		tmpxy=tmpxy+(sxnth[itheta]*cxsth[itheta]*ee[i+itheta*nx*ny]);
		rolxx=rolxx+((cxsth[itheta]*cxsth[itheta])*rr[i+itheta*nx*ny]);
		rolyy=rolyy+((sxnth[itheta]*sxnth[itheta])*rr[i+itheta*nx*ny]);
		rolxy=rolxy+(sxnth[itheta]*cxsth[itheta]*rr[i+itheta*nx*ny]);

	}
	

	Sxx[i]=(n[tx][ty]*tmpxx-0.5f*sume)*dtheta+rolxx*dtheta;
	Syy[i]=(n[tx][ty]*tmpyy-0.5f*sume)*dtheta+rolyy*dtheta;
	Sxy[i]=n[tx][ty]*tmpxy*dtheta+rolxy*dtheta;
}

__global__ void wavforce(int nx,int ny, int ntheta,float dx,float dtheta,float * Sxx,float * Sxy,float * Syy,float * Fx,float * Fy,float * hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	
	//float hmin=0.000002f;
	
	__shared__ float  Sxxi[16][16];
	__shared__ float  Sxxr[16][16];
	__shared__ float  Sxyi[16][16];
	__shared__ float  Sxyt[16][16];
	__shared__ float  Sxytr[16][16];
	__shared__ float  Sxyb[16][16];
	__shared__ float  Sxybr[16][16];
	
	__shared__ float  Syyi[16][16];
	__shared__ float  Syyt[16][16];
	__shared__ float  Sxyr[16][16];
	__shared__ float  Sxyl[16][16];
	__shared__ float  Sxytl[16][16];

	__shared__ float FFx[16][16];
	__shared__ float FFy[16][16] ;
	
	


	Sxxi[tx][ty]=Sxx[i];
	Sxxr[tx][ty]=Sxx[xplus+iy*nx];
	Sxyi[tx][ty]=Sxy[i];
	Sxyt[tx][ty]=Sxy[ix+yplus*nx];
	Sxytr[tx][ty]=Sxy[xplus+yplus*nx];
	Sxyb[tx][ty]=Sxy[ix+yminus*nx];
	Sxybr[tx][ty]=Sxy[xplus+yminus*nx];

	Syyi[tx][ty]=Syy[i];
	Syyt[tx][ty]=Syy[ix+yplus*nx];
	Sxyr[tx][ty]=Sxy[xplus+iy*nx];
	Sxyl[tx][ty]=Sxy[xminus+iy*nx];
	Sxytl[tx][ty]=Sxy[xminus+yplus*nx];

	FFx[tx][ty]=0.0f;
	FFy[tx][ty]=0.0f;


	__syncthreads;
	
	
		//if(hh[i]>hmin)
		//{
			FFx[tx][ty]=-1*(Sxxr[tx][ty]-Sxxi[tx][ty])/dx-0.5f*(Sxyt[tx][ty]+Sxytr[tx][ty]-Sxyb[tx][ty]-Sxybr[tx][ty])/(2.0f*dx);
		
			FFy[tx][ty]=-1*(Syyt[tx][ty]-Syyi[tx][ty])/dx-0.5f*(Sxyr[tx][ty]+Sxytr[tx][ty]-Sxyl[tx][ty]-Sxytl[tx][ty])/(2.0f*dx);
		//}
	
		
		Fx[i]=FFx[tx][ty];

		Fy[i]=FFy[tx][ty];
		

}

__global__ void breakerdelay(int nx,int ny,int ntheta,float dtheta,float g, float rho,float Trep, float eps,float * urms,float *ust,float *H,float *E,float *c,float *k, float *hh,float *rr)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	
	float hmin=0.2;
	float delta=0.0f;
	float Hi=H[i];
	float ki=k[i];
	float hhi=max(hh[i],hmin);
	float ci=max(c[i],sqrt(hmin*g));
	float R;
	float ustw,uwf,vwf,ustr,usd,uorb;
	
	R=rr[i];

	


	 
	uorb=pi*Hi/Trep/sinhf(min(max(ki,0.01f)*hhi,10.0f));

	urms[i]=uorb/(sqrtf(2.0f));
	ustw=E[i]/ci/rho/hhi; 
	//uwf = ustw*cosf(tm[i]);
	//vwf = ustw*sinf(tm[i]);
	// roller contribution
	
	ustr=2.0f*R/ci/rho/hhi;
	// introduce breaker delay
	//usd=ustr; //I don't like breaker delay


	ust[i]=ustw+ustr;

}

__global__ void calcwci(int nx, int ny, float wci, float hwci, float * hh, float *wcig)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	
	
	wcig[i]=wci*min(hh[i]/hwci,1.0f);
	
}



