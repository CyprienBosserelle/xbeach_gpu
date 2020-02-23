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

__device__ int mminus(int ix, int nx)
{
	int xminus;
	if (ix == 0)
	{
		xminus = 0;
	}
	else
	{
		xminus = ix - 1;
	}
	return(xminus);
}
__device__ int pplus(int ix, int nx)
{
	int xplus;
	if (ix == nx - 1)
	{
		xplus = nx - 1;
	}
	else
	{
		xplus = ix + 1;
	}
	return(xplus);

}

__device__ int mminus2(int ix, int nx)
{
	int xminus;
	if (ix <= 1)
	{
		xminus = 0;
	}
	else
	{
		xminus = ix - 2;
	}
	return(xminus);
}
__device__ int pplus2(int ix, int nx)
{
	int xplus;
	if (ix >= nx - 2)
	{
		xplus = nx - 1;
	}
	else
	{
		xplus = ix + 2;
	}
	return(xplus);

}

__device__ int sign(DECNUM x)
{
	return((x > 0.0f) - (x < 0.0f));
}


__global__ void addavg_var(int nx, int ny, DECNUM * Varmean, DECNUM * Var)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	__shared__ DECNUM mvari[16][16];
	__shared__ DECNUM vari[16][16];

	if (ix < nx && iy < ny)
	{

		mvari[tx][ty] = Varmean[i];
		vari[tx][ty] = Var[i];

		Varmean[i] = mvari[tx][ty] + vari[tx][ty];
	}


}


__global__ void divavg_var(int nx, int ny, DECNUM ntdiv, DECNUM * Varmean)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	__shared__ DECNUM mvari[16][16];
	if (ix < nx && iy < ny)
	{
		mvari[tx][ty] = Varmean[i];
		Varmean[i] = mvari[tx][ty] / ntdiv;
	}


}

__global__ void resetavg_var(int nx, int ny, DECNUM * Varmean)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	if (ix < nx && iy < ny)
	{
		Varmean[i] = 0.0f;
	}
}

__global__ void max_var(int nx, int ny, DECNUM * Varmax, DECNUM * Var)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	__shared__ DECNUM mvari[16][16];
	__shared__ DECNUM vari[16][16];

	if (ix < nx && iy < ny)
	{

		mvari[tx][ty] = Varmax[i];
		vari[tx][ty] = Var[i];

		Varmax[i] = max(mvari[tx][ty] , vari[tx][ty]);
	}


}

__global__ void FLOWDT(int nx, int ny, DECNUM dx, DECNUM cfl, DECNUM *dtflow, DECNUM *hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	if (ix < nx && iy < ny)
	{
		dtflow[i] = cfl*dx / (sqrtf(9.81f*hh[i]));
	}
}

__global__ void minmaxKernel(int ntot, DECNUM *max, DECNUM *min, DECNUM *a) {
	__shared__ double maxtile[32];
	__shared__ double mintile[32];

	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < ntot)
	{
		maxtile[tid] = a[i];
		mintile[tid] = a[i];
		__syncthreads();

		// strided index and non-divergent branch
		for (unsigned int s = 1; s < blockDim.x; s *= 2) {
			int index = 2 * s * tid;
			if (index < blockDim.x) {
				if (maxtile[tid + s] > maxtile[tid])
					maxtile[tid] = maxtile[tid + s];
				if (mintile[tid + s] < mintile[tid])
					mintile[tid] = mintile[tid + s];
			}
			__syncthreads();
		}

		if (tid == 0) {
			max[blockIdx.x] = maxtile[0];
			min[blockIdx.x] = mintile[0];
		}
	}
}

__global__ void finalminmaxKernel(DECNUM *max, DECNUM *min) {
	__shared__ double maxtile[32];
	__shared__ double mintile[32];

	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	maxtile[tid] = max[i];
	mintile[tid] = min[i];
	__syncthreads();

	// strided index and non-divergent branch
	for (unsigned int s = 1; s < blockDim.x; s *= 2) {
		int index = 2 * s * tid;
		if (index < blockDim.x) {
			if (maxtile[tid + s] > maxtile[tid])
				maxtile[tid] = maxtile[tid + s];
			if (mintile[tid + s] < mintile[tid])
				mintile[tid] = mintile[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0) {
		max[blockIdx.x] = maxtile[0];
		min[blockIdx.x] = mintile[0];
	}
}


__global__ void WAVEDT(int nx, int ny, int ntheta, DECNUM cfl, DECNUM dtheta, DECNUM *dtwave, DECNUM *ctheta)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int iz = blockIdx.z*blockDim.z + threadIdx.z;
	unsigned int i = ix + iy*nx + iz*nx*ny;

	float mindt=9999999.9f;
	if (ix < nx && iy < ny)
	{
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			mindt = min(cfl*dtheta / abs(ctheta[i + itheta*nx*ny]), mindt);
		}

		dtwave[i] = mindt;
	}
}

__global__ void offshorebndWav(int nx, int ny, int ntheta, DECNUM totaltime, DECNUM Trep, DECNUM *St, DECNUM *sigm, DECNUM *ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	DECNUM taper = min(totaltime / 100.0f, 1.0f);


	if (ix < nx && iy < ny)
	{
		sigm[i] = 2 * pi / Trep;
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			ee[0 + iy*nx + itheta*ny*nx] = St[iy + itheta*ny] * taper;

		}
	}

}

__global__ void set_bnd(int nx, int ny, DECNUM Trep, int ntheta, DECNUM * theta, DECNUM *sigm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;




	//sigmm=2*pi/Trep;

	if (ix < nx && iy < ny)
	{
		sigm[i] = 2.0f*pi / Trep;
	}


	//for (int itheta=0; itheta<ntheta; itheta++)
	//{
	//	sigt[i+itheta*nx*ny] = sigmm;
	//thet[i+itheta*nx*ny] = theta[itheta];

	//}

}

__global__ void inituv(int nx, int ny, DECNUM *uu, DECNUM *vv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix < nx && iy < ny)
	{

		int yminus = mminus(iy, ny);


		DECNUM dummy;


		uu[i] = 0.0f;
		vv[i] = 0.0f;

		//if (iy>0)
		//{
		dummy = vv[ix + (yminus)*nx];
		vv[ix + (iy)*nx] = dummy;
		//}
	}

}

__global__ void sanity(int nx, int ny, DECNUM eps, DECNUM * hh, DECNUM * sigm, int ntheta, DECNUM * ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	//int tx = threadIdx.x;
	//int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		//__shared__ DECNUM hhi[16][16];
		//__shared__ DECNUM ssigm[16][16];
		//ssigm[tx][ty]=0.0f;
		//hhi[tx][ty]=hh[i];
		//hh[i]=max(hhi[tx][ty],eps);


		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			//ssigm[tx][ty]=ssigm[tx][ty]+sigt[i+itheta*nx*ny]/ntheta;
			ee[i + itheta*nx*ny] = max(ee[i + itheta*nx*ny], 0.0f);
		}
		//sigm[i]=max(ssigm[tx][ty],0.0001f);
	}



}

__global__ void dispersion(int nx, int ny, DECNUM twopi, DECNUM g, DECNUM aphi, DECNUM bphi, DECNUM * sigm, DECNUM * hh, DECNUM * k, DECNUM * c, DECNUM * kh, DECNUM * sinh2kh, DECNUM * cg)
{
	DECNUM L0, L1, L2, errdisp;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;


	__shared__ DECNUM sigmi[16][16];
	__shared__ DECNUM hhi[16][16];

	if (ix < nx && iy < ny)
	{
		sigmi[tx][ty] = sigm[i];
		hhi[tx][ty] = hh[i];


		L0 = twopi*g / (sigmi[tx][ty] * sigmi[tx][ty]);
		L1 = L0;

		errdisp = 1000.0f;
		//while (errdisp > 0.0001f)
		for (int p = 1; p < 200; p++)
		{
			L2 = L0*tanh(2 * pi*hhi[tx][ty] / L1);
			errdisp = abs(L2 - L1);
			L1 = (L1*aphi + L2*bphi);//          ! Golden ratio
			if (errdisp <= 0.0001f)
			{
				break;
			}
			if (p == 199)
			{
				L1 = L0*powf(tanh(powf(sigmi[tx][ty] * sigmi[tx][ty] * hhi[tx][ty] / g, 3.0f / 4.0f)), 2.0f / 3.0f);
				break;
			}
		}


		//L1=L0*powf(tanh(powf(sigmi[tx][ty]*sigmi[tx][ty]*hhi[tx][ty]/g,3/4)),2/3);

		DECNUM kk = 2 * pi / L1;
		k[i] = kk;
		DECNUM cc = sigmi[tx][ty] / kk;
		c[i] = cc;
		DECNUM kkhh = min(kk*hhi[tx][ty], 10.0f);
		kh[i] = kkhh;
		DECNUM s2kh = sinhf(2.0f*kkhh);
		sinh2kh[i] = s2kh;
		cg[i] = cc*(0.5f + kkhh / s2kh);
	}

}

__global__ void dispersion_init(int nx, int ny, DECNUM twopi, DECNUM g, DECNUM aphi, DECNUM bphi, DECNUM * sigm, DECNUM * hh, DECNUM * cg)
{
	DECNUM L0, L1, L2, errdisp;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;


	__shared__ DECNUM sigmi[16][16];
	__shared__ DECNUM hhi[16][16];

	if (ix < nx && iy < ny)
	{
		sigmi[tx][ty] = sigm[i];
		hhi[tx][ty] = hh[i];


		L0 = twopi*g / (sigmi[tx][ty] * sigmi[tx][ty]);
		L1 = L0;

		errdisp = 1000.0f;
		//while (errdisp > 0.0001f)
		for (int k = 1; k < 200; k++)
		{
			L2 = L0*tanh(2 * pi*hhi[tx][ty] / L1);
			errdisp = abs(L2 - L1);
			L1 = (L1*aphi + L2*bphi);//          ! Golden ratio
			if (errdisp <= 0.0001f)
			{
				break;
			}
			if (k == 199)
			{
				L1 = L0*powf(tanh(powf(sigmi[tx][ty] * sigmi[tx][ty] * hhi[tx][ty] / g, 3.0f / 4.0f)), 2.0f / 3.0f);
				break;
			}
		}


		//L1=L0*powf(tanh(powf(sigmi[tx][ty]*sigmi[tx][ty]*hhi[tx][ty]/g,3/4)),2/3);

		DECNUM kk = 2 * pi / L1;
		//k[i]  = kk;
		DECNUM cc = sigmi[tx][ty] / kk;
		//c[i]  = cc;
		DECNUM kkhh = min(kk*hhi[tx][ty], 10.0f);
		// kh[i]   = kkhh;
		DECNUM s2kh = sinhf(2.0f*kkhh);
		//sinh2kh[i]=s2kh;
		cg[i] = cc*(0.5f + kkhh / s2kh);
	}

}

__device__ DECNUM slopes2Dx(int nx, DECNUM dx, int i, int ix, int iy, DECNUM * hh)
{


	DECNUM dhdx = 0.0f;

	if (ix < nx)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);




		dhdx = (hh[xplus + iy*nx] - hh[xminus + iy*nx]) / ((xplus - xminus)*dx);
	}
	return(dhdx);
}
__device__ DECNUM slopes2Dy(int nx, int ny, DECNUM dx, int i, int ix, int iy, DECNUM * hh)
{
	DECNUM dhdy = 0.0f;
	if (ix < nx && iy < ny)
	{
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		dhdy = (hh[ix + yplus*nx] - hh[ix + yminus*nx]) / ((yplus - yminus)*dx);
	}
	return(dhdy);
}

__global__ void slopes(int nx, int ny, DECNUM dx, DECNUM * hh, DECNUM * uu, DECNUM * vv, DECNUM * dhdx, DECNUM * dhdy, DECNUM * dudx, DECNUM * dudy, DECNUM * dvdx, DECNUM * dvdy)//
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);



		dhdx[i] = slopes2Dx(nx, dx, i, ix, iy, hh);
		dudx[i] = slopes2Dx(nx, dx, i, ix, iy, uu);
		dvdx[i] = slopes2Dx(nx, dx, i, ix, iy, vv);

		dhdy[i] = slopes2Dy(nx, ny, dx, i, ix, iy, hh);
		dudy[i] = slopes2Dy(nx, ny, dx, i, ix, iy, uu);
		dvdy[i] = slopes2Dy(nx, ny, dx, i, ix, iy, vv);

	}



}


__global__ void propagtheta(int nx, int ny, int ntheta, DECNUM * wci, DECNUM *ctheta,/*DECNUM *c,DECNUM * cx,DECNUM *cy,*/DECNUM *cxsth, DECNUM *sxnth,/*DECNUM *uu,DECNUM *vv,*/DECNUM *dhdx, DECNUM *dhdy, DECNUM *dudx, DECNUM *dudy, DECNUM *dvdx, DECNUM *dvdy, DECNUM *sigm, DECNUM *kh)//
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	DECNUM ctht;

	__shared__ DECNUM sigmiosinh2kh[16][16];
	__shared__ DECNUM dhdxi[16][16];
	__shared__ DECNUM dhdyi[16][16];
	__shared__ DECNUM dudxi[16][16];
	__shared__ DECNUM dudyi[16][16];
	__shared__ DECNUM dvdxi[16][16];
	__shared__ DECNUM dvdyi[16][16];


	if (ix < nx && iy < ny)
	{
		//DECNUM ci=c[i];
		//DECNUM uui=uu[i]*wci;
		//DECNUM vvi=vv[i]*wci;
		sigmiosinh2kh[tx][ty] = sigm[i] / sinhf(2.0f*kh[i]);

		dhdxi[tx][ty] = dhdx[i];
		dhdyi[tx][ty] = dhdy[i];
		dudxi[tx][ty] = dudx[i];
		dudyi[tx][ty] = dudy[i];
		dvdxi[tx][ty] = dvdx[i];
		dvdyi[tx][ty] = dvdy[i];
		__syncthreads();

		for (int itheta = 0; itheta < ntheta; itheta++)
		{

			//cx[i+itheta*nx*ny] =ci*cxsth[itheta]+uui;
			//cy[i+itheta*nx*ny] =ci*sxnth[itheta]+vvi;
			ctht = (sigmiosinh2kh[tx][ty])*(dhdxi[tx][ty] * sxnth[itheta] - dhdyi[tx][ty] * cxsth[itheta]) + wci[i] * (cxsth[itheta] * (sxnth[itheta] * dudxi[tx][ty] - cxsth[itheta] * dudyi[tx][ty]) + sxnth[itheta] * (sxnth[itheta] * dvdxi[tx][ty] - cxsth[itheta] * dvdyi[tx][ty]));
			ctheta[i + itheta*nx*ny] = min(max(ctht, -0.25*sigm[i]), 0.25*sigm[i]);
		}
	}

}

__global__ void action(int ntheta, int nx, int ny, DECNUM * ee, DECNUM * sigm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	__shared__ DECNUM ssigm[16][16];

	if (ix < nx && iy < ny)
	{
		ssigm[tx][ty] = sigm[i];


		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			ee[i + itheta*nx*ny] = ee[i + itheta*nx*ny] / ssigm[tx][ty];
		}
	}

}


__global__ void xadvecupwind(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *cxsth, DECNUM *uu, DECNUM * xadvec)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2 * dx);
	DECNUM xxadvec;
	DECNUM costhet;
	DECNUM arrinx, arrminx, arrmaxx;
	DECNUM cgx, cgxmin;
	__shared__ DECNUM ccg[16][16];
	__shared__ DECNUM ccgxmin[16][16];
	__shared__ DECNUM ccgxmax[16][16];
	__shared__ DECNUM uui[16][16];

	__shared__ DECNUM uuixmin[16][16];
	__shared__ DECNUM uuixmax[16][16];
	//__shared__ DECNUM eet[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);//max(ix-1,0);
		unsigned int xplus = pplus(ix, nx);//min(ix+1,nx-1);
		unsigned int yminus = mminus(iy, ny);//max(iy-1,0);
		unsigned int yplus = pplus(iy, ny);//min(iy+1,ny-1);


		ccg[tx][ty] = cg[i];
		ccgxmin[tx][ty] = cg[xminus + iy*nx];
		ccgxmax[tx][ty] = cg[xplus + iy*nx];


		uui[tx][ty] = uu[i] * wci[i];
		uuixmin[tx][ty] = uu[xminus + iy*nx] * wci[i];
		uuixmax[tx][ty] = uu[xplus + iy*nx] * wci[i];

		__syncthreads();






		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			costhet = cxsth[itheta];

			cgx = 0.5f*(ccg[tx][ty] * costhet + uui[tx][ty] + ccgxmax[tx][ty] * costhet + uuixmax[tx][ty]);
			cgxmin = 0.5f*(ccg[tx][ty] * costhet + uui[tx][ty] + ccgxmin[tx][ty] * costhet + uuixmin[tx][ty]);
			xxadvec = 0;
			//eet[tx][ty]=ee[i+itheta*nx*ny];

			arrinx = ee[i + itheta*nx*ny] * max(cgx, 0.0f) + ee[xplus + iy*nx + itheta*nx*ny] * min(cgx, 0.0f);
			arrminx = ee[xminus + iy*nx + itheta*nx*ny] * max(cgxmin, 0.0f) + ee[i + itheta*nx*ny] * min(cgxmin, 0.0f);

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

			xxadvec = (arrinx - arrminx)*dxplus_i;
			xadvec[i + itheta*nx*ny] = xxadvec;
		}
	}
}




__global__ void xadvecupwind2(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *cxsth, DECNUM *uu, DECNUM * xadvec)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2 * dx);
	DECNUM xxadvec;
	DECNUM costhet;
	DECNUM arrinx, arrminx, arrmaxx;
	DECNUM cgx, cgxmin;
	__shared__ DECNUM ccg[16][16];
	__shared__ DECNUM ccgxmin[16][16];
	__shared__ DECNUM ccgxmax[16][16];
	__shared__ DECNUM uui[16][16];

	__shared__ DECNUM uuixmin[16][16];
	__shared__ DECNUM uuixmax[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);//max(ix-1,0);
		unsigned int xplus = pplus(ix, nx);//min(ix+1,nx-1);
		unsigned int xplus2 = pplus2(ix, nx);//min(ix+1,nx-1);
		unsigned int xminus2 = mminus2(ix, nx);//max(ix-1,0);
		unsigned int yminus = mminus(iy, ny);//max(iy-1,0);
		unsigned int yplus = pplus(iy, ny);//min(iy+1,ny-1);


		ccg[tx][ty] = cg[i];
		ccgxmin[tx][ty] = cg[xminus + iy*nx];
		ccgxmax[tx][ty] = cg[xplus + iy*nx];


		uui[tx][ty] = uu[i] * wci[i];
		uuixmin[tx][ty] = uu[xminus + iy*nx] * wci[i];
		uuixmax[tx][ty] = uu[xplus + iy*nx] * wci[i];

		__syncthreads();






		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			costhet = cxsth[itheta];

			cgx = 0.5f*(ccg[tx][ty] * costhet + uui[tx][ty] + ccgxmax[tx][ty] * costhet + uuixmax[tx][ty]);
			cgxmin = 0.5f*(ccg[tx][ty] * costhet + uui[tx][ty] + ccgxmin[tx][ty] * costhet + uuixmin[tx][ty]);
			xxadvec = 0;



			if (cgx > 0.0f)
			{
				arrinx = (1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[xminus + iy*nx + itheta*nx*ny]);
				if (arrinx < 0.0f)
				{
					arrinx = ee[i + itheta*nx*ny];
				}
				arrinx = arrinx*cgx;
			}
			else
			{
				arrinx = 1.5f*ee[xplus + iy*nx + itheta*nx*ny] - 0.5f*ee[xplus2 + iy*nx + itheta*nx*ny];
				if (arrinx < 0.0f)
				{
					arrinx = ee[xplus + iy*nx + itheta*nx*ny];
				}
				arrinx = arrinx*cgx;
			}
			if (cgxmin > 0.0f)
			{
				arrminx = 1.5f*ee[xminus + iy*nx + itheta*nx*ny] - 0.5f*ee[xminus2 + iy*nx + itheta*nx*ny];
				if (arrminx < 0.0f)
				{
					arrminx = ee[xminus + iy*nx + itheta*nx*ny];
				}
				arrminx = arrminx*cgxmin;
			}
			else
			{
				arrminx = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[xplus + iy*nx + itheta*nx*ny];
				if (arrminx < 0.0f)
				{
					arrminx = ee[i + itheta*nx*ny];
				}
				arrminx = arrminx*cgxmin;
			}

			xxadvec = (arrinx - arrminx)*dxplus_i;
			xadvec[i + itheta*nx*ny] = xxadvec;
		}
	}

}

__global__ void xadvecupwind2SD(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt,DECNUM * thetamean, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *cxsth, DECNUM *uu, DECNUM * xadvec)
{
	//Same as xadvecupwind2 but using single dir formulation which is much faster
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2 * dx);
	DECNUM xxadvec;
	DECNUM costhet;
	DECNUM thetmi;
	DECNUM arrinx, arrminx, arrmaxx;
	DECNUM cgx, cgxmin;
	__shared__ DECNUM ccg[16][16];
	__shared__ DECNUM ccgxmin[16][16];
	__shared__ DECNUM ccgxmax[16][16];
	__shared__ DECNUM uui[16][16];

	__shared__ DECNUM uuixmin[16][16];
	__shared__ DECNUM uuixmax[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);//max(ix-1,0);
		unsigned int xplus = pplus(ix, nx);//min(ix+1,nx-1);
		unsigned int xplus2 = pplus2(ix, nx);//min(ix+1,nx-1);
		unsigned int xminus2 = mminus2(ix, nx);//max(ix-1,0);
		unsigned int yminus = mminus(iy, ny);//max(iy-1,0);
		unsigned int yplus = pplus(iy, ny);//min(iy+1,ny-1);


		ccg[tx][ty] = cg[i];
		ccgxmin[tx][ty] = cg[xminus + iy*nx];
		ccgxmax[tx][ty] = cg[xplus + iy*nx];


		uui[tx][ty] = uu[i] * wci[i];
		uuixmin[tx][ty] = uu[xminus + iy*nx] * wci[i];
		uuixmax[tx][ty] = uu[xplus + iy*nx] * wci[i];

		__syncthreads();

		//thetmi = thetamean[i];

		costhet = cos(thetamean[i]);
		cgx = 0.5f*(ccg[tx][ty] * costhet + uui[tx][ty] + ccgxmax[tx][ty] * costhet + uuixmax[tx][ty]);
		cgxmin = 0.5f*(ccg[tx][ty] * costhet + uui[tx][ty] + ccgxmin[tx][ty] * costhet + uuixmin[tx][ty]);

		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			

			
			xxadvec = 0;



			if (cgx > 0.0f)
			{
				arrinx = (1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[xminus + iy*nx + itheta*nx*ny]);
				if (arrinx < 0.0f)
				{
					arrinx = ee[i + itheta*nx*ny];
				}
				arrinx = arrinx*cgx;
			}
			else
			{
				arrinx = 1.5f*ee[xplus + iy*nx + itheta*nx*ny] - 0.5f*ee[xplus2 + iy*nx + itheta*nx*ny];
				if (arrinx < 0.0f)
				{
					arrinx = ee[xplus + iy*nx + itheta*nx*ny];
				}
				arrinx = arrinx*cgx;
			}
			if (cgxmin > 0.0f)
			{
				arrminx = 1.5f*ee[xminus + iy*nx + itheta*nx*ny] - 0.5f*ee[xminus2 + iy*nx + itheta*nx*ny];
				if (arrminx < 0.0f)
				{
					arrminx = ee[xminus + iy*nx + itheta*nx*ny];
				}
				arrminx = arrminx*cgxmin;
			}
			else
			{
				arrminx = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[xplus + iy*nx + itheta*nx*ny];
				if (arrminx < 0.0f)
				{
					arrminx = ee[i + itheta*nx*ny];
				}
				arrminx = arrminx*cgxmin;
			}

			xxadvec = (arrinx - arrminx)*dxplus_i;
			xadvec[i + itheta*nx*ny] = xxadvec;
		}
	}

}


__global__ void yadvecupwind(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *sxnth, DECNUM *vv, DECNUM * yadvec){

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM yyadvec;
	DECNUM sinthet;
	DECNUM  arriny, arrminy, arrmaxy;
	DECNUM cgy, cgymin;
	__shared__ DECNUM ccg[16][16];
	__shared__ DECNUM ccgymin[16][16];
	__shared__ DECNUM ccgymax[16][16];

	__shared__ DECNUM vvi[16][16];
	__shared__ DECNUM vviymin[16][16];
	__shared__ DECNUM vviymax[16][16];
	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		ccg[tx][ty] = cg[i];
		ccgymin[tx][ty] = cg[ix + (yminus)*nx];
		ccgymax[tx][ty] = cg[ix + (yplus)*nx];

		vvi[tx][ty] = wci[i] * vv[i];
		vviymin[tx][ty] = wci[i] * vv[ix + (yminus)*nx];
		vviymax[tx][ty] = wci[i] * vv[ix + (yplus)*nx];
		__syncthreads();

		for (int itheta = 0; itheta < ntheta; itheta++)
		{

			sinthet = sxnth[itheta];
			yyadvec = 0;
			cgy = 0.5f*(ccg[tx][ty] * sinthet + vvi[tx][ty] + ccgymax[tx][ty] * sinthet + vviymax[tx][ty]);
			cgymin = 0.5f*(ccg[tx][ty] * sinthet + vvi[tx][ty] + ccgymin[tx][ty] * sinthet + vviymin[tx][ty]);

			arriny = ee[i + itheta*nx*ny] * max(cgy, 0.0f) + ee[ix + yplus*nx + itheta*nx*ny] * min(cgy, 0.0f);
			arrminy = ee[ix + yminus*nx + itheta*nx*ny] * max(cgymin, 0.0f) + ee[i + itheta*nx*ny] * min(cgymin, 0.0f);
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
			yyadvec = (arriny - arrminy)*dxplus_i;
			yadvec[i + itheta*nx*ny] = yyadvec;
		}

	}

}



__global__ void yadvecupwind2(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *sxnth, DECNUM *vv, DECNUM * yadvec)
{

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM yyadvec;
	DECNUM sinthet;
	DECNUM  arriny, arrminy, arrmaxy;
	DECNUM cgy, cgymin;
	__shared__ DECNUM ccg[16][16];
	__shared__ DECNUM ccgymin[16][16];
	__shared__ DECNUM ccgymax[16][16];

	__shared__ DECNUM vvi[16][16];
	__shared__ DECNUM vviymin[16][16];
	__shared__ DECNUM vviymax[16][16];


	if (ix < nx && iy < ny)
	{

		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);



		ccg[tx][ty] = cg[i];
		ccgymin[tx][ty] = cg[ix + (yminus)*nx];
		ccgymax[tx][ty] = cg[ix + (yplus)*nx];

		vvi[tx][ty] = wci[i] * vv[i];
		vviymin[tx][ty] = wci[i] * vv[ix + (yminus)*nx];
		vviymax[tx][ty] = wci[i] * vv[ix + (yplus)*nx];
		__syncthreads();

		for (int itheta = 0; itheta < ntheta; itheta++)
		{

			sinthet = sxnth[itheta];
			yyadvec = 0;
			cgy = 0.5f*(ccg[tx][ty] * sinthet + vvi[tx][ty] + ccgymax[tx][ty] * sinthet + vviymax[tx][ty]);
			cgymin = 0.5f*(ccg[tx][ty] * sinthet + vvi[tx][ty] + ccgymin[tx][ty] * sinthet + vviymin[tx][ty]);


			if (cgy > 0.0f)
			{
				arriny = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[ix + yminus*nx + itheta*nx*ny];
				if (arriny < 0.0f)
				{
					arriny = ee[i + itheta*nx*ny];
				}
				arriny = arriny*cgy;
			}
			else
			{
				arriny = 1.5f*ee[ix + yplus*nx + itheta*nx*ny] - 0.5f*ee[ix + yplus2*nx + itheta*nx*ny];
				if (arriny < 0.0f)
				{
					arriny = ee[ix + yplus*nx + itheta*nx*ny];
				}
				arriny = arriny*cgy;
			}
			if (cgymin > 0.0f)
			{
				arrminy = 1.5f*ee[ix + yminus*nx + itheta*nx*ny] - 0.5f*ee[ix + yminus2*nx + itheta*nx*ny];
				if (arrminy < 0.0f)
				{
					arrminy = ee[ix + yminus*nx + itheta*nx*ny];
				}
				arrminy = arrminy*cgymin;
			}
			else
			{
				arrminy = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[ix + yplus*nx + itheta*nx*ny];
				if (arrminy < 0.0f)
				{
					arrminy = ee[i + itheta*nx*ny];
				}
				arrminy = arrminy*cgymin;
			}
			yyadvec = (arriny - arrminy)*dxplus_i;
			yadvec[i + itheta*nx*ny] = yyadvec;
		}
	}


}

__global__ void yadvecupwind2SD(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * thetamean, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *sxnth, DECNUM *vv, DECNUM * yadvec)
{

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM yyadvec;
	DECNUM sinthet;
	DECNUM  arriny, arrminy, arrmaxy;
	DECNUM cgy, cgymin;
	__shared__ DECNUM ccg[16][16];
	__shared__ DECNUM ccgymin[16][16];
	__shared__ DECNUM ccgymax[16][16];

	__shared__ DECNUM vvi[16][16];
	__shared__ DECNUM vviymin[16][16];
	__shared__ DECNUM vviymax[16][16];


	if (ix < nx && iy < ny)
	{

		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);



		ccg[tx][ty] = cg[i];
		ccgymin[tx][ty] = cg[ix + (yminus)*nx];
		ccgymax[tx][ty] = cg[ix + (yplus)*nx];

		vvi[tx][ty] = wci[i] * vv[i];
		vviymin[tx][ty] = wci[i] * vv[ix + (yminus)*nx];
		vviymax[tx][ty] = wci[i] * vv[ix + (yplus)*nx];
		__syncthreads();

		sinthet = sinf(thetamean[i]);
		cgy = 0.5f*(ccg[tx][ty] * sinthet + vvi[tx][ty] + ccgymax[tx][ty] * sinthet + vviymax[tx][ty]);
		cgymin = 0.5f*(ccg[tx][ty] * sinthet + vvi[tx][ty] + ccgymin[tx][ty] * sinthet + vviymin[tx][ty]);

		for (int itheta = 0; itheta < ntheta; itheta++)
		{

			
			yyadvec = 0;
			


			if (cgy > 0.0f)
			{
				arriny = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[ix + yminus*nx + itheta*nx*ny];
				if (arriny < 0.0f)
				{
					arriny = ee[i + itheta*nx*ny];
				}
				arriny = arriny*cgy;
			}
			else
			{
				arriny = 1.5f*ee[ix + yplus*nx + itheta*nx*ny] - 0.5f*ee[ix + yplus2*nx + itheta*nx*ny];
				if (arriny < 0.0f)
				{
					arriny = ee[ix + yplus*nx + itheta*nx*ny];
				}
				arriny = arriny*cgy;
			}
			if (cgymin > 0.0f)
			{
				arrminy = 1.5f*ee[ix + yminus*nx + itheta*nx*ny] - 0.5f*ee[ix + yminus2*nx + itheta*nx*ny];
				if (arrminy < 0.0f)
				{
					arrminy = ee[ix + yminus*nx + itheta*nx*ny];
				}
				arrminy = arrminy*cgymin;
			}
			else
			{
				arrminy = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[ix + yplus*nx + itheta*nx*ny];
				if (arrminy < 0.0f)
				{
					arrminy = ee[i + itheta*nx*ny];
				}
				arrminy = arrminy*cgymin;
			}
			yyadvec = (arriny - arrminy)*dxplus_i;
			yadvec[i + itheta*nx*ny] = yyadvec;
		}
	}


}

__global__ void eectheta(int nx, int ny, int ntheta, DECNUM *ee, DECNUM *ctheta, DECNUM *eect)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			eect[i + itheta*nx*ny] = ee[i + itheta*nx*ny] * ctheta[i + itheta*nx*ny];
		}
	}

}

__global__ void thetaadvecuw(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM *eect, DECNUM * thetaadvec)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	DECNUM arrint, arrmint, arrmaxt;

	DECNUM tthetaadvec;

	if (ix < nx && iy < ny)
	{
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			unsigned int Tminus = mminus(itheta, ntheta);
			unsigned int Tplus = pplus(itheta, ntheta);
			tthetaadvec = 0;

			arrint = eect[i + itheta*nx*ny];
			arrmint = eect[i + Tminus*nx*ny] * (itheta > 1);
			arrmaxt = eect[i + Tplus*nx*ny] * (itheta < ntheta - 2);

			tthetaadvec = ((arrint - arrmint)*(arrint>0) + (arrmaxt - arrint)*(arrint < 0)) / dtheta;


			thetaadvec[i + itheta*nx*ny] = tthetaadvec;
		}
	}
}

__global__ void thetaadvecupwind(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *ctheta, DECNUM * thetaadvec){

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM tthetaadvec;
	DECNUM costhet, sinthet;
	DECNUM arrint, arrmint, arrmaxt;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		//tthetaadvec=((arrint-arrmint)*max(sign(arrint),0.0f)+(arrint-arrmaxt)*min(sign(arrint),0.0f))/dtheta;



		__syncthreads();






		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			unsigned int Tminus = mminus(itheta, ntheta);
			unsigned int Tplus = pplus(itheta, ntheta);
			tthetaadvec = 0;

			arrint = ee[i + itheta*nx*ny] * ctheta[i + itheta*nx*ny];

			if (arrint > 0)
			{
				if (itheta == 0)
				{
					arrmint = 0.0f;
				}
				else
				{
					arrmint = ee[i + (Tminus)*nx*ny] * ctheta[i + (Tminus)*nx*ny];
				}
				tthetaadvec = (arrint - arrmint) / dtheta;
			}
			else
			{
				if (arrint < 0)
				{
					if (itheta == ntheta - 1)
					{
						arrmaxt = 0.0f;
					}
					else
					{

						arrmaxt = ee[i + (Tplus)*nx*ny] * ctheta[i + (Tplus)*nx*ny];
					}
					tthetaadvec = (arrmaxt - arrint) / dtheta;
				}
				else
				{
					arrmint = ee[i + (Tminus)*nx*ny] * ctheta[i + (Tminus)*nx*ny];
					arrmaxt = ee[i + (Tplus)*nx*ny] * ctheta[i + (Tplus)*nx*ny];
					tthetaadvec = (arrmaxt - arrmint) / (2 * dtheta);
				}
			}

			thetaadvec[i + itheta*nx*ny] = tthetaadvec;
		}
	}


}
__global__ void thetaadvecuw1ho(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *ctheta, DECNUM * thetaadvec){

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM tthetaadvec, cthetab;
	DECNUM costhet, sinthet;
	DECNUM arrint, arrmint, arrmaxt;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		if (ntheta > 1)
		{

			for (int itheta = 0; itheta < ntheta; itheta++)
			{
				unsigned int Tminus = mminus(itheta, ntheta);
				unsigned int Tplus = pplus(itheta, ntheta);

				tthetaadvec = 0.0;


				cthetab = 0.5f*(ctheta[i + itheta*nx*ny] + ctheta[i + Tplus*nx*ny]);
				arrint = max(cthetab, 0.00f)*ee[i + itheta*nx*ny] + min(cthetab, 0.0f)*ee[i + Tplus*nx*ny];





				cthetab = 0.5f*(ctheta[i + Tminus*nx*ny] + ctheta[i + itheta*nx*ny]);
				arrmint = max(cthetab, 0.0f)*ee[i + Tminus*nx*ny] + min(cthetab, 0.0f)*ee[i + itheta*nx*ny];



				tthetaadvec = (arrint - arrmint) / dtheta;


				if (itheta == ntheta - 1)
				{
					tthetaadvec = (0.0f - arrmint) / dtheta;
				}
				if (itheta == 0)
				{
					tthetaadvec = (arrint - 0.0f) / dtheta;
				}


				thetaadvec[i + itheta*nx*ny] = tthetaadvec;
			}
		}
		else
		{
			thetaadvec[i + 0 * nx*ny] = 0.0f;
		}
	}

}

__global__ void thetaadvecuw2ho(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *ctheta, DECNUM * thetaadvec){

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM tthetaadvec, cthetab;
	DECNUM costhet, sinthet;
	DECNUM arrint, arrmint, eupw;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		if (ntheta > 1)
		{

			for (int itheta = 0; itheta < ntheta; itheta++)
			{
				unsigned int Tminus = mminus(itheta, ntheta);
				unsigned int Tplus = pplus(itheta, ntheta);
				unsigned int Tminus2 = mminus2(itheta, ntheta);
				unsigned int Tplus2 = pplus2(itheta, ntheta);


				tthetaadvec = 0.0;


				cthetab = 0.5f*(ctheta[i + itheta*nx*ny] + ctheta[i + Tplus*nx*ny]);
				if (cthetab > 0.0f)
				{
					eupw = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[i + Tminus*nx*ny];
					if (eupw < 0.0f)
					{
						eupw = ee[i + itheta*nx*ny];
					}
				}
				else
				{
					eupw = 1.5f*ee[i + Tplus*nx*ny] - 0.5f*ee[i + Tplus2*nx*ny];
					if (eupw < 0.0f)
					{
						eupw = ee[i + Tplus*nx*ny];
					}
				}
				arrint = cthetab*eupw;






				cthetab = 0.5f*(ctheta[i + Tminus*nx*ny] + ctheta[i + itheta*nx*ny]);
				if (cthetab > 0.0f)
				{
					eupw = 1.5f*ee[i + Tminus*nx*ny] - 0.5f*ee[i + Tminus2*nx*ny];
					if (eupw < 0.0f)
					{
						eupw = ee[i + Tminus*nx*ny];
					}
				}
				else
				{
					eupw = 1.5f*ee[i + itheta*nx*ny] - 0.5f*ee[i + Tplus*nx*ny];
					if (eupw < 0.0f)
					{
						eupw = ee[i + itheta*nx*ny];
					}
				}
				arrmint = cthetab*eupw;

				tthetaadvec = (arrint - arrmint) / dtheta;


				if (itheta == ntheta - 1)
				{
					tthetaadvec = (0.0f - arrmint) / dtheta;
				}
				if (itheta == 0)
				{
					tthetaadvec = (arrint - 0.0f) / dtheta;
				}


				thetaadvec[i + itheta*nx*ny] = tthetaadvec;
			}
		}
		else
		{
			thetaadvec[i + 0 * nx*ny] = 0.0f;
		}
	}

}


__global__ void thetaadvecuw2hoSLOW(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *ctheta, DECNUM * thetaadvec){

	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int iz = blockIdx.z*blockDim.z + threadIdx.z;
	unsigned int i = ix + iy*nx + iz*nx*ny;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int tz = threadIdx.z;
	

	DECNUM dxplus_i = 1.0f / dx;
	DECNUM dxcent_i = 1.0f / (2.0f*dx);
	DECNUM tthetaadvec, cthetab;
	DECNUM costhet, sinthet;
	DECNUM arrint, arrmint, eupwP, eupwM;

	__shared__ DECNUM eek[16][16];
	__shared__ DECNUM eekp[16][16];
	__shared__ DECNUM eekp2[16][16];
	__shared__ DECNUM eekm[16][16];
	__shared__ DECNUM eekm2[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		




		if (ntheta > 1)
		{

			for (int itheta = 0; itheta < ntheta; itheta++)
			{
				unsigned int Tminus = mminus(itheta, ntheta);
				unsigned int Tplus = pplus(itheta, ntheta);
				unsigned int Tminus2 = mminus2(itheta, ntheta);
				unsigned int Tplus2 = pplus2(itheta, ntheta);


				tthetaadvec = 0.0;

				eek[tx][ty] = ee[i+itheta*nx*ny];
				eekp[tx][ty] = ee[i + Tplus*nx*ny];
				eekp2[tx][ty] = ee[i + Tplus2*nx*ny];
				eekm[tx][ty] = ee[i + Tminus*nx*ny];
				eekm2[tx][ty] = ee[i + Tplus*nx*ny];


				cthetab = 0.5f*(ctheta[i + itheta*nx*ny] + ctheta[i + Tplus*nx*ny]);
				eupwM = eek[tx][ty] + max(eekm[tx][ty] - 3.0f*eek[tx][ty], 0.0f) / (eekm[tx][ty] - 3.0f*eek[tx][ty]) * (0.5f*eek[tx][ty] - 0.5f*eekm[tx][ty]);
				eupwP = eekp[tx][ty] + max(eekp2[tx][ty] - 3.0f*eekp[tx][ty], 0.0f) / (eekp2[tx][ty] - 3.0f*eekp[tx][ty]) * (0.5f*eekp[tx][ty] - 0.5f*eekp2[tx][ty]);
				
				arrint = max(cthetab, 0.0f)*eupwM + min(cthetab, 0.0f)*eupwP;

				
				cthetab = 0.5f*(ctheta[i + Tminus*nx*ny] + ctheta[i + itheta*nx*ny]);

				eupwM = eekm[tx][ty] + max(eekm2[tx][ty] - 3.0f*eekm[tx][ty], 0.0f) / (eekm2[tx][ty] - 3.0f*eekm[tx][ty]) * (0.5f*eekm[tx][ty] - 0.5f*eekm2[tx][ty]);
				eupwP = eek[tx][ty] + max(eekp[tx][ty] - 3.0f*eek[tx][ty], 0.0f) / (eekp[tx][ty] - 3.0f*eek[tx][ty]) * (0.5f*eek[tx][ty] - 0.5f*eekp[tx][ty]);

				arrmint = max(cthetab, 0.0f)*eupwM + min(cthetab, 0.0f)*eupwP;


				tthetaadvec = (arrint - arrmint) / dtheta;


				if (itheta == ntheta - 1)
				{
					tthetaadvec = (0.0f - arrmint) / dtheta;
				}
				if (itheta == 0)
				{
					tthetaadvec = (arrint - 0.0f) / dtheta;
				}


				thetaadvec[i + itheta*nx*ny] = tthetaadvec;
			}
		}
		else
		{
			thetaadvec[i + 0 * nx*ny] = 0.0f;
		}
	}

}

__global__ void eulerupwind(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *xadvec, DECNUM *yadvec, DECNUM * thetaadvec){



	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			ee[i + itheta*nx*ny] = ee[i + itheta*nx*ny] - dt*(xadvec[i + itheta*nx*ny] + yadvec[i + itheta*nx*ny] + thetaadvec[i + itheta*nx*ny]);
		}

	}
}

__global__ void energy(int nx, int ny, int ntheta, DECNUM * ee, DECNUM * sigm)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	__shared__ DECNUM eetmp[16][16];
	__shared__ DECNUM ssigm[16][16];


	if (ix < nx && iy < ny)
	{
		ssigm[tx][ty] = sigm[i];

		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			eetmp[tx][ty] = ee[i + itheta*nx*ny] * ssigm[tx][ty];

			ee[i + itheta*nx*ny] = max(eetmp[tx][ty], 0.0f);

		}
	}
}

__global__ void energint(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM rho, DECNUM g, DECNUM gammax, DECNUM * E, DECNUM * H, DECNUM * hh, DECNUM * ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;
	__shared__ DECNUM Etmp[16][16];
	__shared__ DECNUM Htmp[16][16];
	__shared__ DECNUM hhtmp[16][16];



	if (ix < nx && iy < ny)
	{
		Etmp[tx][ty] = 0.0f;
		hhtmp[tx][ty] = hh[i];
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			Etmp[tx][ty] = Etmp[tx][ty] + ee[i + itheta*nx*ny];
		}
		Etmp[tx][ty] = Etmp[tx][ty] * dtheta;
		Htmp[tx][ty] = sqrtf(Etmp[tx][ty] / (rho*g / 8.0f));//Hrms
		DECNUM idiv = max(1.0f, powf(Htmp[tx][ty] / (gammax*hhtmp[tx][ty]), 2.0f));

		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			ee[i + itheta*nx*ny] = ee[i + itheta*nx*ny] / idiv;
		}
		Htmp[tx][ty] = min(Htmp[tx][ty], gammax*hhtmp[tx][ty]);
		E[i] = (rho*g*Htmp[tx][ty] * Htmp[tx][ty]) / 8.0f;
		H[i] = Htmp[tx][ty];

	}
}

__global__ void calctm(int nx, int ny, int ntheta, DECNUM * tm, DECNUM * theta, DECNUM * ee)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	//DECNUM sume=0;
	//DECNUM sumthet=0;

	__shared__ DECNUM sume[16][16];
	__shared__ DECNUM sumthet[16][16];


	if (ix < nx && iy < ny)
	{
		sume[tx][ty] = 0.0f;
		sumthet[tx][ty] = 0.0f;




		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			sume[tx][ty] = sume[tx][ty] + ee[i + itheta*nx*ny];
			sumthet[tx][ty] = sumthet[tx][ty] + ee[i + itheta*nx*ny] * theta[itheta];
		}

		tm[i] = (sumthet[tx][ty] / ntheta) / (max(sume[tx][ty], 0.00001) / ntheta);
	}
}

__global__ void roelvink(int nx, int ny, DECNUM rho, DECNUM g, DECNUM gamma, DECNUM alpha, DECNUM n, DECNUM Trep, DECNUM * fwm, DECNUM * cfm, DECNUM *hh, DECNUM *H, DECNUM *E, DECNUM *D, DECNUM *k)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	DECNUM fac, hroelvink, Qb;
	DECNUM htmp, etmp;
	//int bk=3; //1 or 3
	if (ix < nx && iy < ny)
	{
		fac = 8.0f / rho / g;
		hroelvink = hh[i];
		etmp = E[i];

		htmp = sqrtf(fac*etmp);


		Qb = 1.0f - exp(max(-1.0f*pow((htmp / gamma / hroelvink), n), -100.0f));

		Qb = min(Qb, 1.0f);


		D[i] = Qb*2.0f*alpha / Trep*etmp*htmp / hroelvink;

		DECNUM fw;
		DECNUM sqrt2 = 1.4142136;
		DECNUM urms = pi*htmp / (Trep*sqrt2*sinh(k[i] * hroelvink));
		DECNUM uorb = pi*htmp / (Trep*sinhf(min(k[i] * hroelvink, 10.0f)));
		/* DECNUM Ab=urms*Trep/(2*pi);
		 DECNUM kb=30*cfm[i];
		 DECNUM Abkb=Ab/kb;

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
		fw = fwm[i];
		D[i] = D[i] + 0.66666666f / pi*rho*fw*uorb*uorb*uorb;
	}

}

__global__ void baldock(int nx, int ny, DECNUM rho, DECNUM g, DECNUM gamma, DECNUM alpha, DECNUM n, DECNUM Trep, DECNUM *fwm, DECNUM * cfm, DECNUM *hh, DECNUM *H, DECNUM *E, DECNUM *D, DECNUM *k)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	DECNUM fac, hbald, Qb;
	DECNUM Htmp, etmp;
	DECNUM kh, Hb, R, ki;
	//int bk=3; //1 or 3
	if (ix < nx && iy < ny)
	{
		//fac=8.0f/rho/g;
		hbald = hh[i];
		//etmp=E[i];
		ki = k[i];

		Htmp = H[i];

		kh = ki*hbald;

		Hb = (0.88f / ki)*tanhf(gamma*kh / 0.88f);

		R = Hb / max(Htmp, 0.00001f);
		Qb = expf(-1.0f*R*R);

		D[i] = 0.25f*alpha*Qb*rho*(1.0f / Trep)*g*(Hb*Hb + Htmp*Htmp);

		//ADD dissipation due to bottom friction

		// DECNUM fw;
		//DECNUM sqrt2=1.4142136f;
		//DECNUM urms=pi*Htmp/(Trep*sqrt2*sinh(min(max(k[i],0.01f)*hbald,10.0f)));
		DECNUM uorb = pi*Htmp / (Trep*sinhf(min(max(k[i], 0.01f)*hbald, 10.0f)));
		//DECNUM urms=uorb/sqrt2;
		/* DECNUM Ab=urms*Trep/(2.0f*pi);
		 DECNUM kb=30.0f*zo[i];
		 DECNUM Abkb=Ab/kb;

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
		D[i] = D[i] + 0.66666666f / pi*rho*fwm[i] * uorb*uorb*uorb;

	}
}



__global__ void rollerlatbnd(int nx, int ny, int ntheta, DECNUM eps, DECNUM *hh, DECNUM *rr)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	unsigned int xminus = mminus(ix, nx);
	unsigned int xplus = pplus(ix, nx);
	unsigned int yminus = mminus(iy, ny);
	unsigned int yplus = pplus(iy, ny);

	__shared__ DECNUM hhi[16][16];
	__shared__ DECNUM rri[16][16];
	__shared__ DECNUM rrt[16][16];
	__shared__ DECNUM rrb[16][16];
	__shared__ DECNUM rrr[16][16];


	if (ix < nx && iy < ny)
	{
		hhi[tx][ty] = hh[i];
		DECNUM wet = 0.0f;

		if (hhi[tx][ty] > eps)
		{
			wet = 1.0f;
		}

		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			//rri[tx][ty]=rr[i+itheta*nx*ny];
			//rrt[tx][ty]=rr[ix+yplus*nx+itheta*nx*ny];
			//rrb[tx][ty]=rr[ix+yminus*nx+itheta*nx*ny];
			//rrr[tx][ty]=rr[xplus+iy*nx+itheta*nx*ny];

			//rr[i+itheta*nx*ny]=rri[tx][ty]*wet;
			if (iy == 0)
			{
				rr[i + itheta*nx*ny] = rr[ix + yplus*nx + itheta*nx*ny] * wet;
			}
			if (iy == ny - 1)
			{
				rr[i + itheta*nx*ny] = rr[ix + yminus*nx + itheta*nx*ny] * wet;
			}
			//if (ix==0)
			//{
			//	rr[i+itheta*nx*ny]=rri[tx][ty]*wet;
			//}
		}
	}




}

__global__ void twodimbndnoix(int nx, int ny, DECNUM eps, DECNUM * hh, DECNUM * F)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	unsigned int xminus = mminus(ix, nx);
	unsigned int xplus = pplus(ix, nx);
	unsigned int yminus = mminus(iy, ny);
	unsigned int yplus = pplus(iy, ny);

	__shared__ DECNUM hhi[16][16];
	__shared__ DECNUM Fi[16][16];
	__shared__ DECNUM Ft[16][16];
	__shared__ DECNUM Fb[16][16];
	__shared__ DECNUM Fr[16][16];



	if (ix < nx && iy < ny)
	{

		hhi[tx][ty] = hh[i];
		DECNUM wet = 0.0f;

		if (hhi[tx][ty] > eps)
		{
			wet = 1.0f;
		}


		Fi[tx][ty] = F[i];
		Ft[tx][ty] = F[ix + yplus*nx];
		Fb[tx][ty] = F[ix + yminus*nx];
		Fr[tx][ty] = F[xplus + iy*nx];

		//F[i]=Fi[tx][ty]*wet;
		if (iy == 0)
		{
			F[i] = Ft[tx][ty] * wet;
		}
		if (iy == ny - 1)
		{
			F[i] = Fb[tx][ty] * wet;
		}
		//if (ix==0)
		//{
		//	F[i]=Fr[tx][ty]*wet;
		//}
	}

}
__global__ void twodimbnd(int nx, int ny, DECNUM eps, DECNUM * hh, DECNUM * F)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		__shared__ DECNUM hhi[16][16];
		__shared__ DECNUM Fi[16][16];
		__shared__ DECNUM Ft[16][16];
		__shared__ DECNUM Fb[16][16];
		__shared__ DECNUM Fr[16][16];





		hhi[tx][ty] = hh[i];
		DECNUM wet = 0.0f;

		if (hhi[tx][ty] > eps)
		{
			wet = 1.0f;
		}


		Fi[tx][ty] = F[i];
		Ft[tx][ty] = F[ix + yplus*nx];
		Fb[tx][ty] = F[ix + yminus*nx];
		Fr[tx][ty] = F[xplus + iy*nx];

		//F[i]=Fi[tx][ty]*wet;
		if (iy == 0)
		{
			F[i] = Ft[tx][ty] * wet;
		}
		if (iy == ny - 1)
		{
			F[i] = Fb[tx][ty] * wet;
		}
		if (ix == 0)
		{
			F[i] = Fr[tx][ty] * wet;
		}
	}

}

__global__ void dissipation(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM eps, DECNUM dt, DECNUM g, DECNUM beta, DECNUM * wci, DECNUM *hh, DECNUM *ee, DECNUM *D, DECNUM *E, DECNUM *rr, DECNUM *c, DECNUM *cxsth, DECNUM *sxnth, DECNUM * uu, DECNUM * vv, DECNUM *DR, DECNUM *R)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	DECNUM dd;
	DECNUM RR = 0.0f;
	DECNUM drr = 0.0f;
	DECNUM cx, cy;
	//DECNUM wci=0.0f;
	DECNUM rrd;
	DECNUM eet, rrt;

	__shared__ DECNUM cc[16][16];
	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM vvi[16][16];


	if (ix < nx && iy < ny)
	{
		cc[tx][ty] = c[i];
		uui[tx][ty] = uu[i] * wci[i];
		vvi[tx][ty] = vv[i] * wci[i];
		__syncthreads();

		for (int itheta = 0; itheta<ntheta; itheta++)
		{
			cx = cc[tx][ty] * cxsth[itheta] + uui[tx][ty];
			cy = cc[tx][ty] * sxnth[itheta] + vvi[tx][ty];

			eet = ee[i + itheta*nx*ny];
			dd = eet*D[i] / max(E[i], 0.0000001f);




			if (hh[i]>eps)
			{
				rrt = rr[i + itheta*nx*ny];
				rrd = 2.0f*g*beta*rrt / sqrtf(cx*cx + cy*cy);
				eet = max(eet - dt*dd, 0.0f);
				rrt = max(rrt + dt*(dd - rrd), 0.0f);
				drr = drr + rrd;
			}
			else
			{
				eet = 0.0f;
				rrt = 0.0f;
				drr = 0.0f;

			}
			RR = RR + rrt;
			ee[i + itheta*nx*ny] = eet;
			rr[i + itheta*nx*ny] = rrt;
		}
		R[i] = RR*dtheta;
		DR[i] = drr*dtheta;
	}
}

__global__ void meandir(int nx, int ny, int ntheta, DECNUM rho, DECNUM g, DECNUM dtheta, DECNUM * ee, DECNUM * thet, DECNUM * thetamean, DECNUM * E, DECNUM * H)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;


	if (ix < nx && iy < ny)
	{
		DECNUM sumethet = 0;

		DECNUM sume = 0;
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			sume = sume + ee[i + itheta*nx*ny];
			sumethet = sumethet + ee[i + itheta*nx*ny] * thet[itheta];
		}
		sume = max(sume, 0.00001f);
		__syncthreads;
		thetamean[i] = (sumethet / ntheta) / (sume / ntheta);
		E[i] = sume*dtheta;

		H[i] = sqrtf(sume*dtheta / (rho*g / 8));//sqrt(E[i]/(1/8*rho*g));
	}
}

__global__ void meanSingledir(int nx, int ny, int ntheta, DECNUM rho, DECNUM g, DECNUM dtheta, DECNUM theta0, DECNUM * c, DECNUM * ee,  DECNUM * thetamean, DECNUM * E, DECNUM * H)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;


	if (ix < nx && iy < ny)
	{
		DECNUM sumethet = 0;
		DECNUM ttm = 0;
		DECNUM sume = 0;
		DECNUM c1;

		c1 = c[0];


		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			sume = sume + ee[i + itheta*nx*ny];
		//	sumethet = sumethet + ee[i + itheta*nx*ny] * thet[itheta];
		}
		sume = max(sume, 0.00001f);
		__syncthreads;
		//ttm = (sumethet / ntheta) / (sume / ntheta);
		//theta0 = (1.5f * pi) - Dp * atanf(1.0f) / 45.0f; // already pbeen pre-calculated
		c1 = c[0];
		thetamean[i] = asin(max(-1.0f,min(1.0f,sinf(theta0)*c[i]/c1)));
		E[i] = sume*dtheta;

		H[i] = sqrtf(sume*dtheta / (rho*g / 8.0f));//sqrt(E[i]/(1/8*rho*g));
				
	}
}

__global__ void radstress(int nx, int ny, int ntheta, DECNUM dx, DECNUM dtheta, DECNUM * ee, DECNUM *rr, DECNUM * cxsth, DECNUM * sxnth, DECNUM * cg, DECNUM * c, DECNUM * Sxx, DECNUM * Sxy, DECNUM * Syy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;



	//__shared__ DECNUM n[16][16];

	if (ix < nx && iy < ny)
	{

		DECNUM n = cg[i] / c[i];
		
		DECNUM sume = 0.0f;
		DECNUM tmpxx = 0.0f;
		DECNUM tmpxy = 0.0f;
		DECNUM tmpyy = 0.0f;
		DECNUM rolxx = 0.0f;
		DECNUM rolxy = 0.0f;
		DECNUM rolyy = 0.0f;

		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			sume = sume + ee[i + itheta*nx*ny];
			tmpxx = tmpxx + ((1.0f + cxsth[itheta] * cxsth[itheta])*ee[i + itheta*nx*ny]);
			tmpyy = tmpyy + ((1.0f + sxnth[itheta] * sxnth[itheta])*ee[i + itheta*nx*ny]);
			tmpxy = tmpxy + (sxnth[itheta] * cxsth[itheta] * ee[i + itheta*nx*ny]);
			rolxx = rolxx + ((cxsth[itheta] * cxsth[itheta])*rr[i + itheta*nx*ny]);
			rolyy = rolyy + ((sxnth[itheta] * sxnth[itheta])*rr[i + itheta*nx*ny]);
			rolxy = rolxy + (sxnth[itheta] * cxsth[itheta] * rr[i + itheta*nx*ny]);

		}


		Sxx[i] = (n * tmpxx - 0.5f*sume)*dtheta + rolxx*dtheta;
		Syy[i] = (n * tmpyy - 0.5f*sume)*dtheta + rolyy*dtheta;
		Sxy[i] = n * tmpxy*dtheta + rolxy*dtheta;
	}
}

__global__ void wavforce(int nx, int ny, int ntheta, DECNUM dx, DECNUM dtheta, DECNUM * Sxx, DECNUM * Sxy, DECNUM * Syy, DECNUM * Fx, DECNUM * Fy, DECNUM * hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;



	//DECNUM hmin=0.000002f;

	__shared__ DECNUM  Sxxi[16][16];
	__shared__ DECNUM  Sxxr[16][16];
	__shared__ DECNUM  Sxyi[16][16];
	__shared__ DECNUM  Sxyt[16][16];
	__shared__ DECNUM  Sxytr[16][16];
	__shared__ DECNUM  Sxyb[16][16];
	__shared__ DECNUM  Sxybr[16][16];

	__shared__ DECNUM  Syyi[16][16];
	__shared__ DECNUM  Syyt[16][16];
	__shared__ DECNUM  Sxyr[16][16];
	__shared__ DECNUM  Sxyl[16][16];
	__shared__ DECNUM  Sxytl[16][16];

	__shared__ DECNUM FFx[16][16];
	__shared__ DECNUM FFy[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		Sxxi[tx][ty] = Sxx[i];
		Sxxr[tx][ty] = Sxx[xplus + iy*nx];
		Sxyi[tx][ty] = Sxy[i];
		Sxyt[tx][ty] = Sxy[ix + yplus*nx];
		Sxytr[tx][ty] = Sxy[xplus + yplus*nx];
		Sxyb[tx][ty] = Sxy[ix + yminus*nx];
		Sxybr[tx][ty] = Sxy[xplus + yminus*nx];

		Syyi[tx][ty] = Syy[i];
		Syyt[tx][ty] = Syy[ix + yplus*nx];
		Sxyr[tx][ty] = Sxy[xplus + iy*nx];
		Sxyl[tx][ty] = Sxy[xminus + iy*nx];
		Sxytl[tx][ty] = Sxy[xminus + yplus*nx];

		FFx[tx][ty] = 0.0f;
		FFy[tx][ty] = 0.0f;


		__syncthreads;


		//if(hh[i]>hmin)
		//{
		FFx[tx][ty] = -1.0f * (Sxxr[tx][ty] - Sxxi[tx][ty]) / dx - (Sxyt[tx][ty] + Sxytr[tx][ty] - Sxyb[tx][ty] - Sxybr[tx][ty]) / (4.0f*dx);

		FFy[tx][ty] = -1.0f * (Syyt[tx][ty] - Syyi[tx][ty]) / dx - (Sxyr[tx][ty] + Sxytr[tx][ty] - Sxyl[tx][ty] - Sxytl[tx][ty]) / (4.0f*dx);
		//}


		Fx[i] = FFx[tx][ty];

		Fy[i] = FFy[tx][ty];
	}


}

__global__ void breakerdelay(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM g, DECNUM rho, DECNUM Trep, DECNUM eps, DECNUM * urms, DECNUM *ust, DECNUM *H, DECNUM *E, DECNUM *c, DECNUM *k, DECNUM *hh, DECNUM *rr)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix < nx && iy < ny)
	{
		DECNUM hmin = 0.2;
		DECNUM delta = 0.0f;
		DECNUM Hi = H[i];
		DECNUM ki = k[i];
		DECNUM hhi = max(hh[i], Hi*delta);
		DECNUM ci = max(c[i], sqrt(hmin*g));
		DECNUM R;
		DECNUM ustw, uwf, vwf, ustr, usd, uorb;

		R = rr[i];





		uorb = pi*Hi / Trep / sinhf(min(max(ki, 0.01f)*hhi, 10.0f));

		urms[i] = uorb / (sqrtf(2.0f));
		ustw = E[i] / ci / rho / hhi;
		//uwf = ustw*cosf(tm[i]);
		//vwf = ustw*sinf(tm[i]);
		// roller contribution

		ustr = 2.0f*R / ci / rho / hhi;
		// introduce breaker delay
		//usd=ustr; //I don't like breaker delay


		ust[i] = ustw + ustr;
	}
}

__global__ void calcwci(int nx, int ny, DECNUM wci, DECNUM hwci, DECNUM * hh, DECNUM *wcig)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;


	if (ix < nx && iy < ny)
	{

		wcig[i] = wci*min(hh[i] / hwci, 1.0f);
	}

}



