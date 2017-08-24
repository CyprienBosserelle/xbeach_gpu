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

__device__ float minmod(float delta1, float delta2)
{
	float minmod = 0.0f;

	if (delta1*delta2 > 0.0f)
	{


		if (delta1 > 0.0f)
		{
			minmod = min(delta1, delta2);

		}
		if (delta1 < 0.0f)
		{
			minmod = max(delta1, delta2);
		}
	}

	return(minmod);
}

__global__ void ubnd(int nx, int ny, DECNUM dx, DECNUM dt, DECNUM g, DECNUM rho, DECNUM totaltime, DECNUM wavbndtime, DECNUM rt, DECNUM slbndtime, DECNUM rtsl, DECNUM zsbndold, DECNUM zsbndnew, DECNUM Trep, DECNUM * qbndold, DECNUM * qbndnew, DECNUM *zs, DECNUM * uu, DECNUM * vv, DECNUM *vu, DECNUM * umean, DECNUM * vmean, DECNUM * zb, DECNUM * cg, DECNUM * hum, DECNUM * zo, DECNUM *Fx, DECNUM *hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		DECNUM ui, vi, thetai, vert;
		DECNUM beta, betar, betat, betab, bnp1, bn;
		DECNUM ht, htr;
		DECNUM theta0 = 0.0f;
		DECNUM alpha2 = -1.0f*theta0;
		DECNUM epsi = 0.005; //Not used!
		DECNUM ur, uumean, vvmean, urr, alphanew;
		DECNUM dbetadx, dbetady, dvudy, dhdx;
		DECNUM qx, qy, zsbnd;
		DECNUM order = 2.0f;
		DECNUM ccg = cg[i];
		DECNUM cats = 4; // number of wave period to average the current from
		DECNUM factime = 0.0f;// 1.0f/cats/Trep*dt;
		DECNUM taper = min(totaltime / 100.0f, 1.0f);

		if (ix == 0)
		{

			qx = (qbndold[iy] + (totaltime - wavbndtime + rt)*(qbndnew[iy] - qbndold[iy]) / rt)*taper;
			qy = (qbndold[iy + ny] + (totaltime - wavbndtime + rt)*(qbndnew[iy + ny] - qbndold[iy + ny]) / rt)*taper;
			zsbnd = zsbndold + (totaltime - rtsl)*(zsbndnew - zsbndold) / (slbndtime - rtsl);

			ht = zsbnd + zb[i];
			htr = zsbnd + zb[xplus + iy*nx];
			ui = qx / ht;
			vi = qy / ht;
			beta = uu[i] - 2.0f*sqrt(g*hum[i]);
			betar = uu[xplus + iy*nx] - 2.0f*sqrtf(g*hum[xplus + iy*nx]);
			betat = uu[ix + yplus*nx] - 2.0f*sqrtf(g*hum[ix + yplus*nx]);
			betab = uu[ix + yminus*nx] - 2.0f*sqrtf(g*hum[ix + yminus*nx]);

			dvudy = (vu[ix + (yminus)*nx] - vu[ix + (yplus)*nx]) / (2.0f*dx);
			dbetadx = (betar - beta) / dx;
			dbetady = (betat - betab) / (2.0f*dx);

			dhdx = (htr - ht) / dx;

			bn = -1.0f*(uu[i] - sqrt(g*hum[i]))*dbetadx - vu[i] * dbetady + sqrtf(g*hum[i])*dvudy + 1 / rho*Fx[i] / hum[i] - zo[i] * sqrtf(uu[i] * uu[i] + vu[i] * vu[i])*uu[i] / hum[i] + g*dhdx;
			bnp1 = beta + bn*dt;

			//WARNING this should be very inefficient. Need to find a better way. possibly inside another kernel
			// not neededd when epsi ==0.0...or factime==0.0
			DECNUM uumm = 0.0f;
			DECNUM vvmm = 0.0f;
			/*for (int jj=0; jj<ny; jj++)
			{
			uumm=uumm+uu[ix+jj*nx];
			vvmm=vvmm+vv[ix+jj*nx];
			}*/



			uumean = factime*uumm + umean[iy] * (1 - factime);
			vvmean = factime*vvmm + vmean[iy] * (1 - factime);
			umean[iy] = uumean;
			vmean[iy] = vvmean;




			thetai = atanf(vi / (ui + 0.0000001f));

			vert = vu[i] - vvmean - vi;

			urr = (bnp1 - uumean + 2.0f*sqrtf(g*0.5f*(ht + htr)) - ui*(ccg*(cosf(thetai)) - sqrtf(g*0.5f*(ht + htr))) / (ccg*cosf(thetai)));

			for (int jj = 0; jj < 50; jj++)
			{
				ur = cosf(alpha2) / (cosf(alpha2) + 1.0f)*urr;
				/*if(ur==0.0f)
				{
				ur=0.0000001f;
				}*/
				alphanew = atanf(vert / (ur + 0.0000001f));
				if (alphanew > pi*0.5f)
				{
					alphanew = alphanew - pi;
				}
				if (alphanew <= -0.5f*pi)
				{
					alphanew = alphanew + pi;
				}
				if (abs(alphanew - alpha2) < 0.001f)
				{
					break;
				}
				alpha2 = alphanew;
			}


			//
			uu[i] = uumean;//(order - 1.0f)*ui + ur + uumean;//2.0f*ui-(sqrtf(g/(zs[i]+zb[i]))*(zs[i]-zsbnd));;//
			zs[i] = 1.5f*((bnp1 - uu[i])*(bnp1 - uu[i]) / (4.0f*g) - 0.5f*(zb[i] + zb[xplus + iy*nx])) - 0.5f*((betar - uu[xplus + iy*nx])*(betar - uu[xplus + iy*nx]) / (4.0f*g) - 0.5f*(zb[xplus + iy*nx] + zb[xplus2 + iy*nx]));
			////
			//zsbnd+qx/(dx*dx)*dt;//

			hh[i] = zs[i] + zb[i];
			vv[i] = vv[xplus + iy*nx];
		}




		//__syncthreads;
	}


}
__global__ void ubnd1D(int nx, int ny, DECNUM dx, DECNUM dt, DECNUM g, DECNUM rho, DECNUM totaltime, DECNUM time, DECNUM timenext, DECNUM zsbnd, DECNUM Trep, DECNUM * qbndold, DECNUM * qbndnew, DECNUM *zs, DECNUM * uu, DECNUM * vv, DECNUM *vu, DECNUM * umean, DECNUM * vmean, DECNUM * zb, DECNUM * cg, DECNUM * hum, DECNUM * zo, DECNUM *Fx, DECNUM *hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix == 0)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		DECNUM ui, vi, thetai, vert;
		
		DECNUM ht, htr;
		DECNUM epsi = 0.005; //Not used!
		DECNUM ur, uumean, vvmean, urr, alphanew;
		
		DECNUM qx, qy;
		
		DECNUM cats = 4.0f; // number of wave period to average the current from
		DECNUM factime = 1.0f/cats/Trep*dt;
		DECNUM taper = min(totaltime / 100.0f, 1.0f);
		DECNUM zsplus = zs[xplus + iy*nx];
		
		//prev + (time) / (timenext)*(next - prev);
		qx = (qbndold[iy] + time/ timenext*(qbndnew[iy] - qbndold[iy]) )*taper;
		qy = (qbndold[iy + ny] + time / timenext*(qbndnew[iy + ny] - qbndold[iy + ny]) )*taper;

			//zsbnd = zsbndold + (totaltime - rtsl)*(zsbndnew - zsbndold) / (slbndtime - rtsl);

			ht = zsbnd + zb[i];
			htr = zsbnd + zb[xplus + iy*nx];
			ui = qx / ht;
			vi = qy / ht;
			ur = -1.0f*sqrtf(g/hh[i])*(zsplus - zsbnd);
			
			float uuavg = 0.0f;
			for (int n = 0; n < ny; n++)
			{
				uuavg = uuavg + uu[ix + n*nx];
			}
			uuavg = uuavg / ny;
			

			//Tidetype=velocity
			uumean = factime*uuavg + umean[iy] * (1 - factime);
			//vvmean = factime*vvmm + vmean[iy] * (1 - factime);
			umean[iy] = uumean;
			//vmean[iy] = vvmean;



			// Freewave==1
			uu[i] =  (2.0f)*ui + ur + uumean;//2.0f*ui-(sqrtf(g/(zs[i]+zb[i]))*(zs[i]-zsbnd));;//
			//zs[i] = 1.5f*((bnp1 - uu[i])*(bnp1 - uu[i]) / (4.0f*g) - 0.5f*(zb[i] + zb[xplus + iy*nx])) - 0.5f*((betar - uu[xplus + iy*nx])*(betar - uu[xplus + iy*nx]) / (4.0f*g) - 0.5f*(zb[xplus + iy*nx] + zb[xplus2 + iy*nx]));
			////
			//zsbnd+qx/(dx*dx)*dt;//
			zs[i] =  zsplus;
			//hh[i] = zsplus + zb[i];
			
			vv[i] = vv[xplus + iy*nx];
						
		




		
	}


}


// I doubt this should be on the GPU but although it is slow it may be ok for now

__global__ void discharge_bnd_h(int nx, int ny, DECNUM dx, DECNUM eps,DECNUM qnow, int istart, int jstart, int iend,int jend, DECNUM *hu, DECNUM *hv, DECNUM *qx, DECNUM *qy, DECNUM * uu, DECNUM *vv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix >= istart && ix <= iend &&iy >= jstart && iy <= jend)
	{
		//discharge should run along the x or y axis
		//TODO modify the algorithm to alow diagonal discharges
		if (istart == iend)//discharge on uu along a vertical line
		{
			float A = 0.0;

			for (int k = jstart; k <= jend; k++)
			{
				A = A+(cbrtf(hu[ix + k*ny])*dx);
			}

			float cst = qnow / max(A, eps);
			uu[i] = cst*sqrt(hu[i]);
			qx[i] = cst*cbrtf(hu[i]);
		}

		if (jstart == jend)//discharge on uu along a vertical line
		{
			float A = 0.0;

			for (int k = istart; k <= iend; k++)
			{
				A = A + (cbrtf(hv[k + iy*ny])*dx);
			}

			float cst = qnow / max(A, eps);
			vv[i] = cst*sqrt(hv[i]);
			qy[i] = cst*cbrtf(hv[i]);
		}




	}



}

__global__ void discharge_bnd_v(int nx, int ny, DECNUM dx, DECNUM eps, DECNUM dt, DECNUM qnow, int istart, int jstart, int iend, int jend, DECNUM * zs, DECNUM *hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix >= istart && ix <= iend &&iy >= jstart && iy <= jend)
	{
		float A = (iend - istart + 1)*dx*(jend - jstart + 1)*dx;

		float dzsdt = qnow*dt / A;
		zs[i] = zs[i] + dzsdt;
		// Do hh[i] too although Im not sure it is worth it
		hh[i] = hh[i] + dzsdt;

	}
}



__global__ void ubnd1Dnowaves(int nx, int ny, DECNUM dx, DECNUM dt, DECNUM g, DECNUM rho, DECNUM totaltime, DECNUM time, DECNUM timenext, DECNUM zsbnd,  DECNUM *zs, DECNUM * uu, DECNUM * vv, DECNUM *vu, DECNUM * umean, DECNUM * vmean, DECNUM * zb, DECNUM * hum, DECNUM * zo,  DECNUM *hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix == 0)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		DECNUM ui, vi, thetai, vert;

		DECNUM ht, htr;
		DECNUM epsi = 0.005; //Not used!
		DECNUM ur, uumean, vvmean, urr, alphanew;

		DECNUM qx, qy;

		DECNUM cats = 4.0f; // number of wave period to average the current from
		DECNUM factime = 1.0f;// 1.0f / cats / Trep*dt;
		DECNUM taper = min(totaltime / 100.0f, 1.0f);
		DECNUM zsplus = zs[xplus + iy*nx];

		//prev + (time) / (timenext)*(next - prev);
		qx = 0.0f;
		qy = 0.0f;

		//zsbnd = zsbndold + (totaltime - rtsl)*(zsbndnew - zsbndold) / (slbndtime - rtsl);

		ht = zsbnd + zb[i];
		htr = zsbnd + zb[xplus + iy*nx];
		ui = qx / ht;
		vi = qy / ht;
		ur = -1.0f*sqrtf(g / hh[i])*(zsplus - zsbnd);

		float uuavg = 0.0f;
		for (int n = 0; n < ny; n++)
		{
			uuavg = uuavg + uu[ix + n*nx];
		}
		uuavg = uuavg / ny;


		//Tidetype=velocity
		uumean = factime*uuavg + umean[iy] * (1 - factime);
		//vvmean = factime*vvmm + vmean[iy] * (1 - factime);
		umean[iy] = uumean;
		//vmean[iy] = vvmean;



		// Freewave==1
		uu[i] =  ur + uumean;//2.0f*ui-(sqrtf(g/(zs[i]+zb[i]))*(zs[i]-zsbnd));;//
										//zs[i] = 1.5f*((bnp1 - uu[i])*(bnp1 - uu[i]) / (4.0f*g) - 0.5f*(zb[i] + zb[xplus + iy*nx])) - 0.5f*((betar - uu[xplus + iy*nx])*(betar - uu[xplus + iy*nx]) / (4.0f*g) - 0.5f*(zb[xplus + iy*nx] + zb[xplus2 + iy*nx]));
										////
										//zsbnd+qx/(dx*dx)*dt;//
		zs[i] = zsplus;
		//hh[i] = zsplus + zb[i];

		vv[i] = vv[xplus + iy*nx];







	}


}


__global__ void wlevslopes(int nx, int ny, DECNUM dx, DECNUM eps, DECNUM *zs, DECNUM * dzsdx, DECNUM *dzsdy, DECNUM*hh)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		__shared__ DECNUM zsi[16][16];
		__shared__ DECNUM zsr[16][16];
		__shared__ DECNUM zst[16][16];
		//int whi;
		//int whr;
		//int wht;


		// Should dzsdx be ==0 near dry cells?

		//whi=0;
		//whr=0;
		//wht=0;


		/*	if (hh[i]>eps)
			{
			whi=1;
			}
			if(hh[xplus+iy*nx]>eps)
			{
			whr=1;
			}
			if(hh[ix+yplus*nx]>eps)
			{
			wht=1;

			}
			*/
		zsi[tx][ty] = zs[i];
		zsr[tx][ty] = zs[xplus + iy*nx];
		zst[tx][ty] = zs[ix + yplus*nx];


		//dzsdx[i]=(zs[ix+1+iy*nx]-zs[ix-1+iy*nx])/(2*dx);
		dzsdx[i] = (zsr[tx][ty] - zsi[tx][ty]) / dx;//*whi*whr;
		dzsdy[i] = (zst[tx][ty] - zsi[tx][ty]) / dx;//*whi*wht;
		
	}


}
__global__ void calcuvvu(int nx, int ny, DECNUM dx, DECNUM *uu, DECNUM *vv, DECNUM *vu, DECNUM *uv, DECNUM * ust, DECNUM *thetamean, DECNUM *ueu_g, DECNUM *vev_g, DECNUM *vmageu, DECNUM *vmagev, DECNUM * uuold, DECNUM * vvold, DECNUM * zs, DECNUM *zsold, int* wetu, int* wetv)
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

		DECNUM vsu, usu, ueu, veu, usv, vsv, vev, uev;

		__shared__ DECNUM usti[16][16];
		__shared__ DECNUM ustr[16][16];
		__shared__ DECNUM ustt[16][16];

		__shared__ DECNUM tmeani[16][16];
		__shared__ DECNUM tmeanr[16][16];
		__shared__ DECNUM tmeant[16][16];
		__shared__ int wetui[16][16];
		__shared__ int wetvi[16][16];

		usti[tx][ty] = ust[i];
		tmeani[tx][ty] = thetamean[i];
		ustr[tx][ty] = ust[xplus + iy*nx];
		tmeanr[tx][ty] = thetamean[xplus + iy*nx];
		ustt[tx][ty] = ust[ix + yplus*nx];
		tmeant[tx][ty] = thetamean[ix + yplus*nx];
		wetui[tx][ty] = wetu[i];
		wetvi[tx][ty] = wetv[i];


		// Uold and Vold used in 2nd order flow correction
		uuold[i] = uu[i];
		vvold[i] = vv[i];
		zsold[i] = zs[i];
		// V-velocities at u-points

		vu[i] = 0.25f*(vv[ix + yminus*nx] + vv[ix + iy*nx] + vv[xplus + yminus*nx] + vv[xplus + iy*nx])*wetui[tx][ty];

		// U-velocities at v-points
		uv[i] = 0.25f*(uu[xminus + iy*nx] + uu[ix + iy*nx] + uu[xminus + yplus*nx] + uu[ix + yplus*nx])*wetvi[tx][ty];



		//Calculate V-stokes at u points
		vsu = 0.5f*(usti[tx][ty] * sinf(tmeani[tx][ty]) + ustr[tx][ty] * sinf(tmeanr[tx][ty]))*wetui[tx][ty];
		//Calculate U-stokes at u points
		usu = 0.5f*(usti[tx][ty] * cosf(tmeani[tx][ty]) + ustr[tx][ty] * cosf(tmeanr[tx][ty]))*wetui[tx][ty];
		//Calculate U-euler at u points
		ueu = uu[i] - usu;
		//Calculate V-euler at u points
		veu = vu[i] - vsu;
		vmageu[i] = sqrtf(ueu*ueu + veu*veu);
		ueu_g[i] = ueu;


		usv = 0.5f*(usti[tx][ty] * cosf(tmeani[tx][ty]) + ustt[tx][ty] * cosf(tmeant[tx][ty]))*wetvi[tx][ty];
		vsv = 0.5f*(usti[tx][ty] * sinf(tmeani[tx][ty]) + ustt[tx][ty] * sinf(tmeant[tx][ty]))*wetvi[tx][ty];
		vev = vv[i] - vsv;
		uev = uv[i] - usv;
		vmagev[i] = sqrtf(uev*uev + vev*vev);
		vev_g[i] = vev;

	}


}




__global__ void udepthmomcont(int nx, int ny, DECNUM dx, DECNUM eps, DECNUM ummn, int* wetu, DECNUM * zs, DECNUM * uu, DECNUM * hh, DECNUM *hum, DECNUM *hu, DECNUM * zb)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		DECNUM humi;
		DECNUM hui;

		__shared__ DECNUM  hhi[16][16];
		__shared__ DECNUM  hhip[16][16];
		//__shared__ DECNUM  hhjp[4][4];

		hui = hu[i];
		hhi[tx][ty] = hh[i];
		hhip[tx][ty] = hh[xplus + iy*nx];
		//hhjp[tx][ty]=hh[ix+(iy+1)*nx];
		__syncthreads;

		//Water depth in u-points do momentum equation: mean
		humi = 0.5f*(hhi[tx][ty] + hhip[tx][ty]);
		// Water depth in u-points do continuity equation: upwind





		__syncthreads;
		if (hui > eps && humi > eps)
		{
			wetu[i] = 1;
		}
		else
		{
			wetu[i] = 0;
		}


		hum[i] = max(humi, eps);
	}
}

__global__ void vdepthmomcont(int nx, int ny, DECNUM dx, DECNUM eps, DECNUM ummn, int* wetv, DECNUM * zs, DECNUM * vv, DECNUM * hh, DECNUM *hvm, DECNUM *hv, DECNUM * zb)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		DECNUM hvmi, hvi;

		__shared__ DECNUM  hhi[16][16];
		//__shared__ DECNUM  hhip[4][4];
		__shared__ DECNUM  hhjp[16][16];


		hhi[tx][ty] = hh[i];
		//hhip[tx][ty]=hh[ix+1+iy*nx];
		hhjp[tx][ty] = hh[ix + yplus*nx];
		__syncthreads;


		//Water depth in u-points do momentum equation: mean
		//hvmi=max(0.5f*(hh[i]+hh[ix+(min(iy,ny-2)+1)*nx]),eps);
		// Water depth in u-points do continuity equation: upwind

		//hvi=0.5f*(hhjp[tx][ty]-hhi[tx][ty])+hhi[tx][ty];
		hvmi = 0.5f*(hhi[tx][ty] + hhjp[tx][ty]);
		//hvm(i,j)=max(.5d0*(hh(i,j)+hh(i,min(ny,j)+1)),par%eps)  


		hvi = hv[i];


		if (hvi > eps && hvmi > eps)
		{
			wetv[i] = 1;
		}
		else
		{
			wetv[i] = 0;
		}
		hvm[i] = max(hvmi, eps);
	}
}

__global__ void depthhu(int nx, int ny, DECNUM dx, DECNUM ummn, DECNUM eps, DECNUM *hh, DECNUM * uu, DECNUM * hu, DECNUM *zs, DECNUM *zb)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		DECNUM hui = 0.0f;

		__shared__ DECNUM  hhi[16][16];
		__shared__ DECNUM  hhip[16][16];
		hhi[tx][ty] = hh[i];
		hhip[tx][ty] = hh[xplus + iy*nx];

		if (uu[i]>ummn)
		{
			//hui=hhi[tx][ty];
			hui = zs[i] - max(-1.0f*zb[i], -1.0f*zb[xplus + iy*nx]);
		}
		else
		{
			if (uu[i] < -1.0f*ummn)
			{
				//hui=hhip[tx][ty];
				hui = zs[xplus + iy*nx] - max(-1.0f*zb[i], -1.0f*zb[xplus + iy*nx]);
			}
			else
			{
				hui = max(max(zs[i], zs[xplus + iy*nx]) - max(-1.0f*zb[i], -1.0f*zb[xplus + iy*nx]), eps);
			}

		}
		//hui=0.5f*(hhip[tx][ty]+hhi[tx][ty]);
		hui = max(hui, 0.0f);
		hu[i] = hui;
	}

}

__global__ void depthhv(int nx, int ny, DECNUM dx, DECNUM ummn, DECNUM eps, DECNUM *hh, DECNUM * vv, DECNUM * hv, DECNUM *zs, DECNUM *zb)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		DECNUM hvi = 0.0f;

		__shared__ DECNUM  hhi[16][16];
		//__shared__ DECNUM  hhip[4][4];
		__shared__ DECNUM  hhjp[16][16];


		hhi[tx][ty] = hh[i];
		//hhip[tx][ty]=hh[ix+1+iy*nx];
		hhjp[tx][ty] = hh[ix + yplus*nx];
		__syncthreads;
		if (vv[i]>ummn)
		{
			//hvi=hhi[tx][ty];
			hvi = zs[i] - max(-1.0f*zb[i], -1.0f*zb[ix + yplus*nx]);
		}
		else
		{
			if (vv[i] < -1 * ummn)
			{
				//hvi=hhjp[tx][ty];
				hvi = zs[ix + yplus*nx] - max(-1.0f*zb[i], -1.0f*zb[ix + yplus*nx]);
			}
			else
			{
				hvi = max(max(zs[i], zs[ix + yplus*nx]) - max(-1.0f*zb[i], -1.0f*zb[ix + yplus*nx]), eps);
				//hv[i]=hvm[i];
			}
		}
		hvi = max(hvi, 0.0f);

		hv[i] = hvi;
	}
}

__global__ void ududx_adv(int nx, int ny, DECNUM dx, DECNUM * hu, DECNUM * hum, DECNUM * uu, DECNUM * ududx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;
	DECNUM qin, uududx, du;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uur[16][16];
	__shared__ DECNUM uul[16][16];


	__shared__ DECNUM hui[16][16];




	__shared__ DECNUM humi[16][16];
	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		uui[tx][ty] = uu[i];
		uur[tx][ty] = uu[xplus + iy*nx];
		uul[tx][ty] = uu[xminus + iy*nx];


		hui[tx][ty] = hu[i];
		humi[tx][ty] = hum[i];



		uududx = 0.0f;
		qin = 0.5f*(hui[tx][ty] * uui[tx][ty] + hu[xminus + iy*nx] * uul[tx][ty]);
		du = uui[tx][ty] - uul[tx][ty];
		//ududx
		if (qin > 0.0f)
		{
			if (du > 0.5f)
			{
				uududx = uududx + 0.5f*(uul[tx][ty] + uui[tx][ty])* (uui[tx][ty] - uul[tx][ty]) / dx;
			}
			else
			{
				uududx = uududx + qin / humi[tx][ty] * (uui[tx][ty] - uul[tx][ty]) / dx;
			}
			
		}
		qin = -0.5f*(hui[tx][ty] * uui[tx][ty] + hu[xplus + iy*nx] * uur[tx][ty]);
		du = uur[tx][ty] - uui[tx][ty];
		if (qin > 0.0f)
		{
			if (du > 0.5f)
			{
				uududx = uududx - 0.5f*(uur[tx][ty] + uui[tx][ty])* (uui[tx][ty] - uur[tx][ty]) / dx;
			}
			else
			{
				uududx = uududx + qin / humi[tx][ty] * (uui[tx][ty] - uur[tx][ty]) / dx;
			}
			
		}
		ududx[i] = uududx;
	}
}


__global__ void ududx_adv2(int nx, int ny, DECNUM dx, DECNUM * hu, DECNUM * hum, DECNUM * uu, DECNUM * ududx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;



	DECNUM qin, uududx;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uur[16][16];
	__shared__ DECNUM uur2[16][16];
	__shared__ DECNUM uul[16][16];
	__shared__ DECNUM uul2[16][16];

	__shared__ DECNUM hui[16][16];
	__shared__ DECNUM hur[16][16];
	__shared__ DECNUM hul[16][16];
	__shared__ DECNUM humi[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		uui[tx][ty] = uu[i];
		uur[tx][ty] = uu[xplus + iy*nx];
		uur2[tx][ty] = uu[xplus2 + iy*nx];
		uul[tx][ty] = uu[xminus + iy*nx];
		uul2[tx][ty] = uu[xminus2 + iy*nx];


		hui[tx][ty] = hu[i];
		hur[tx][ty] = hu[xplus + iy*nx];
		hul[tx][ty] = hu[xminus + iy*nx];
		humi[tx][ty] = hum[i];



		uududx = 0.0f;
		qin = 0.5f*(hui[tx][ty] * uui[tx][ty] + hul[tx][ty] * uul[tx][ty]);
		//ududx
		if (qin > 0.0f)
		{
			uududx = uududx + qin / humi[tx][ty] * (3 * uui[tx][ty] - 4 * uul[tx][ty] + uul2[tx][ty]) / (2 * dx);
		}
		qin = -0.5f*(hui[tx][ty] * uui[tx][ty] + hur[tx][ty] * uur[tx][ty]);
		if (qin > 0.0f)
		{
			uududx = uududx + qin / humi[tx][ty] * (3 * uui[tx][ty] - 4 * uur[tx][ty] + uur2[tx][ty]) / (2 * dx);
		}
		ududx[i] = uududx;
	}
}



__global__ void vdudy_adv(int nx, int ny, DECNUM dx, DECNUM * hv, DECNUM * hum, DECNUM * uu, DECNUM *vv, DECNUM * vdudy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	DECNUM qin, dv, vvdudy;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uut[16][16];
	__shared__ DECNUM uub[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		uui[tx][ty] = uu[i];
		uut[tx][ty] = uu[ix + yplus*nx];
		uub[tx][ty] = uu[ix + yminus*nx];

		vvdudy = 0.0f;

		qin = 0.5f*(vv[ix + yminus*nx] * hv[ix + yminus*nx] + vv[xplus + yminus*nx] * hv[xplus + yminus*nx]);
		dv = vv[i] - vv[ix + yminus*nx];
		if (qin > 0.0f)
		{
			if (dv > 0.5f)
			{
				// Conservation of energy head
				vvdudy = vvdudy + 0.5f*(vv[ix + yminus*nx] + vv[xplus + yminus*nx])*(uui[tx][ty]-uub[tx][ty])/dx;//
			}
			else
			{
				// Conservation of momentum
				vvdudy = vvdudy + qin / hum[i] * (uui[tx][ty] - uub[tx][ty]) / dx;
			}
		}
		qin = -0.5f*(vv[i] * hv[i] + vv[xplus + iy*nx] * hv[xplus + iy*nx]);
		dv = vv[ix+yplus*nx] - vv[i];
		if (qin > 0.0f)
		{
			if (dv > 0.5f)
			{
				// Conservation of energy head
				vvdudy = vvdudy - 0.5f*(vv[i] + vv[xplus + iy*nx])*(uui[tx][ty] - uut[tx][ty]) / dx;;
			}
			else
			{
				//Conservation of momentum
				vvdudy = vvdudy + qin / hum[i] * (uui[tx][ty] - uut[tx][ty]) / dx;

			}
		}
		vdudy[i] = vvdudy;
	}

}

__global__ void vdudy_adv2(int nx, int ny, DECNUM dx, DECNUM * hv, DECNUM * hum, DECNUM * uu, DECNUM *vv, DECNUM * vdudy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	DECNUM qin, vvdudy;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uut[16][16];
	__shared__ DECNUM uut2[16][16];
	__shared__ DECNUM uub[16][16];
	__shared__ DECNUM uub2[16][16];

	__shared__ DECNUM vvi[16][16];
	__shared__ DECNUM vvr[16][16];
	__shared__ DECNUM vvb[16][16];
	__shared__ DECNUM vvbr[16][16];
	__shared__ DECNUM hvi[16][16];
	__shared__ DECNUM hvr[16][16];
	__shared__ DECNUM hvb[16][16];
	__shared__ DECNUM hvbr[16][16];
	__shared__ DECNUM humi[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);


		uui[tx][ty] = uu[i];
		uut[tx][ty] = uu[ix + yplus*nx];
		uub[tx][ty] = uu[ix + yminus*nx];
		uut2[tx][ty] = uu[ix + yplus2*nx];
		uub2[tx][ty] = uu[ix + yminus2*nx];

		vvi[tx][ty] = vv[i];
		vvr[tx][ty] = vv[xplus + iy*nx];
		vvb[tx][ty] = vv[ix + yminus*nx];
		vvbr[tx][ty] = vv[xplus + yminus*nx];
		hvi[tx][ty] = hv[i];
		hvr[tx][ty] = hv[xplus + iy*nx];
		hvb[tx][ty] = hv[ix + yminus*nx];
		hvbr[tx][ty] = hv[xplus + yminus*nx];
		humi[tx][ty] = hum[i];


		vvdudy = 0.0f;

		qin = 0.5f*(vvb[tx][ty] * hvb[tx][ty] + vvbr[tx][ty] * hvbr[tx][ty]);
		if (qin > 0.0f)
		{
			vvdudy = vvdudy + qin / humi[tx][ty] * (3 * uui[tx][ty] - 4 * uub[tx][ty] + uub2[tx][ty]) / (2 * dx);
		}
		qin = -0.5f*(vvi[tx][ty] * hvi[tx][ty] + vvr[tx][ty] * hvr[tx][ty]);
		if (qin > 0.0f)
		{
			vvdudy = vvdudy + qin / humi[tx][ty] * (3 * uui[tx][ty] - 4 * uut[tx][ty] + uut2[tx][ty]) / (2 * dx);
		}
		vdudy[i] = vvdudy;
	}

}



__global__ void vdvdy_adv(int nx, int ny, DECNUM dx, DECNUM * hv, DECNUM * hvm, DECNUM * vv, DECNUM * vdvdy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	DECNUM qin,dv, vvdvdy;
	vvdvdy = 0.0f;
	if (ix < nx && ix>0 && iy>0 && iy < (ny-2) )
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		

		qin = 0.5f*(vv[i] * hv[i] + vv[ix + yminus*nx] * hv[ix + yminus*nx]);
		dv = vv[i] - vv[ix + yminus*nx];
		if (qin > 0.0f)
		{
			if (dv > 0.5f)
			{
				vvdvdy = vvdvdy + 0.5f*(vv[ix + (yminus)*nx] - vv[i]) * (vv[i] - vv[ix + (yminus)*nx]) / (dx);
			}
			else
			{
				vvdvdy = vvdvdy + qin / hvm[i] * (vv[i] - vv[ix + (yminus)*nx]) / (dx);
			}
			
		}
		qin = -0.5f*(hv[i] * vv[i] + hv[ix + (yplus)*nx] * vv[ix + (yplus)*nx]);
		dv = vv[ix + (yplus)*nx] - vv[i];
		if (qin > 0.0f)
		{
			if (dv > 0.5f)
			{
				vvdvdy = vvdvdy - 0.5f*(vv[ix + (yplus)*nx] + vv[i]) * (vv[i] - vv[ix + (yplus)*nx]) / (dx);
			}
			else
			{
				vvdvdy = vvdvdy + qin / hvm[i] * (vv[i] - vv[ix + (yplus)*nx]) / (dx);
			}
			
		}
		
	}
	vdvdy[i] = vvdvdy;
}

__global__ void vdvdy_fixbnd(int nx, int ny, DECNUM dx, DECNUM * hv, DECNUM * hvm, DECNUM * vv, DECNUM * vdvdy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	DECNUM qin, dv, vvdvdy;
	
	if (ix < nx && ix>0 && iy > 0 && iy == (ny - 1))
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		qin = 0.5f*(vv[i] * hv[i] + vv[ix + yminus*nx] * hv[ix + yminus*nx]);
		if (qin > 0.0f)
		{
			vvdvdy = vdvdy[i] + qin / hvm[i] * (vv[i] - vv[ix + (yminus)*nx]) / (dx);
		}

	}
	vdvdy[i] = vvdvdy;
}


__global__ void vdvdy_adv2(int nx, int ny, DECNUM dx, DECNUM * hv, DECNUM * hvm, DECNUM * vv, DECNUM * vdvdy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;



	__shared__ DECNUM vvi[16][16];
	__shared__ DECNUM vvb[16][16];
	__shared__ DECNUM vvb2[16][16];
	__shared__ DECNUM vvt[16][16];
	__shared__ DECNUM vvt2[16][16];
	__shared__ DECNUM hvi[16][16];
	__shared__ DECNUM hvb[16][16];
	__shared__ DECNUM hvt[16][16];
	__shared__ DECNUM hvmi[16][16];
	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);
		vvi[tx][ty] = vv[i];
		vvb[tx][ty] = vv[ix + yminus*nx];
		vvb2[tx][ty] = vv[ix + yminus2*nx];
		vvt[tx][ty] = vv[ix + yplus*nx];
		vvt2[tx][ty] = vv[ix + yplus2*nx];
		hvi[tx][ty] = hv[i];
		hvb[tx][ty] = hv[ix + yminus*nx];
		hvt[tx][ty] = hv[ix + yplus*nx];
		hvmi[tx][ty] = hvm[i];


		DECNUM qin, vvdvdy;

		vvdvdy = 0.0f;

		qin = 0.5*(vvi[tx][ty] * hvi[tx][ty] + vvb[tx][ty] * hvb[tx][ty]);
		if (qin > 0.0f)
		{
			vvdvdy = vvdvdy + qin / hvmi[tx][ty] * (3.0f*vvi[tx][ty] - 4.0f*vvb[tx][ty] + vvb2[tx][ty]) / (2 * dx);
		}
		qin = -0.5f*(hvi[tx][ty] * vvi[tx][ty] + hvt[tx][ty] * vvt[tx][ty]);
		if (qin > 0.0f)
		{
			vvdvdy = vvdvdy + qin / hvmi[tx][ty] * (3.0f*vvi[tx][ty] - 4.0f*vvt[tx][ty] + vvt2[tx][ty]) / (2 * dx);
		}
		vdvdy[i] = vvdvdy;
	}
}

__global__ void udvdx_adv(int nx, int ny, DECNUM dx, DECNUM * hu, DECNUM * hvm, DECNUM * uu, DECNUM * vv, DECNUM * udvdx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	DECNUM qin,du, uudvdx;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		uudvdx = 0.0f;
		qin = 0.5*(uu[xminus + iy*nx] * hu[xminus + iy*nx] + uu[xminus + yplus*nx] * hu[xminus + yplus*nx]);
		du = uu[i] - uu[xminus + iy*nx];
		if (qin > 0.0f)
		{
			if (du > 0.5f)
			{
				uudvdx = uudvdx + 0.5f* (uu[xminus+iy*nx]+uu[xminus+yplus*nx]) * (vv[i] - vv[xminus + iy*nx]) / (dx);
			}
			else
			{
				uudvdx = uudvdx + qin / hvm[i] * (vv[i] - vv[xminus + iy*nx]) / (dx);
			}
			
		}
		qin = -0.5*(uu[i] * hu[i] + uu[ix + yplus*nx] * hu[ix + yplus*nx]);
		du = uu[xplus + iy*nx]-uu[i];
		if (qin > 0.0f)
		{
			if (du > 0.5f)
			{
				uudvdx = uudvdx -0.5f*(uu[i]+uu[xplus+iy*nx]) * (vv[i] - vv[xplus + iy*nx]) / (dx);
			}
			else
			{
				uudvdx = uudvdx + qin / hvm[i] * (vv[i] - vv[xplus + iy*nx]) / (dx);
			}
			
		}

		udvdx[i] = uudvdx;
	}


}

__global__ void udvdx_adv2(int nx, int ny, DECNUM dx, DECNUM * hu, DECNUM * hvm, DECNUM * uu, DECNUM * vv, DECNUM * udvdx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;



	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uut[16][16];
	__shared__ DECNUM uul[16][16];
	__shared__ DECNUM uutl[16][16];
	__shared__ DECNUM vvi[16][16];
	__shared__ DECNUM vvl[16][16];
	__shared__ DECNUM vvl2[16][16];
	__shared__ DECNUM vvr[16][16];
	__shared__ DECNUM vvr2[16][16];
	__shared__ DECNUM hui[16][16];
	__shared__ DECNUM hut[16][16];
	__shared__ DECNUM hul[16][16];
	__shared__ DECNUM hutl[16][16];
	__shared__ DECNUM hvmi[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		uui[tx][ty] = uu[i];
		uut[tx][ty] = uu[ix + yplus*nx];
		uul[tx][ty] = uu[xminus + iy*nx];
		uutl[tx][ty] = uu[xminus + yplus*nx];
		vvi[tx][ty] = vv[i];
		vvl[tx][ty] = vv[xminus + iy*nx];
		vvl2[tx][ty] = vv[xminus2 + iy*nx];
		vvr[tx][ty] = vv[xplus + iy*nx];
		vvr2[tx][ty] = vv[xplus2 + iy*nx];
		hui[tx][ty] = hu[i];
		hut[tx][ty] = hu[ix + yplus*nx];
		hul[tx][ty] = hu[xminus + iy*nx];
		hutl[tx][ty] = hu[xminus + yplus*nx];
		hvmi[tx][ty] = hvm[i];


		DECNUM qin, uudvdx;

		uudvdx = 0.0f;
		qin = 0.5*(uul[tx][ty] * hul[tx][ty] + uutl[tx][ty] * hutl[tx][ty]);
		if (qin > 0.0f)
		{
			uudvdx = uudvdx + qin / hvmi[tx][ty] * (3 * vvi[tx][ty] - 4 * vvl[tx][ty] + vvl2[tx][ty]) / (2 * dx);
		}
		qin = -0.5*(uui[tx][ty] * hui[tx][ty] + uut[tx][ty] * hut[tx][ty]);
		if (qin > 0.0f)
		{
			uudvdx = uudvdx + qin / hvmi[tx][ty] * (3 * vvi[tx][ty] - 4 * vvr[tx][ty] + vvr2[tx][ty]) / (2 * dx);
		}

		udvdx[i] = uudvdx;
	}

}


__global__ void smago(int nx, int ny, DECNUM dx, DECNUM rho, DECNUM nuhfac, DECNUM *hh, DECNUM * uu, DECNUM * vv, DECNUM nuh, DECNUM * nuhgrid, DECNUM * DR, int usesmago)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	DECNUM dudx, dudy, dvdx, dvdy, tau;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		if (usesmago == 1)
		{
			dudx = (uu[i] - uu[xminus + iy*nx]) / dx;
			dudy = 0.50f*(uu[ix + yplus*nx] - uu[ix + yminus*nx] + uu[xminus + yplus*nx] - uu[xminus + yminus*nx]) / dx;
			dvdy = (vv[i] - vv[ix + yminus*nx]) / dx;
			dvdx = 0.50f*(vv[xplus + iy*nx] - vv[xminus + iy*nx] + vv[xplus + yminus*nx] - vv[xminus + yminus*nx]) / dx;
			tau = sqrt(2.0f*dudx*dudx + 2.0f*dvdy*dvdy + powf(dudy + dvdx, 2.0f));
			nuhgrid[i] = max(nuh*nuh*dx*dx*tau, nuhfac*hh[i] * powf(DR[i] / rho, 1.0f / 3.0f));
		}
		else
		{
			//nuhgrid[i] = nuh;
			nuhgrid[i] = max(nuh, nuhfac*hh[i] * powf(DR[i] / rho, 1.0f / 3.0f));
		}

	}


}

__global__ void viscou1(int nx, int ny, DECNUM dx, DECNUM rho, DECNUM eps, DECNUM nuhfac, DECNUM * nuhgrid, DECNUM *hh, DECNUM *hum, DECNUM *hvm, DECNUM * DR, DECNUM *uu, int * wetu, DECNUM * viscu)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;



	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		DECNUM nuh1, nuh2;
		nuh1= nuhgrid[xplus+iy*nx];
		nuh2 = nuhgrid[i];
		DECNUM dudx1, dudx2, dudy1, dudy2;

		//if(ix>3)
		//{
		//nnuh = max(nuh, nuhfac*hh[i] * powf(DR[i] / rho, 1.0f / 3.0f));
		//}
		//else
		//{
		//	nnuh=nuh*10;
		//}
		dudx1 = nuh1 * hh[xplus + iy*nx] * (uu[xplus + iy*nx] - uu[i]) / dx;
		dudx2 = nuh2 * hh[i] * (uu[i] - uu[xminus + iy*nx]) / dx;
		//dudy1 = 0.5f*(hvm[i] + hvm[xplus + iy*nx])*(uu[ix + yplus*nx] - uu[i]) / dx;
		//dudy2 = 0.5f*(hvm[ix + yminus*nx] + hvm[xplus + yminus*nx])*(uu[i] - uu[ix + yminus*nx]) / dx;
		//viscu[i] = (1.0f / hum[i]) * (2.0f*(dudx1 - dudx2) / (2.0f*dx));
		viscu[i] = (1.0f / hum[i]) * ((dudx1 - dudx2) / (dx));
		//*wetu[xplus+iy*nx]*wetu[xplus+iy*nx]
	}
}

__global__ void viscou2(int nx, int ny, DECNUM dx, DECNUM rho, DECNUM eps, DECNUM nuhfac, DECNUM * nuhgrid, DECNUM *hh, DECNUM *hum, DECNUM *hvm, DECNUM * DR, DECNUM *uu, int * wetu, DECNUM * viscu)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;



	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		DECNUM nnuh = nuhgrid[i];

		DECNUM nuh1, nuh2, dudx1, dudx2, dudy1, dudy2;

		nuh1 = 0.25f*(nuhgrid[i] + nuhgrid[xplus + iy*nx] + nuhgrid[xplus + yplus*nx] + nuhgrid[ix + yplus*nx]);
		nuh2 = 0.25f*(nuhgrid[i] + nuhgrid[xplus + iy*nx] + nuhgrid[xplus + yminus*nx] + nuhgrid[ix + yminus*nx]);

		//if(ix>3)
		//{
		//nnuh = max(nuh, nuhfac*hh[i] * powf(DR[i] / rho, 1.0f / 3.0f));
		//}
		//else
		//{
		//	nnuh=nuh*10;
		//}
		//dudx1 = hh[xplus + iy*nx] * (uu[xplus + iy*nx] - uu[i]) / dx;
		//dudx2 = hh[i] * (uu[i] - uu[xminus + iy*nx]) / dx;
		dudy1 = nuh1*0.5f*(hvm[i] + hvm[xplus + iy*nx])*(uu[ix + yplus*nx] - uu[i]) / dx;
		dudy2 = nuh2*0.5f*(hvm[ix + yminus*nx] + hvm[xplus + yminus*nx])*(uu[i] - uu[ix + yminus*nx]) / dx;
		viscu[i] = viscu[i] +(1.0f/ hum[i]) * ((dudy1 - dudy2) / dx*wetu[ix + yplus*nx] * wetu[ix + yminus*nx]);

		//*wetu[xplus+iy*nx]*wetu[xplus+iy*nx]
	}
}

__global__ void viscov1(int nx, int ny, DECNUM dx, DECNUM rho, DECNUM eps, DECNUM nuhfac, DECNUM * nuhgrid, DECNUM *hh, DECNUM *hum, DECNUM *hvm, DECNUM * DR, DECNUM *vv, int * wetv, DECNUM * viscv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		DECNUM nuh1, nuh2;
		nuh1 = nuhgrid[ix+yplus*nx];
		nuh2 = nuhgrid[i];

		DECNUM nnuh, dvdx1, dvdx2, dvdy1, dvdy2;

		//nnuh = max(nuh, nuhfac*hh[i] * powf(DR[i] / rho, 1.0f / 3.0f));

		//dvdx1 = 0.5f*(hum[i] + hum[ix + yplus*nx])*(vv[xplus + iy*nx] - vv[i]) / dx;
		//dvdx2 = 0.5f*(hum[xminus + iy*nx] + hum[xminus + yplus*nx])*(vv[i] - vv[xminus + iy*nx]) / dx;
		dvdy1 = nuh1*hh[ix + yplus*nx] * (vv[ix + yplus*nx] - vv[i]) / dx;
		dvdy2 = nuh2*hh[i] * (vv[i] - vv[ix + yminus*nx]) / dx;
		//viscv[i] = nnuh / hvm[i] * ((dvdx1 - dvdx2) / (dx)*wetv[xplus + iy*nx] * wetv[xminus + iy*nx] + (dvdy1 - dvdy2) / dx*wetv[ix + yplus*nx] * wetv[ix + yminus*nx]);
		viscv[i] = (1.0f / hvm[i]) * (dvdy1 - dvdy2) / dx*wetv[ix + yplus*nx] * wetv[ix + yminus*nx];
	}
}
__global__ void viscov2(int nx, int ny, DECNUM dx, DECNUM rho, DECNUM eps, DECNUM nuhfac, DECNUM * nuhgrid, DECNUM *hh, DECNUM *hum, DECNUM *hvm, DECNUM * DR, DECNUM *vv, int * wetv, DECNUM * viscv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		DECNUM nuh1, nuh2;
		nuh1 = 0.25f*(nuhgrid[i] + nuhgrid[xplus + iy*nx] + nuhgrid[xplus + yplus*nx] + nuhgrid[ix + yplus*nx]);
		nuh2 = 0.25f*(nuhgrid[i] + nuhgrid[xminus + iy*nx] + nuhgrid[xminus + yplus*nx] + nuhgrid[ix + yplus*nx]);

		DECNUM nnuh, dvdx1, dvdx2, dvdy1, dvdy2;

		//nnuh = max(nuh, nuhfac*hh[i] * powf(DR[i] / rho, 1.0f / 3.0f));

		dvdx1 = nuh1*0.5f*(hum[i] + hum[ix + yplus*nx])*(vv[xplus + iy*nx] - vv[i]) / dx;
		dvdx2 = nuh2*0.5f*(hum[xminus + iy*nx] + hum[xminus + yplus*nx])*(vv[i] - vv[xminus + iy*nx]) / dx;
		//dvdy1 = hh[ix + yplus*nx] * (vv[ix + yplus*nx] - vv[i]) / dx;
		//dvdy2 = hh[i] * (vv[i] - vv[ix + yminus*nx]) / dx;
		viscv[i] = viscv[i] + (1.0f / hvm[i]) * ((dvdx1 - dvdx2) / (dx)*wetv[xplus + iy*nx] * wetv[xminus + iy*nx] );
	}
}

__global__ void viscovbnd(int nx, int ny, DECNUM * viscv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		if (iy == ny - 2)
		{
			viscv[i] = viscv[ix + yminus*nx];
		}
		if (iy == 0)
		{
			viscv[i] = viscv[ix + yplus*nx];
		}
	}
}



__global__ void eulerustep(int nx, int ny, DECNUM dx, DECNUM dt, DECNUM g, DECNUM rho, DECNUM * zo, DECNUM fc, DECNUM windv, DECNUM Cd, DECNUM *uu, DECNUM * urms, DECNUM *ududx, DECNUM *vdudy, DECNUM *viscu, DECNUM *dzsdx, DECNUM *hu, DECNUM *hum, DECNUM *Fx, DECNUM *vu, DECNUM * ueu_g, DECNUM * vmageu, int *wetu)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	DECNUM ueu;
	DECNUM taubx;
	DECNUM hmin = 0.2;


	__shared__ DECNUM  uui[16][16];


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);




		uui[tx][ty] = uu[i];
		ueu = ueu_g[i];

		__syncthreads;

		int signwind = (windv>0.0) - (windv < 0.0);


		//&& ix>0
		if (wetu[i] == 1)
		{
			taubx = zo[i] * rho*ueu*sqrtf(1.3456f*urms[i] * urms[i] + vmageu[i] * vmageu[i]);

			uui[tx][ty] = uui[tx][ty] - dt*(ududx[i] + vdudy[i] - viscu[i] + g*dzsdx[i] + taubx / (rho*hu[i]) - Fx[i] / (rho*max(hum[i], hmin)) - 1.25f*Cd* signwind *windv*windv / (rho*hum[i]) - fc*vu[i]);

			//viscu[i]=taubx;

		}
		else
		{
			uui[tx][ty] = 0.0f;
			//viscu[i] = 0.0f;

		}
		if (ix > 0)
		{
			uu[i] = uui[tx][ty];

		}
	}

}

__global__ void eulervstep(int nx, int ny, DECNUM dx, DECNUM dt, DECNUM g, DECNUM rho, DECNUM * zo, DECNUM fc, DECNUM windv, DECNUM Cd, DECNUM *vv, DECNUM * urms, DECNUM *udvdx, DECNUM *vdvdy, DECNUM *viscv, DECNUM *dzsdy, DECNUM *hv, DECNUM *hvm, DECNUM *Fy, DECNUM *uv, DECNUM * vev_g, DECNUM * vmagev, int *wetv)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	__shared__ DECNUM  vvi[16][16];
	__shared__ DECNUM  urmsi[16][16];
	__shared__ DECNUM  vmagvi[16][16];
	__shared__ DECNUM  hvmi[16][16];

	if (ix < nx && iy < ny-2)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);



		DECNUM tauby, vev;

		DECNUM hmin = 0.2;

		vvi[tx][ty] = vv[i];
		urmsi[tx][ty] = urms[i];
		vmagvi[tx][ty] = vmagev[i];
		hvmi[tx][ty] = hvm[i];
		int signwind = (windv>0.0)-(windv < 0.0);
		// && ix>0
		if (wetv[i] == 1)
		{
			vev = vev_g[i];

			tauby = zo[i] * rho*vev*sqrtf(1.3456f*urmsi[tx][ty] * urmsi[tx][ty] + vmagvi[tx][ty] * vmagvi[tx][ty]);
			vvi[tx][ty] = vvi[tx][ty] - dt*(udvdx[i] + vdvdy[i] - viscv[i] + g*dzsdy[i] + tauby / (rho*hv[i]) - Fy[i] / (rho*max(hvmi[tx][ty], hmin)) + fc*uv[i] - 1.25f*Cd * signwind * windv * windv / (rho*hvmi[tx][ty]));

			//viscv[i]=tauby;

		}
		else
		{
			vvi[tx][ty] = 0.0f;
			//viscv[i] = 0.0f;
		}
		if (ix > 0)// && iy>0 && iy<ny)
		{
			vv[i] = vvi[tx][ty];
		}//vdvdy[i]=tauby;

	}
}

__global__ void flowsecO_advUV(int nx, int ny, DECNUM *uuold, DECNUM *vvold, DECNUM *uu, DECNUM *vv)
{
	float mindepth;
	
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		
		//
		mindepth = 0.0;//
	}
}


__global__ void continuity(int nx, int ny, DECNUM dx, DECNUM dt, DECNUM eps, DECNUM* qx,  DECNUM* qy, DECNUM* zs, DECNUM *hh, DECNUM *zb, DECNUM * dzsdt)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;



	//DECNUM qxi, qyi, qxm, qym, dzdt;
	DECNUM zz,dzdt;

	__shared__ DECNUM qxi[16][16];
	__shared__ DECNUM qxl[16][16];
	__shared__ DECNUM qyi[16][16];
	__shared__ DECNUM qyb[16][16];
	

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		qxi[tx][ty] = qx[i];
		qyi[tx][ty] = qy[i];
		qxl[tx][ty] = qx[xminus + iy*nx];
		qyb[tx][ty] = qy[ix + yminus*nx];
		




		zz = zs[i];

		
		dzdt = (qxl[tx][ty] - qxi[tx][ty] + qyb[tx][ty] - qyi[tx][ty]) / dx;


		

		if (ix > 0)
		{
			dzsdt[i] = dzdt;


			zs[i] = zz + dzdt*dt;

			//hh[i]=max(hh[i]+dzdt*dt,eps);
		}
	}

}

__global__ void hsbnd(int nx, int ny, DECNUM eps, DECNUM * hh, DECNUM *zb, DECNUM *zs)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;




	__shared__ DECNUM Fi[16][16];
	__shared__ DECNUM Ft[16][16];
	__shared__ DECNUM Fb[16][16];
	__shared__ DECNUM Fr[16][16];

	if (ix < nx && iy < ny)
	{
		
		//hh[i]=Fi[tx][ty];



		hh[i] = max(zb[i] + zs[i], eps);
	}


}

__global__ void uvlatbnd(int nx, int ny, DECNUM * vu, DECNUM * uv, DECNUM * ueu, DECNUM * vev, DECNUM * vmageu, DECNUM * vmagev)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;



	__shared__ DECNUM vut[16][16];
	__shared__ DECNUM vub[16][16];
	__shared__ DECNUM uvt[16][16];
	__shared__ DECNUM uvb[16][16];
	__shared__ DECNUM ueub[16][16];
	__shared__ DECNUM ueut[16][16];
	__shared__ DECNUM vevt[16][16];
	__shared__ DECNUM vevb[16][16];
	__shared__ DECNUM vmagevt[16][16];
	__shared__ DECNUM vmagevb[16][16];
	__shared__ DECNUM vmageut[16][16];
	__shared__ DECNUM vmageub[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		uvt[tx][ty] = uv[ix + yplus*nx];
		uvb[tx][ty] = uv[ix + yminus*nx];
		vut[tx][ty] = vu[ix + yplus*nx];
		vub[tx][ty] = vu[ix + yminus*nx];
		ueut[tx][ty] = ueu[ix + yplus*nx];
		ueub[tx][ty] = ueu[ix + yminus*nx];
		vevt[tx][ty] = vev[ix + yplus*nx];
		vevb[tx][ty] = vev[ix + yminus*nx];
		vmageut[tx][ty] = vmageu[ix + yplus*nx];
		vmageub[tx][ty] = vmageu[ix + yminus*nx];
		vmagevt[tx][ty] = vmagev[ix + yplus*nx];
		vmagevb[tx][ty] = vmagev[ix + yminus*nx];

		if (iy == 0)
		{
			uv[i] = uvt[tx][ty];
			vu[i] = vut[tx][ty];
			ueu[i] = ueut[tx][ty];
			vev[i] = vevt[tx][ty];
			vmageu[i] = vmageut[tx][ty];
			vmagev[i] = vmagevt[tx][ty];
		}
		if (iy == ny - 1)
		{
			uv[i] = uvb[tx][ty];
			vu[i] = vub[tx][ty];
			ueu[i] = ueub[tx][ty];
			vev[i] = vevb[tx][ty];
			vmageu[i] = vmageub[tx][ty];
			vmagev[i] = vmagevb[tx][ty];
		}
	}
}


__global__ void latbnd(int nx, int ny, DECNUM * uu)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	
	__shared__ DECNUM uut[16][16];
	__shared__ DECNUM uub[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int yminus = mminus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		uut[tx][ty] = uu[ix + yplus*nx];
		uub[tx][ty] = uu[ix + yminus*nx];

		if (iy == 0)
		{
			uu[i] = uut[tx][ty];
		}
		if (iy == ny - 1)
		{
			uu[i] = uub[tx][ty];
		}
	}
	//
}


__global__ void vvlatbnd(int nx, int ny, DECNUM * uu)
{
	// Neumann_v boundary
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	__shared__ DECNUM uut[16][16];
	__shared__ DECNUM uub[16][16];

	if (ix < nx && iy < ny)
	{
		unsigned int yminus = mminus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		uut[tx][ty] = uu[ix + yplus*nx];
		uub[tx][ty] = uu[ix + yminus*nx];

		if (iy == 0)
		{
			uu[i] = uut[tx][ty];
		}
		if (iy == ny - 2)
		{
			uu[i] = uub[tx][ty];
		}
	}
	//
}



__global__ void uuvvzslatbnd(int nx, int ny, DECNUM * uu, DECNUM * vv, DECNUM *zs)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	int tx = threadIdx.x;
	int ty = threadIdx.y;



	__shared__ DECNUM vvr[16][16];
	__shared__ DECNUM vvb[16][16];
	__shared__ DECNUM vvb2[16][16];
	__shared__ DECNUM vvt[16][16];
	__shared__ DECNUM uut[16][16];
	__shared__ DECNUM uub[16][16];
	__shared__ DECNUM zst[16][16];
	__shared__ DECNUM zsb[16][16];
	__shared__ DECNUM zsl[16][16];

	if (ix < nx && iy < ny)
	{

		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus = pplus(iy, ny);


		uut[tx][ty] = uu[ix + yplus*nx];
		uub[tx][ty] = uu[ix + yminus*nx];
		vvr[tx][ty] = vv[xplus + iy*nx];
		vvt[tx][ty] = vv[ix + yplus*nx];
		vvb[tx][ty] = vv[ix + yminus*nx];
		vvb2[tx][ty] = vv[ix + yminus2*nx];
		zst[tx][ty] = zs[ix + yplus*nx];
		zsb[tx][ty] = zs[ix + yminus*nx];
		zsl[tx][ty] = zs[xminus + iy*nx];

		//F[i]=Fi[tx][ty]*wet;
		if (iy == 0)
		{
			uu[i] = uut[tx][ty];
			vv[i] = vvt[tx][ty];
			zs[i] = zst[tx][ty];
		}
		if (iy == ny - 1)
		{
			uu[i] = uub[tx][ty];
			vv[i] = vvb[tx][ty];
			zs[i] = zsb[tx][ty];
		}
		if (iy==ny-2)
		{
			//vv[i]=vvb[tx][ty];// THis is to follow XBeach definition although I don't really agree with it 
							  // It should be that vv(i,ny-1)=vv(i,ny-2) end of story
		}
		if (ix == 0)
		{
			//vv[i]=vvr[tx][ty];//Imcompatible with abs_2d front boundary 
		}
		if (ix == nx - 1)
		{
			//zs[i]=zsl[tx][ty];//Need to fix 
		}
	}

}
__global__ void CalcQFlow(int nx, int ny, DECNUM * uu, DECNUM* hu, DECNUM* vv, DECNUM* hv, DECNUM *qx, DECNUM * qy)

{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;



	__shared__ DECNUM uui[16][16];
	
	__shared__ DECNUM vvi[16][16];
	
	__shared__ DECNUM hui[16][16];
	
	__shared__ DECNUM hvi[16][16];
	

	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);

		uui[tx][ty] = uu[i];
		vvi[tx][ty] = vv[i];
		
		hui[tx][ty] = hu[i];
		
		hvi[tx][ty] = hv[i];
		

		qx[i] = uui[tx][ty] * hui[tx][ty];
		qy[i] = vvi[tx][ty] * hvi[tx][ty];

	}

}

__global__ void wrkuu2Ocorr(int nx, int ny, float dx, float dt,float eps, DECNUM* uu, DECNUM* uuold, DECNUM* hum,DECNUM *zs, DECNUM *zb, DECNUM *qx, DECNUM * qy,DECNUM *wrk1, DECNUM*wrk2)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	
	float delta1, delta2;
	float mindepth,minzs,minzb;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uoi[16][16];

	

	

	if ( ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);

		wrk1[i] = 0.0f;
		wrk2[i] = 0.0f;

		minzs = min(min(min(min(zs[xminus2 + iy*nx], zs[xminus + iy*nx]), zs[i]), zs[xplus + iy*nx]), zs[xplus2 + iy*nx]);
		minzb = min(min(min(min(zb[xminus2 + iy*nx], zb[xminus + iy*nx]), zb[i]), zb[xplus + iy*nx]), zb[xplus2 + iy*nx]);

		mindepth = minzs + minzb;

		if (ix>0 && mindepth > eps*2.0f)
		{

			if ((qx[i] + qx[xminus + iy*nx]) > 0.0f)
			{
				delta1 = (uu[i] - uuold[xminus + iy*nx]);
				delta2 = (uu[xminus + iy*nx] - uuold[xminus2 + iy*nx]);
				wrk1[i] = 0.5f*minmod(delta1, delta2);
			}
			if ((qx[i] + qx[xminus + iy*nx]) < 0.0f)
			{
				delta1 = (uuold[i] - uu[xminus + iy*nx]);
				delta2 = (uuold[xplus + iy*nx] - uu[i]);
				wrk1[i] = -0.5f*minmod(delta1, delta2);
			}
		}
		minzs = min(min(min(zs[ix + iy*nx], zs[ix + yminus*nx]), zs[ix + yplus*nx]), zs[ix + yplus2*nx]);
		minzs = min(minzs, min(min(min(zs[xplus + iy*nx], zs[xplus + yminus*nx]), zs[xplus + yplus*nx]), zs[xplus + yplus2*nx]));

		minzb = min(min(min(zb[ix + iy*nx], zb[ix + yminus*nx]), zb[ix + yplus*nx]), zb[ix + yplus2*nx]);
		minzb = min(minzb, min(min(min(zb[xplus + iy*nx], zb[xplus + yminus*nx]), zb[xplus + yplus*nx]), zb[xplus + yplus2*nx]));

		mindepth = minzs + minzb;

		
		if (ix>0 && iy>0 && mindepth > eps*2.0f)
		{

			if ((qy[xplus + iy*nx] + qy[i]) > 0.0f)
			{
				delta1 = (uu[ix + yplus*nx] - uuold[i]);
				delta2 = (uu[i] - uuold[ix + yminus*nx]);
				wrk2[i] = 0.5f*minmod(delta1, delta2);
			}
			if ((qy[xplus + iy*nx] + qy[i]) < 0.0f)
			{
				delta1 = (uuold[ix + yplus*nx] - uu[i]);
				delta2 = (uuold[ix + yplus2*nx] - uu[ix + yplus*nx]);
				wrk2[i] = -0.5f*minmod(delta1, delta2);
			}

		}
		
	}

}

__global__ void wrkvv2Ocorr(int nx, int ny, float dx, float dt, float eps, DECNUM* vv, DECNUM* vvold, DECNUM* hvm, DECNUM *zs, DECNUM *zb, DECNUM *qx, DECNUM * qy, DECNUM *wrk1, DECNUM*wrk2)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	float delta1, delta2;
	float mindepth, minzs, minzb;
	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uoi[16][16];



	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);

		wrk1[i] = 0.0f;
		wrk2[i] = 0.0f;

		minzs = min(min(min(min(zs[ix + yminus2*nx], zs[ix + yminus*nx]), zs[ix + iy*nx]), zs[ix + yplus*nx]), zs[ix + yplus2*nx]);
		minzb = min(min(min(min(zb[ix + yminus2*nx], zb[ix + yminus*nx]), zb[ix + iy*nx]), zb[ix + yplus*nx]), zb[ix + yplus2*nx]);

		mindepth = minzs + minzb;

		if (ix>0 && iy>0 && iy <ny-2 && mindepth > eps*2.0f)
		{

			if ((qy[i] + qy[ix + yminus*nx]) > 0.0f)
			{
				delta1 = (vv[i] - vvold[ix + yminus*nx]);
				delta2 = (vv[ix + yminus*nx] - vvold[ix + yminus2*nx]);
				wrk1[i] = 0.5f*minmod(delta1, delta2);
			}
			if ((qy[i] + qy[ix + yminus*nx]) < 0.0f)
			{
				delta1 = (vvold[i] - vv[ix + yminus*nx]);
				delta2 = (vvold[ix + yplus*nx] - vv[i]);
				wrk1[i] = -0.5f*minmod(delta1, delta2);
			}

		}

		minzs = min(min(min(zs[xminus + iy*nx], zs[ix + iy*nx]), zs[xplus + iy*nx]), zs[xplus2 + iy*nx]);
		minzs = min(minzs, min(min(min(zs[xminus + yplus*nx], zs[ix + yplus*nx]), zs[xplus + yplus*nx]), zs[xplus2 + yplus*nx]));

		minzb = min(min(min(zb[xminus + iy*nx], zb[ix + iy*nx]), zb[xplus + iy*nx]), zb[xplus2 + iy*nx]);
		minzb = min(minzb, min(min(min(zb[xminus + yplus*nx], zb[ix + yplus*nx]), zb[xplus + yplus*nx]), zb[xplus2 + yplus*nx]));

		mindepth = minzs + minzb;

		if (ix>0 && iy>0 && iy <ny - 2 && mindepth > eps*2.0f)
		{

			if ((qx[ix + yplus*nx] + qx[i]) > 0.0f)
			{
				delta1 = (vv[xplus + iy*nx] - vvold[i]);
				delta2 = (vv[i] - vvold[xminus + iy*nx]);
				wrk2[i] = 0.5f*minmod(delta1, delta2);
			}
			if ((qx[ix + yplus*nx] + qx[i]) < 0.0f)
			{
				delta1 = (vvold[xplus + iy*nx] - vv[i]);
				delta2 = (vvold[xplus2 + iy*nx] - vv[xplus + iy*nx]);
				wrk2[i] = -0.5f*minmod(delta1, delta2);
			}

		}

	}

}

__global__ void wrkzs2Ocorr(int nx, int ny, float dx, float dt,float eps, DECNUM* zs, DECNUM* zsold,  DECNUM *uu, DECNUM * vv, DECNUM *zb, DECNUM *wrk1, DECNUM*wrk2)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	float delta1, delta2;
	float mindepth, minzs, minzb;
	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uoi[16][16];



	if (ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);

		wrk1[i] = 0.0f;
		wrk2[i] = 0.0f;

		minzs = min(min(min( zs[xminus + iy*nx], zs[i]), zs[xplus + iy*nx]), zs[xplus2 + iy*nx]);
		minzb = min(min(min( zb[xminus + iy*nx], zb[i]), zb[xplus + iy*nx]), zb[xplus2 + iy*nx]);

		mindepth = minzs + minzb;

		
		if (ix>0 && mindepth > eps*2.0f)
		{
			if ((uu[i]) > 0.0f)
			{
				delta1 = (zs[xplus + iy*nx] - zsold[i]);
				delta2 = (zs[i] - zsold[xminus + iy*nx]);
				wrk1[i] = uu[i] * 0.5f*minmod(delta1, delta2);
			}
			if ((uu[i]) < 0.0f)
			{
				delta1 = (zsold[xplus2 + iy*nx] - zs[xplus + iy*nx]);
				delta2 = (zsold[xplus + iy*nx] - zs[i]);
				wrk1[i] = -1.0f*uu[i] * 0.5f*minmod(delta1, delta2);
			}
		}

		minzs = min(min(min(zs[ix + yminus*nx], zs[i]), zs[ix + yplus*nx]), zs[ix + yplus2*nx]);
		minzb = min(min(min(zb[ix + yminus*nx], zb[i]), zb[ix + yplus*nx]), zb[ix + yplus2*nx]);

		mindepth = minzs + minzb;


		if (ix>0 && iy>0 && mindepth > eps*2.0f)
		{
			if ((vv[i]) > 0.0f)
			{
				delta1 = (zs[ix + yplus*nx] - zsold[i]);
				delta2 = (zs[i] - zsold[ix + yminus*nx]);
				wrk2[i] = vv[i] * 0.5f*minmod(delta1, delta2);
			}
			if ((vv[i]) < 0.0f)
			{
				delta1 = (zsold[ix + yplus2*nx] - zs[ix + yplus*nx]);
				delta2 = (zsold[ix + yplus*nx] - zs[i]);
				wrk2[i] = vv[i] * 0.5f*minmod(delta1, delta2);
			}
		}


	}

}

__global__ void uu2Ocorr(int nx, int ny, float dx, float dt, DECNUM* uu, DECNUM* uuold, DECNUM* hum, DECNUM *qx, DECNUM * qy, DECNUM *wrk1, DECNUM*wrk2)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	float qe, qw, qs, qn;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uoi[16][16];



	if (ix>0 && iy>0 && ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);
		qe = 0.5f*(qx[xplus + iy*nx] + qx[i]);
		qw = 0.5f*(qx[i] + qx[xminus + iy*nx]);
		qn = 0.5f*(qy[ix + yminus*nx] + qy[xplus + yminus*nx]);
		qs = 0.5f*(qy[i]+qy[xminus+iy*nx]);

		uu[i] = uu[i] - dt / hum[i] * ((qe*wrk1[xplus + iy*nx] - qw*wrk1[i]) / dx + (qs*wrk2[i] - qn*wrk2[ix + yminus*nx]) / dx);

	}
}

__global__ void vv2Ocorr(int nx, int ny, float dx, float dt, DECNUM* vv, DECNUM* vvold, DECNUM* hvm, DECNUM *qx, DECNUM * qy, DECNUM *wrk1, DECNUM*wrk2)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	float qe, qw, qs, qn;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uoi[16][16];



	if (ix>0 && iy>0 && ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);
		qe = 0.5f*(qx[ix + yplus*nx] + qx[i]);
		qw = 0.5f*(qx[xminus+iy*nx] + qx[xminus + yplus*nx]);
		qn = 0.5f*(qy[i] + qy[ix + yminus*nx]);
		qs = 0.5f*(qy[ix+yplus*nx] + qy[i]);

		vv[i] = vv[i] - dt / hvm[i] * ((qe*wrk2[i] - qw*wrk2[xminus+iy*nx]) / dx + (qs*wrk1[ix+yplus*nx] - qn*wrk1[i]) / dx);

	}
}

__global__ void zs2Ocorr(int nx, int ny, float dx, float dt, DECNUM* zs, DECNUM *qx, DECNUM *qy, DECNUM *wrk1, DECNUM*wrk2)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;


	float qe, qw, qs, qn;

	__shared__ DECNUM uui[16][16];
	__shared__ DECNUM uoi[16][16];



	if (ix>0 && iy>0 && ix < nx && iy < ny)
	{
		unsigned int xminus = mminus(ix, nx);
		unsigned int xminus2 = mminus2(ix, nx);
		unsigned int xplus = pplus(ix, nx);
		unsigned int xplus2 = pplus2(ix, nx);
		unsigned int yminus = mminus(iy, ny);
		unsigned int yplus = pplus(iy, ny);
		unsigned int yminus2 = mminus2(iy, ny);
		unsigned int yplus2 = pplus2(iy, ny);
		

		zs[i] = zs[i] - dt * ((wrk1[i] - wrk1[xminus+iy*nx]) / dx + (wrk2[i] - wrk2[ix + yminus*nx]) / dx);

		//Update fluxes although Im not sure if there is a point to do it since Qx and Qy are not used elsewhere 
		qx[i] = qx[i] + wrk1[i];
		qy[i] = qy[i] + wrk2[i];
	}
}