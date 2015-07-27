#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>
#include "XBeachGPU.h"



template <class T> const T& max (const T& a, const T& b) {
  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class T> const T& min(const T& a, const T& b) {
	return !(b<a) ? a : b;     // or: return comp(a,b)?b:a; for version (2)
}

extern "C" int signC(DECNUM x)
{
	return((x>0.0f) - (x<0.0f));
}


extern "C" void offshorebndWavCPU(int nx, int ny,int ntheta,DECNUM totaltime,DECNUM Trep,DECNUM *St,DECNUM *&sigm, DECNUM *&ee)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			DECNUM taper=min(totaltime/100.0,1.0);
	
			sigm[i]=2*pi/Trep;
			for (int itheta=0; itheta<ntheta; itheta++)
			{
				ee[0+iy*nx+itheta*ny*nx]=St[iy+itheta*ny]*taper;
			}
		}
	}

}

extern "C" void set_bndCPU(int nx, int ny,DECNUM Trep,int ntheta,DECNUM * theta,DECNUM *&sigm)
{
	for (int ix=0; ix<nx; ix++)
		{
			for (int iy=0; iy<ny; iy++)
				{
					int i=ix+iy*nx;
					sigm[i]=2.0*pi/Trep;
			}
	}
	
			
}


extern "C" void sanityCPU(int nx, int ny,DECNUM eps,DECNUM * hh,DECNUM * sigm, int ntheta,DECNUM * &ee)
{
	//printf("%f\t",ee[0+16*nx+6*nx*ny]);
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
	

		    for (int itheta=0; itheta<ntheta; itheta++)
			{
				//ssigm[tx][ty]=ssigm[tx][ty]+sigt[i+itheta*nx*ny]/ntheta;
				ee[i+itheta*nx*ny]=max(ee[i+itheta*nx*ny],0.000f);
			}
			//sigm[i]=max(ssigm[tx][ty],0.0001f);
		}
	}


        


}


extern "C" void dispersion_initCPU(int nx,int ny,DECNUM twopi,DECNUM g,DECNUM aphi,DECNUM bphi,DECNUM * sigm,DECNUM * hh,DECNUM * &cg)
{
	DECNUM L0, L1,L2, errdisp;
	for (int ix=0; ix<nx; ix++)
		{
			for (int iy=0; iy<ny; iy++)
				{
					
					unsigned int i=ix+iy*nx;
					//int tx =threadIdx.x;
					//int ty= threadIdx.y;
	

					//__shared__ DECNUM sigmi[16][16];
					//__shared__ DECNUM hhi[16][16];
	//sigmi[tx][ty]=sigm[i];
	//hhi[tx][ty]=hh[i];
  
        
        L0=twopi*g/(sigm[i]*sigm[i]);
        L1=L0;
        
        errdisp=1000.0f;
        //while (errdisp > 0.0001f)
        for (int k=1; k<200; k++)
        {
           L2        = L0*tanh(2*pi*hh[i]/L1);        
            errdisp       = abs(L2 - L1);
            L1 = (L1*aphi + L2*bphi);//          ! Golden ratio
            if (errdisp <= 0.0001f)
            {
				break;
			}
			if(k==199)
			{
				L1=L0*powf(tanh(powf(sigm[i]*sigm[i]*hh[i]/g,3.0f/4.0f)),2.0f/3.0f);
				break;
			}
        }
        

		//L1=L0*powf(tanh(powf(sigmi[tx][ty]*sigmi[tx][ty]*hhi[tx][ty]/g,3/4)),2/3);
        
		DECNUM kk=2*pi/L1;
		//k[i]  = kk;
		DECNUM cc=sigm[i]/kk;
        //c[i]  = cc;
		DECNUM kkhh=min(kk*hh[i],10.0f);
       // kh[i]   = kkhh;
		DECNUM s2kh=sinhf(2.0f*kkhh);
        //sinh2kh[i]=s2kh;
		cg[i] = cc*(0.5f+kkhh/s2kh);
			}
	}
 
}

extern "C" void dispersionCPU(int nx,int ny,DECNUM twopi,DECNUM g,DECNUM aphi,DECNUM bphi,DECNUM * sigm,DECNUM * hh,DECNUM * &k,DECNUM * &c,DECNUM * &kh,DECNUM * &sinh2kh,DECNUM * &cg)
{
	DECNUM L0, L1,L2, errdisp;
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			//int tx =threadIdx.x;
			//int ty= threadIdx.y;
	
			
        
			L0=twopi*g/(sigm[i]*sigm[i]);
			L1=L0;
        
			errdisp=1000.0f;
			//while (errdisp > 0.0001f)
			for (int n=1; n<200; n++)
			{	
				L2= L0*tanh(2*pi*hh[i]/L1);        
				errdisp= abs(L2 - L1);
				L1 = (L1*aphi + L2*bphi);//          ! Golden ratio
				if (errdisp <= 0.0001f)
				{
					break;
				}
				if(n==199)
				{
					L1=L0*powf(tanh(powf(sigm[i]*sigm[i]*hh[i]/g,3.0f/4.0f)),2.0f/3.0f);
					break;
				}
			}
        

			//L1=L0*powf(tanh(powf(sigmi[tx][ty]*sigmi[tx][ty]*hhi[tx][ty]/g,3/4)),2/3);
        
			DECNUM kk=2*pi/L1;
			k[i]  = kk;
			DECNUM cc=sigm[i]/kk;
			c[i]  = cc;
			DECNUM kkhh=min(kk*hh[i],10.0f);
			kh[i]   = kkhh;
			DECNUM s2kh=sinhf(2.0f*kkhh);
			sinh2kh[i]=s2kh;
			cg[i] = cc*(0.5f+kkhh/s2kh);
		}
	}
 
}


extern "C" void calcwciCPU(int nx, int ny, DECNUM wci, DECNUM hwci, DECNUM * hh, DECNUM * &wcig)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
	
	
			wcig[i]=wci*min(hh[i]/hwci,1.0f);
		}
	}
	
}

extern "C" void slopesCPU(int nx,int ny,DECNUM dx,DECNUM * hh,DECNUM * uu,DECNUM * vv,DECNUM * &dhdx,DECNUM * &dhdy,DECNUM * &dudx,DECNUM * &dudy,DECNUM * &dvdx,DECNUM * &dvdy)//
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
	

			
            dhdx[i]=slopes2DxCPU(nx,dx,i,ix,iy,hh);
            dudx[i]=slopes2DxCPU(nx,dx,i,ix,iy,uu);
            dvdx[i]=slopes2DxCPU(nx,dx,i,ix,iy,vv);

			dhdy[i]=slopes2DyCPU(nx,ny,dx,i,ix,iy,hh);
       	    dudy[i]=slopes2DyCPU(nx,ny,dx,i,ix,iy,uu);
			dvdy[i]=slopes2DyCPU(nx,ny,dx,i,ix,iy,vv);
		}
	}
      
        
	
}
extern "C" DECNUM slopes2DxCPU(int nx,DECNUM dx,int i,int ix, int iy,DECNUM * hh)
{


	DECNUM dhdx;
	unsigned int xminus=mminusC(ix,nx);
	unsigned int xplus=pplusC(ix,nx);
	



	dhdx=(hh[xplus+iy*nx]-hh[xminus+iy*nx])/((xplus-xminus)*dx);
	return(dhdx);
}
extern "C" DECNUM slopes2DyCPU(int nx,int ny,DECNUM dx,int i,int ix, int iy,DECNUM * hh)
{
	DECNUM dhdy;
	
	unsigned int yminus=mminusC(iy,ny);
	unsigned int yplus=pplusC(iy,ny);



	dhdy=(hh[ix+yplus*nx]-hh[ix+yminus*nx])/((yplus-yminus)*dx);
	return(dhdy);
}

extern "C" void propagthetaCPU(int nx,int ny,int ntheta,DECNUM * wci,DECNUM * &ctheta,DECNUM *cxsth,DECNUM *sxnth,DECNUM *dhdx,DECNUM *dhdy,DECNUM *dudx,DECNUM *dudy,DECNUM *dvdx,DECNUM *dvdy,DECNUM *sigm,DECNUM *kh)//
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			//int tx =threadIdx.x;
			//int ty= threadIdx.y;
			DECNUM ctht;

			
			DECNUM sigmiosinh2kh=sigm[i]/sinhf(2.0f*kh[i]);
	 
			
			

			for (int itheta=0; itheta<ntheta; itheta++)
			{

				//cx[i+itheta*nx*ny] =ci*cxsth[itheta]+uui;
				//cy[i+itheta*nx*ny] =ci*sxnth[itheta]+vvi;
				ctht= (sigmiosinh2kh)*(dhdx[i]*sxnth[itheta]-dhdy[i]*cxsth[itheta]) + wci[i]*(cxsth[itheta]*(sxnth[itheta]*dudx[i] - cxsth[itheta]*dudy[i]) + sxnth[itheta]*(sxnth[itheta]*dvdx[i] - cxsth[itheta]*dvdy[i]));
				ctheta[i+itheta*nx*ny]=min(max(ctht,-0.25f*sigm[i]),0.25f*sigm[i]);
			}
		}
	}

}

extern "C" void actionCPU(int ntheta,int nx,int ny,DECNUM * &ee,DECNUM * sigm)
{
	
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
			
			
	
	
			for (int itheta=0; itheta<ntheta; itheta++)
			{	
				ee[i+itheta*nx*ny] = ee[i+itheta*nx*ny]/sigm[i]; 
			}
		}
	}

}

extern "C" void xadvecupwindCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM * wci,DECNUM *ee,DECNUM *cg,DECNUM *cxsth,DECNUM *uu,DECNUM * &xadvec)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
	
	
			DECNUM dxplus_i = 1.0/dx;
			DECNUM dxcent_i = 1.0/(2*dx);
			DECNUM xxadvec;
			DECNUM costhet;
			DECNUM arrinx, arrminx, arrmaxx;
			DECNUM cgx,cgxmin;
	
	
			unsigned int xminus=mminusC(ix,nx);//max(ix-1,0);
			unsigned int xplus=pplusC(ix,nx);//min(ix+1,nx-1);
			unsigned int yminus=mminusC(iy,ny);//max(iy-1,0);
			unsigned int yplus=pplusC(iy,ny);//min(iy+1,ny-1);

			DECNUM ccg, ccgxmin, ccgxmax,uui, uuixmax,uuixmin;

			ccg=cg[i];
			ccgxmin=cg[xminus+iy*nx];
			ccgxmax=cg[xplus+iy*nx];
	
	
			uui=uu[i]*wci[i];
			uuixmin = uu[xminus + iy*nx] * wci[xminus + iy*nx];
			uuixmax = uu[xplus + iy*nx] * wci[xplus + iy*nx];
	

			for (int itheta=0; itheta<ntheta; itheta++)
			{
				costhet=cxsth[itheta];
		
				cgx=0.5*(ccg*costhet+uui+ccgxmax*costhet+uuixmax);
				cgxmin=0.5*(ccg*costhet+uui+ccgxmin*costhet+uuixmin);
				xxadvec=0;
		
		
				arrinx=	ee[i+itheta*nx*ny]*max(cgx,0.0f)+ee[xplus+iy*nx+itheta*nx*ny]*min(cgx,0.0f);
				arrminx=ee[xminus+iy*nx+itheta*nx*ny]*max(cgxmin,0.0f)+ee[i+itheta*nx*ny]*min(cgxmin,0.0f);
		
		

				xxadvec=(arrinx-arrminx)*dxplus_i;
				xadvec[i+itheta*nx*ny]=xxadvec;
			}
		}
	}
	
}

extern "C" void xadvecupwind2CPU(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *cxsth, DECNUM *uu, DECNUM * xadvec)
{
	for (int ix = 0; ix < nx; ix++)
	{
		for (int iy = 0; iy < ny; iy++)
		{
			unsigned int i = ix + iy*nx;


			DECNUM dxplus_i = 1.0f / dx;
			DECNUM dxcent_i = 1.0f / (2 * dx);
			DECNUM xxadvec;
			DECNUM costhet;
			DECNUM arrinx, arrminx, arrmaxx;
			DECNUM cgx, cgxmin;
			DECNUM ccg;
			DECNUM ccgxmin;
			DECNUM ccgxmax;
			DECNUM uui;

			DECNUM uuixmin;
			DECNUM uuixmax;

			unsigned int xminus = mminusC(ix, nx);//max(ix-1,0);
			unsigned int xplus = pplusC(ix, nx);//min(ix+1,nx-1);
			unsigned int xplus2 = pplus2C(ix, nx);//min(ix+1,nx-1);
			unsigned int xminus2 = mminus2C(ix, nx);//max(ix-1,0);
			unsigned int yminus = mminusC(iy, ny);//max(iy-1,0);
			unsigned int yplus = pplusC(iy, ny);//min(iy+1,ny-1);


			ccg = cg[i];
			ccgxmin = cg[xminus + iy*nx];
			ccgxmax = cg[xplus + iy*nx];


			uui = uu[i] * wci[i];
			uuixmin = uu[xminus + iy*nx] * wci[i];
			uuixmax = uu[xplus + iy*nx] * wci[i];

			






			for (int itheta = 0; itheta < ntheta; itheta++)
			{
				costhet = cxsth[itheta];

				cgx = 0.5f*(ccg * costhet + uui + ccgxmax * costhet + uuixmax);
				cgxmin = 0.5f*(ccg * costhet + uui + ccgxmin * costhet + uuixmin);
				xxadvec = 0;



				if (cgx > 0.0f)
				{
					arrinx = (1.5*ee[i + itheta*nx*ny] - 0.5*ee[xminus + iy*nx + itheta*nx*ny]);
					if (arrinx < 0.0f)
					{
						arrinx = ee[i + itheta*nx*ny];
					}
					arrinx = arrinx*cgx;
				}
				else
				{
					arrinx = 1.5*ee[xplus + iy*nx + itheta*nx*ny] - 0.5*ee[xplus2 + iy*nx + itheta*nx*ny];
					if (arrinx < 0.0f)
					{
						arrinx = ee[xplus + iy*nx + itheta*nx*ny];
					}
					arrinx = arrinx*cgx;
				}
				if (cgxmin > 0.0f)
				{
					arrminx = 1.5*ee[xminus + iy*nx + itheta*nx*ny] - 0.5*ee[xminus2 + iy*nx + itheta*nx*ny];
					if (arrminx < 0.0f)
					{
						arrminx = ee[xminus + iy*nx + itheta*nx*ny];
					}
					arrminx = arrminx*cgxmin;
				}
				else
				{
					arrminx = 1.5*ee[i + itheta*nx*ny] - 0.5*ee[xplus + iy*nx + itheta*nx*ny];
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


}


extern "C" void yadvecupwind2CPU(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *sxnth, DECNUM *vv, DECNUM * yadvec)
{
	for (int ix = 0; ix < nx; ix++)
	{
		for (int iy = 0; iy < ny; iy++)
		{

			unsigned int i = ix + iy*nx;


			DECNUM dxplus_i = 1.0f / dx;
			DECNUM dxcent_i = 1.0f / (2.0f*dx);
			DECNUM yyadvec;
			DECNUM sinthet;
			DECNUM  arriny, arrminy, arrmaxy;
			DECNUM cgy, cgymin;
			DECNUM ccg;
			DECNUM ccgymin;
			DECNUM ccgymax;

			DECNUM vvi;
			DECNUM vviymin;
			DECNUM vviymax;

			unsigned int xminus = mminusC(ix, nx);
			unsigned int xplus = pplusC(ix, nx);
			unsigned int yminus = mminusC(iy, ny);
			unsigned int yplus = pplusC(iy, ny);
			unsigned int yminus2 = mminus2C(iy, ny);
			unsigned int yplus2 = pplus2C(iy, ny);



			ccg = cg[i];
			ccgymin = cg[ix + (yminus)*nx];
			ccgymax = cg[ix + (yplus)*nx];

			vvi = wci[i] * vv[i];
			vviymin = wci[i] * vv[ix + (yminus)*nx];
			vviymax = wci[i] * vv[ix + (yplus)*nx];
			

			for (int itheta = 0; itheta < ntheta; itheta++)
			{

				sinthet = sxnth[itheta];
				yyadvec = 0;
				cgy = 0.5f*(ccg * sinthet + vvi + ccgymax * sinthet + vviymax);
				cgymin = 0.5f*(ccg * sinthet + vvi + ccgymin * sinthet + vviymin);


				if (cgy > 0.0f)
				{
					arriny = 1.5*ee[i + itheta*nx*ny] - 0.5*ee[ix + yminus*nx + itheta*nx*ny];
					if (arriny < 0.0f)
					{
						arriny = ee[i + itheta*nx*ny];
					}
					arriny = arriny*cgy;
				}
				else
				{
					arriny = 1.5*ee[ix + yplus*nx + itheta*nx*ny] - 0.5*ee[ix + yplus2*nx + itheta*nx*ny];
					if (arriny < 0.0f)
					{
						arriny = ee[ix + yplus*nx + itheta*nx*ny];
					}
					arriny = arriny*cgy;
				}
				if (cgymin > 0.0f)
				{
					arrminy = 1.5*ee[ix + yminus*nx + itheta*nx*ny] - 0.5*ee[ix + yminus2*nx + itheta*nx*ny];
					if (arrminy < 0.0f)
					{
						arrminy = ee[ix + yminus*nx + itheta*nx*ny];
					}
					arrminy = arrminy*cgymin;
				}
				else
				{
					arrminy = 1.5*ee[i + itheta*nx*ny] - 0.5*ee[ix + yplus*nx + itheta*nx*ny];
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

}

extern "C" void yadvecupwindCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM * wci,DECNUM *ee,DECNUM *cg,DECNUM *sxnth,DECNUM *vv,DECNUM * &yadvec)
{

	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
	
	
			DECNUM dxplus_i = 1.0f/dx;
			DECNUM dxcent_i = 1.0f/(2.0f*dx);
			DECNUM yyadvec;
			DECNUM sinthet;
			DECNUM  arriny, arrminy, arrmaxy;
			DECNUM cgy,cgymin;
			DECNUM ccg;
			DECNUM ccgymin;
			DECNUM ccgymax;
	
			DECNUM vvi;
			DECNUM vviymin;
			DECNUM vviymax;

			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);


		
			ccg=cg[i];
			ccgymin=cg[ix+(yminus)*nx];
			ccgymax=cg[ix+(yplus)*nx];
		
			vvi=wci[i]*vv[i];
			vviymin=wci[i]*vv[ix+(yminus)*nx];
			vviymax=wci[i]*vv[ix+(yplus)*nx];
	

			for (int itheta=0; itheta<ntheta; itheta++)
			{
		
				sinthet=sxnth[itheta];
				yyadvec=0;
				cgy=0.5f*(ccg*sinthet+vvi+ccgymax*sinthet+vviymax);
				cgymin=0.5f*(ccg*sinthet+vvi+ccgymin*sinthet+vviymin);
		
				arriny=ee[i+itheta*nx*ny]*max(cgy,0.0f)+ee[ix+yplus*nx+itheta*nx*ny]*min(cgy,0.0f);
				arrminy=ee[ix+yminus*nx+itheta*nx*ny]*max(cgymin,0.0f)+ee[i+itheta*nx*ny]*min(cgymin,0.0f);
				//if (iy > 0 && iy<ny)
				{
					yyadvec = (arriny - arrminy)*dxplus_i;
					yadvec[i + itheta*nx*ny] = yyadvec;
				}
			}
		}
	}
	
	

}

extern "C" void thetaadvecuw1hoCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM wci,DECNUM *ee,DECNUM *ctheta,DECNUM * &thetaadvec)
{
	
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			DECNUM dxplus_i = 1.0f/dx;
			DECNUM dxcent_i = 1.0f/(2.0f*dx);
			DECNUM tthetaadvec,cthetab;
			DECNUM costhet,sinthet;
			DECNUM arrint,arrmint,arrmaxt;
			
	
			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);
			if (ntheta>1)
			{
				for (int itheta = 0; itheta < ntheta; itheta++)
				{
					unsigned int Tminus = mminusC(itheta, ntheta);
					unsigned int Tplus = pplusC(itheta, ntheta);

					tthetaadvec = 0.0;


					cthetab = 0.5f*(ctheta[i + itheta*nx*ny] + ctheta[i + Tplus*nx*ny]);
					if (cthetab > 0.0f)
					{
						arrint = ee[i + itheta*nx*ny] * cthetab;
					}
					else
					{
						arrint = ee[i + Tplus*nx*ny] * cthetab;
					}



					cthetab = 0.5f*(ctheta[i + Tminus*nx*ny] + ctheta[i + itheta*nx*ny]);
					if (cthetab > 0.0f)
					{
						arrmint = ee[i + Tminus*nx*ny] * cthetab;
					}
					else
					{
						arrmint = ee[i + itheta*nx*ny] * cthetab;
					}



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
				thetaadvec[i + 0*nx*ny] = 0.0f;
			}
		}
	}
	
	
}

extern "C" void thetaadvecuw2hoCPU(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *ctheta, DECNUM * thetaadvec)
{

	for (int ix = 0; ix < nx; ix++)
	{
		for (int iy = 0; iy < ny; iy++)
		{
			unsigned int i = ix + iy*nx;


			DECNUM dxplus_i = 1.0f / dx;
			DECNUM dxcent_i = 1.0f / (2.0f*dx);
			DECNUM tthetaadvec, cthetab;
			DECNUM costhet, sinthet;
			DECNUM arrint, arrmint, eupw;


			unsigned int xminus = mminusC(ix, nx);
			unsigned int xplus = pplusC(ix, nx);
			unsigned int yminus = mminusC(iy, ny);
			unsigned int yplus = pplusC(iy, ny);
			if (ntheta > 1)
			{

				for (int itheta = 0; itheta < ntheta; itheta++)
				{
					unsigned int Tminus = mminusC(itheta, ntheta);
					unsigned int Tplus = pplusC(itheta, ntheta);
					unsigned int Tminus2 = mminus2C(itheta, ntheta);
					unsigned int Tplus2 = pplus2C(itheta, ntheta);


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


}

extern "C" void eulerupwindCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM wci,DECNUM *&ee,DECNUM *xadvec,DECNUM *yadvec,DECNUM * thetaadvec)
{

	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			for (int itheta=0; itheta<ntheta; itheta++)
			{
				ee[i+itheta*nx*ny]=ee[i+itheta*nx*ny]-dt*(xadvec[i+itheta*nx*ny]+yadvec[i+itheta*nx*ny]+thetaadvec[i+itheta*nx*ny]);
			}
		}
	}
	

}

extern "C" void rollerlatbndCPU(int nx,int ny,int ntheta,DECNUM eps,DECNUM *hh,DECNUM *&rr)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);
			
			DECNUM hhi;



			hhi=hh[i];
			DECNUM wet=0.0f;

			if (hhi>eps)
			{
				wet=1.0f;
			}
			
			for (int itheta=0; itheta<ntheta; itheta++)
			{
								
				if (iy==0)
				{
					rr[i+itheta*nx*ny]=rr[ix+yplus*nx+itheta*nx*ny]*wet;
				}
				if (iy==ny-1)
				{
					rr[i+itheta*nx*ny]=rr[ix+yminus*nx+itheta*nx*ny]*wet;
				}
				
			}
		}
	}
}

extern "C" void energyCPU(int nx,int ny,int ntheta,DECNUM * &ee,DECNUM * sigm)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			DECNUM eetmp;
			DECNUM ssigm;
			ssigm=sigm[i];

			for (int itheta=0; itheta<ntheta; itheta++)
			{
				eetmp = ee[i+itheta*nx*ny]*ssigm;
				
				ee[i+itheta*nx*ny]=max(eetmp,0.0f);
				
			}
		}
	}
}

extern "C" void energintCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM rho,DECNUM g,DECNUM gammax,DECNUM * &E,DECNUM * &H,DECNUM * hh,DECNUM * &ee)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
			DECNUM Etmp;
			DECNUM Htmp;
			DECNUM hhtmp;
			
			Etmp=0.0f;
			hhtmp=hh[i];
			for (int itheta=0; itheta<ntheta; itheta++)
			{
				Etmp=Etmp+ee[i+itheta*nx*ny];
			}
			Etmp=Etmp*dtheta;
			Htmp=sqrtf(Etmp/(rho*g/8.0f));//Hrms
			DECNUM idiv=max(1.0f,powf(Htmp/(gammax*hhtmp),2.0f));
			
			for (int itheta=0; itheta<ntheta; itheta++)
			{
				ee[i+itheta*nx*ny]=ee[i+itheta*nx*ny]/idiv;
			}
			Htmp=min(Htmp,gammax*hhtmp);
			E[i]=(rho*g*Htmp*Htmp)/8.0f;
			H[i]=Htmp;
		}
	}

	
}

extern "C" void roelvinkCPU(int nx, int ny,DECNUM rho,DECNUM g,DECNUM gamma,DECNUM alpha,DECNUM n,DECNUM Trep,DECNUM * fwm,DECNUM * cfm,DECNUM *hh,DECNUM *H,DECNUM *E,DECNUM *&D, DECNUM *k)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			DECNUM fac,hroelvink,Qb;
			DECNUM Htmp,etmp;
			//int bk=3; //1 or 3
			
			fac=8.0/rho/g;
			hroelvink=hh[i];
			etmp=E[i];

			Htmp = sqrt(etmp *fac);


			Qb=1.0-exp(max(-1.0*pow((Htmp/gamma/hroelvink),n),-100.0));

			Qb=min(Qb,1.0f);


 			D[i]=Qb*2.0*alpha/Trep*etmp*Htmp/hroelvink;
			
			DECNUM fw;
			DECNUM sqrt2 = sqrt(2);// 1.4142136;
			DECNUM urms=pi*Htmp/(Trep*sqrt2*sinh(k[i]*hroelvink));
			DECNUM uorb=pi*Htmp/(Trep*sinhf(min(k[i]*hroelvink,10.0f)));
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
			fw=fwm[i];
			D[i]=D[i]+0.66666666f/pi*rho*fw*uorb*uorb*uorb;
		}
	}

}

extern "C" void baldockCPU(int nx, int ny,DECNUM rho,DECNUM g,DECNUM gamma,DECNUM alpha,DECNUM n,DECNUM Trep,DECNUM *fwm,DECNUM * cfm,DECNUM *hh,DECNUM *H,DECNUM *E,DECNUM *&D, DECNUM *k)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			DECNUM fac,hbald,Qb;
			DECNUM Htmp,etmp;
			DECNUM kh,Hb,R,ki;	
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
			
			// DECNUM fw;
			//DECNUM sqrt2=1.4142136f;
			//DECNUM urms=pi*Htmp/(Trep*sqrt2*sinh(min(max(k[i],0.01f)*hbald,10.0f)));
			//DECNUM uorb = pi*Htmp / (Trep*sinhf(min(max(k[i],0.01f)*hbald,10.0f)));
			DECNUM uorb = pi*Htmp / (Trep*sinhf(min(k[i] * hbald, 10.0f)));
			//DECNUM urms=uorb/sqrt2;
			D[i]=D[i]+0.66666666f/pi*rho*fwm[i]*uorb*uorb*uorb;
		}
	}
			
	
}

extern "C" void dissipationCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM eps,DECNUM dt,DECNUM g,DECNUM beta,DECNUM * wci,DECNUM *hh,DECNUM *&ee,DECNUM *D,DECNUM *E,DECNUM *&rr,DECNUM *c,DECNUM *cxsth,DECNUM *sxnth,DECNUM * uu,DECNUM * vv,DECNUM *&DR,DECNUM *&R)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
			DECNUM dd;
			DECNUM RR=0.0f;
			DECNUM drr=0.0f;
			DECNUM cx,cy;
			//DECNUM wci=0.0f;
			DECNUM rrd;
			DECNUM eet,rrt;

			DECNUM cc;
			DECNUM uui;
			DECNUM vvi;

			cc=c[i];
			uui=uu[i]*wci[i];
			vvi=vv[i]*wci[i];
			

			for (int itheta=0; itheta<ntheta; itheta++)
			{
					cx=cc*cxsth[itheta]+uui;
					cy=cc*sxnth[itheta]+vvi;

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
	}
}


extern "C" void meandirCPU(int nx,int ny,int ntheta,DECNUM rho,DECNUM g,DECNUM dtheta,DECNUM * ee,DECNUM * thet, DECNUM * &thetamean, DECNUM * &E,DECNUM * &H)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
			DECNUM sumethet=0;
			
			DECNUM sume=0;
			for (int itheta=0; itheta<ntheta; itheta++)
			{
				sume=sume+ee[i+itheta*nx*ny];
				sumethet=sumethet+ee[i+itheta*nx*ny]*thet[itheta];
			}
			sume=max(sume,0.00001f);
			
			thetamean[i]=(sumethet/ntheta)/(sume/ntheta);
			E[i]=sume*dtheta;

			H[i]=sqrtf(sume*dtheta/(rho*g/8));//sqrt(E[i]/(1/8*rho*g));
		}
	}
}

extern "C" void radstressCPU(int nx,int ny, int ntheta,DECNUM dx,DECNUM dtheta,DECNUM * ee,DECNUM *rr,DECNUM * cxsth,DECNUM * sxnth,DECNUM * cg,DECNUM * c,DECNUM * &Sxx,DECNUM * &Sxy,DECNUM * &Syy)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;

			DECNUM n;
			 

			n=cg[i]/c[i];
			
			DECNUM sume=0.0f;
			DECNUM tmpxx=0.0f;
			DECNUM tmpxy=0.0f;
			DECNUM tmpyy=0.0f;
			DECNUM rolxx=0.0f;
			DECNUM rolxy=0.0f;
			DECNUM rolyy=0.0f;

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
			

			Sxx[i]=(n*tmpxx-0.5f*sume)*dtheta+rolxx*dtheta;
			Syy[i]=(n*tmpyy-0.5f*sume)*dtheta+rolyy*dtheta;
			Sxy[i]=n*tmpxy*dtheta+rolxy*dtheta;
		}
	}
}

extern "C" void wavforceCPU(int nx,int ny, int ntheta,DECNUM dx,DECNUM dtheta,DECNUM * Sxx,DECNUM * Sxy,DECNUM * Syy,DECNUM * &Fx,DECNUM * &Fy,DECNUM * hh)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);
			
			//DECNUM hmin=0.000002f;
			
			DECNUM  Sxxi;
			DECNUM  Sxxr;
			DECNUM  Sxyi;
			DECNUM  Sxyt;
			DECNUM  Sxytr;
			DECNUM  Sxyb;
			DECNUM  Sxybr;
			
			DECNUM  Syyi;
			DECNUM  Syyt;
			DECNUM  Sxyr;
			DECNUM  Sxyl;
			DECNUM  Sxytl;

			DECNUM FFx;
			DECNUM FFy;
			
			


			Sxxi=Sxx[i];
			Sxxr=Sxx[xplus+iy*nx];
			Sxyi=Sxy[i];
			Sxyt=Sxy[ix+yplus*nx];
			Sxytr=Sxy[xplus+yplus*nx];
			Sxyb=Sxy[ix+yminus*nx];
			Sxybr=Sxy[xplus+yminus*nx];

			Syyi=Syy[i];
			Syyt=Syy[ix+yplus*nx];
			Sxyr=Sxy[xplus+iy*nx];
			Sxyl=Sxy[xminus+iy*nx];
			Sxytl=Sxy[xminus+yplus*nx];

			FFx=0.0f;
			FFy=0.0f;


			FFx=-1*(Sxxr-Sxxi)/dx-0.5f*(Sxyt+Sxytr-Sxyb-Sxybr)/(2.0f*dx);
				
			FFy=-1*(Syyt-Syyi)/dx-0.5f*(Sxyr+Sxytr-Sxyl-Sxytl)/(2.0f*dx);
			
			Fx[i]=FFx;
			Fy[i]=FFy;
		}
	}
				

}
extern "C" void twodimbndnoixCPU(int nx,int ny,DECNUM eps,DECNUM * hh,DECNUM * &F)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);
			
			DECNUM hhi;
			DECNUM Fi;
			DECNUM Ft;
			DECNUM Fb;
			DECNUM Fr;

			
			hhi=hh[i];
			DECNUM wet=0.0f;

			if (hhi>eps)
			{
				wet=1.0f;
			}
			

				Fi=F[i];
				Ft=F[ix+yplus*nx];
				Fb=F[ix+yminus*nx];
				Fr=F[xplus+iy*nx];

				//F[i]=Fi[tx][ty]*wet;
				if (iy==0)
				{
					F[i]=Ft*wet;
				}
				if (iy==ny-1)
				{
					F[i]=Fb*wet;
				}
				//if (ix==0)
				//{
				//	F[i]=Fr[tx][ty]*wet;
				//}
		}
	}

}

extern "C" void twodimbndCPU(int nx,int ny,DECNUM eps,DECNUM * hh,DECNUM * &F)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);
			
			DECNUM hhi;
			DECNUM Fi;
			DECNUM Ft;
			DECNUM Fb;
			DECNUM Fr;

			
			


			hhi=hh[i];
			DECNUM wet=0.0f;

			if (hhi>eps)
			{
				wet=1.0f;
			}
			

				Fi=F[i];
				Ft=F[ix+yplus*nx];
				Fb=F[ix+yminus*nx];
				Fr=F[xplus+iy*nx];

				//F[i]=Fi[tx][ty]*wet;
				if (iy==0)
				{
					F[i]=Ft*wet;
				}
				if (iy==ny-1)
				{
					F[i]=Fb*wet;
				}
				if (ix==0)
				{
					F[i]=Fr*wet;
				}
		}
	}

}


extern "C" void breakerdelayCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM g, DECNUM rho,DECNUM Trep, DECNUM eps,DECNUM * &urms,DECNUM *&ust,DECNUM *H,DECNUM *E,DECNUM *c,DECNUM *k, DECNUM *hh,DECNUM *rr)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
			DECNUM hmin=0.2;
			DECNUM delta=0.0f;
			DECNUM Hi=H[i];
			DECNUM ki=k[i];
			DECNUM hhi=max(hh[i],hmin);
			DECNUM ci=max(c[i],sqrt(hmin*g));
			DECNUM R;
			DECNUM ustw,uwf,vwf,ustr,usd,uorb;
			
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
	}

}

extern "C" void addavg_varCPU(int nx, int ny,DECNUM * &Varmean,DECNUM * Var)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			

			DECNUM mvari;
			DECNUM vari;

			mvari=Varmean[i];
			vari=Var[i];

			Varmean[i]=mvari+vari;
		}
	}


}

extern "C" void divavg_varCPU(int nx, int ny,DECNUM ntdiv,DECNUM * &Varmean)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			
			DECNUM mvari;
			mvari=Varmean[i];
			Varmean[i]=mvari/ntdiv;
		}
	}
	

}

extern "C" void resetavg_varCPU(int nx, int ny,DECNUM * &Varmean)
{
	for (int ix=0; ix<nx; ix++)
	{
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			Varmean[i]=0.0f;
		}
	}
}