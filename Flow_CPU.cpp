#include <stdio.h>


#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>

#define pi 3.14159265
using DECNUM = float;

template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

extern "C" int mminusC(int ix,int nx)
{
	int xminus;
	if (ix<=0)
	{
		xminus=0;
	}	
	else
	{
		xminus=ix-1;	
	}
	return(xminus);
}
extern "C" int pplusC(int ix, int nx)
{
	int xplus;
	if (ix>=nx-1)
	{
		xplus=nx-1;
	}	
	else
	{
		xplus=ix+1;	
	}
	return(xplus);

}

extern "C" int mminus2C(int ix,int nx)
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
extern "C" int pplus2C(int ix, int nx)
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



extern "C" void ubndCPU(int nx, int ny, DECNUM dx, DECNUM dt,DECNUM g, DECNUM rho,DECNUM totaltime,DECNUM wavbndtime,DECNUM rt,DECNUM slbndtime, DECNUM rtsl,DECNUM zsbndold,DECNUM zsbndnew,DECNUM Trep,DECNUM * qbndold, DECNUM * qbndnew,DECNUM *&zs, DECNUM * &uu,DECNUM * &vv, DECNUM *vu, DECNUM * umean, DECNUM * vmean,DECNUM * zb,DECNUM * cg,DECNUM * hum, DECNUM * zo, DECNUM *Fx,DECNUM *&hh)
{
		int ix=0;
		
		for (int iy=0; iy<ny; iy++)
		{
			unsigned int i=ix+iy*nx;
			unsigned int xminus=mminusC(ix,nx);
			unsigned int xplus=pplusC(ix,nx);
			unsigned int xplus2=pplus2C(ix,nx);
			unsigned int yminus=mminusC(iy,ny);
			unsigned int yplus=pplusC(iy,ny);

			DECNUM ui,vi,thetai,vert;
			DECNUM beta, betar,betat,betab,bnp1,bn;
			DECNUM ht,htr;
			DECNUM theta0=0.0f;
			DECNUM alpha2=-1.0f*theta0;
			DECNUM epsi=0.005; //Not used!
			DECNUM ur,uumean,vvmean,urr,alphanew;
			DECNUM dbetadx,dbetady,dvudy,dhdx;
			DECNUM qx,qy,zsbnd;
			DECNUM order=2.0f;
			DECNUM ccg=cg[i];
			DECNUM cats=4; // number of wave period to average the current from
			DECNUM factime=0.0f;// 1.0f/cats/Trep*dt;
			DECNUM taper=min(totaltime/100.0,1.0);

				

			qx=(qbndold[iy]+(totaltime-wavbndtime+rt)*(qbndnew[iy]-qbndold[iy])/rt)*taper;
			qy=(qbndold[iy+ny]+(totaltime-wavbndtime+rt)*(qbndnew[iy+ny]-qbndold[iy+ny])/rt)*taper;
			zsbnd=zsbndold+(totaltime-rtsl)*(zsbndnew-zsbndold)/(slbndtime-rtsl);
	
			ht=zsbnd+zb[i];
			htr=zsbnd+zb[xplus+iy*nx];
			ui=qx/ht;
			vi=qy/ht;
			beta=uu[i]-2.0f*sqrt(g*hum[i]);
			betar=uu[xplus+iy*nx]-2.0f*sqrtf(g*hum[xplus+iy*nx]);
			betat=uu[ix+yplus*nx]-2.0f*sqrtf(g*hum[ix+yplus*nx]);
			betab=uu[ix+yminus*nx]-2.0f*sqrtf(g*hum[ix+yminus*nx]);
	
			dvudy=(vu[ix+(yminus)*nx]-vu[ix+(yplus)*nx])/(2.0f*dx);
			dbetadx=(betar-beta)/dx;
			dbetady=(betat-betab)/(2.0f*dx);

			dhdx=(htr-ht)/dx;

			bn=-1.0f*(uu[i]-sqrt(g*hum[i]))*dbetadx-vu[i]*dbetady+sqrtf(g*hum[i])*dvudy+1/rho*Fx[i]/hum[i]-zo[i]*sqrtf(uu[i]*uu[i]+vu[i]*vu[i])*uu[i]/hum[i]+g*dhdx;
			bnp1=beta+bn*dt;
	
			//WARNING this should be very inefficient. Need to find a better way. possibly inside another kernel
			// not neededd when epsi ==0.0...or factime==0.0
			DECNUM uumm=0.0f;
			DECNUM vvmm=0.0f;
			/*for (int jj=0; jj<ny; jj++)
			{
				uumm=uumm+uu[ix+jj*nx];
				vvmm=vvmm+vv[ix+jj*nx];
			}*/



			uumean=factime*uumm+umean[iy]*(1-factime);
			vvmean=factime*vvmm+vmean[iy]*(1-factime);
			umean[iy]=uumean;
			vmean[iy]=vvmean;
	
	


			thetai=atanf(vi/(ui+0.0000001f));
	
			vert=vu[i]-vvmean-vi;
	
			urr=(bnp1-uumean+2.0f*sqrtf(g*0.5f*(ht+htr))-ui*(ccg*(cosf(thetai))-sqrtf(g*0.5f*(ht+htr)))/(ccg*cosf(thetai)));
	
			for (int jj=0; jj<50; jj++)
			{
				ur=cosf(alpha2)/(cosf(alpha2)+1.0f)*urr;
				/*if(ur==0.0f)
				{
					ur=0.0000001f;
				}*/
				alphanew=atanf(vert/(ur+0.0000001f));
				if (alphanew>pi*0.5f)
				{
					alphanew=alphanew-pi;
				}
				if (alphanew<=-0.5f*pi)
				{
					alphanew=alphanew+pi;
				}
				if(abs(alphanew-alpha2)<0.001f)
				{
					break;
				}
				alpha2=alphanew;
			}


			//
			uu[i]=(order-1.0f)*ui+ur+uumean;//2.0f*ui-(sqrtf(g/(zs[i]+zb[i]))*(zs[i]-zsbnd));;//
			zs[i]=1.5f*((bnp1-uu[i])*(bnp1-uu[i])/(4.0f*g)-0.5f*(zb[i]+zb[xplus+iy*nx]))-0.5f*((betar-uu[xplus+iy*nx])*(betar-uu[xplus+iy*nx])/(4.0f*g)-0.5f*(zb[xplus+iy*nx]+zb[xplus2+iy*nx]));
			////
			//zsbnd+qx/(dx*dx)*dt;//
	
			hh[i]=zs[i]+zb[i];
			vv[i]=vv[xplus+iy*nx];
		}
	
	



}