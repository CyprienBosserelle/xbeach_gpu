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



__global__ void ubnd(int nx, int ny, float dx, float dt,float g, float rho,float totaltime,float wavbndtime,float rt,float slbndtime, float rtsl,float zsbndold,float zsbndnew,float Trep,float * qbndold, float * qbndnew,float *zs, float * uu,float * vv, float *vu, float * umean, float * vmean,float * zb,float * cg,float * hum, float * zo, float *Fx,float *hh)
{	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int xplus2=pplus2(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	float ui,vi,thetai,vert;
	float beta, betar,betat,bnp1,bn;
	float ht,htr;
	float theta0=0.0f;
	float alpha2=-1*theta0;
	float epsi=0.005; //??
	float ur,uumean,vvmean,urr,alphanew;
	float dbetadx,dbetady,dvudy,dhdx;
	float qx,qy,zsbnd;
	float order=1;
	float ccg=cg[i];
	float cats=12;//4; // number of wave period to average the current from
	float factime= 1.0f/cats/Trep*dt;


	if (ix==0)
	{

	qx=0.0f;//(qbndold[iy]+(totaltime-wavbndtime+rt)*(qbndnew[iy]-qbndold[iy])/rt);
	qy=0.0f;//(qbndold[iy+ny]+(totaltime-wavbndtime+rt)*(qbndnew[iy+ny]-qbndold[iy+ny])/rt);
	zsbnd=zsbndold+(totaltime-rtsl)*(zsbndnew-zsbndold)/(slbndtime-rtsl);
	
	ht=zsbnd+zb[i];
	htr=zsbnd+zb[xplus+iy*nx];
	ui=qx/ht;
	vi=qy/ht;
	beta=uu[i]-2*sqrt(g*hum[i]);
	betar=uu[xplus+iy*nx]-2*sqrt(g*hum[xplus+iy*nx]);
	betat=uu[ix+yplus*nx]-2*sqrt(g*hum[ix+yplus*nx]);
	
	dvudy=(vu[ix+(yminus)*nx]-vu[ix+(yplus)*nx])/(2*dx);
    dbetadx=(betar-beta)/dx;
    dbetady=(betat-beta)/dx;

	dhdx=(htr-ht)/dx;

	bn=-1*(uu[i]-sqrt(g*hum[i]))*dbetadx-vu[i]*dbetady+sqrt(g*hum[i])*dvudy+1/rho*Fx[i]/hum[i]-zo[i]*sqrt(uu[i]*uu[i]+vu[i]*vu[i])*uu[i]/hum[i]+g*dhdx;
	bnp1=beta+bn*dt;
	uumean=factime*uu[i];//+umean[iy]*(1-factime);
	vvmean=factime*vv[i];//+vmean[iy]*(1-factime);
	umean[iy]=uumean;
	vmean[iy]=vvmean;
	if (ui==0.0f)
	{
		ui=0.00001f;
	}
	


	thetai=atanf(vi/ui);
	
	vert=vu[i]-vvmean-vi;
	
	urr=(bnp1-uumean+2*sqrtf(g*0.5*(ht+htr))-ui*(ccg*(cosf(thetai))-sqrtf(g*0.5*(ht+htr)))/(ccg*cosf(thetai)));
	
	for (int jj=0; jj<50; jj++)
	{
		ur=cosf(alpha2)/(cosf(alpha2)+1)*urr;
		alphanew=atanf(vert/max(ur,0.00000001f));
		if (alphanew>pi*0.5)
		{
			alphanew=alphanew-pi;
		}
		if (alphanew<-0.5*pi)
		{
			alphanew=alphanew+pi;
		}
		if((alphanew-alpha2)<0.001)
		{
			break;
		}
		alpha2=alphanew;
	}


    //
	uu[i]=uumean;//(order-1)*ui+ur+uumean;//2.0f*ui-(sqrtf(g/(zs[i]+zb[i]))*(zs[i]-zsbnd));;//
	zs[i]=zsbnd;//1.5*((bnp1-uu[i])*(bnp1-uu[i])/(4*g)-0.5*(zb[i]+zb[xplus+iy*nx]))-0.5*((betar-uu[xplus+iy*nx])*(betar-uu[xplus+iy*nx])/(4*g)-0.5*(zb[xplus+iy*nx]+zb[xplus2+iy*nx]));
	////
	//zsbnd+qx/(dx*dx)/dt;//hh[i]=zsbnd+zb[i];
	vv[i]=vv[xplus+iy*nx];
	}
	

	

	__syncthreads;



}
__global__ void wlevslopes(int nx, int ny,float dx,float eps,float *zs,float * dzsdx,float *dzsdy,float*hh)
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

	__shared__ float zsi[16][16];
	__shared__ float zsr[16][16];
	__shared__ float zst[16][16];
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
	zsi[tx][ty]=zs[i];
	zsr[tx][ty]=zs[xplus+iy*nx];
	zst[tx][ty]=zs[ix+yplus*nx];

	   
	//dzsdx[i]=(zs[ix+1+iy*nx]-zs[ix-1+iy*nx])/(2*dx);
	dzsdx[i]=(zsr[tx][ty]-zsi[tx][ty])/dx;//*whi*whr;
	dzsdy[i]=(zst[tx][ty]-zsi[tx][ty])/dx;//*whi*wht;
        
        	
}
__global__ void calcuvvu(int nx,int ny,float dx,float *uu,float *vv,float *vu,float *uv,float * ust,float *thetamean,float *ueu_g,float *vev_g,float *vmageu,float *vmagev,int* wetu, int* wetv)
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

	float vsu,usu,ueu,veu,usv,vsv,vev,uev;

	__shared__ float usti[16][16];
	__shared__ float ustr[16][16];
	__shared__ float ustt[16][16];

	__shared__ float tmeani[16][16];
	__shared__ float tmeanr[16][16];
	__shared__ float tmeant[16][16];
	__shared__ int wetui[16][16];
	__shared__ int wetvi[16][16];

	usti[tx][ty]=ust[i];
	tmeani[tx][ty]=thetamean[i];
	ustr[tx][ty]=ust[xplus+iy*nx];
	tmeanr[tx][ty]=thetamean[xplus+iy*nx];
	ustt[tx][ty]=ust[ix+yplus*nx];
	tmeant[tx][ty]=thetamean[ix+yplus*nx];
	wetui[tx][ty]=wetu[i];
	wetvi[tx][ty]=wetv[i];


	// V-velocities at u-points

    vu[i]=0.25f*(vv[ix+yminus*nx]+vv[ix+iy*nx]+vv[xplus+yminus*nx]+vv[xplus+iy*nx])*wetui[tx][ty];

	// U-velocities at v-points
	uv[i]=0.25f*(uu[xminus+iy*nx]+uu[ix+iy*nx]+uu[xminus+yplus*nx]+uu[ix+yplus*nx])*wetvi[tx][ty];
	

	
	//Calculate V-stokes at u points
	vsu=0.5f*(usti[tx][ty]*sinf(tmeani[tx][ty])+ustr[tx][ty]*sinf(tmeanr[tx][ty]))*wetui[tx][ty];
	//Calculate U-stokes at u points
	usu=0.5f*(usti[tx][ty]*cosf(tmeani[tx][ty])+ustr[tx][ty]*cosf(tmeanr[tx][ty]))*wetui[tx][ty];
	//Calculate U-euler at u points
	ueu=uu[i]-usu;
	//Calculate V-euler at u points
	veu=vu[i]-vsu;
	vmageu[i]=sqrtf(ueu*ueu+veu*veu);
	ueu_g[i]=ueu;


	usv=0.5f*(usti[tx][ty]*cosf(tmeani[tx][ty])+ustt[tx][ty]*cosf(tmeant[tx][ty]))*wetvi[tx][ty];
	vsv=0.5f*(usti[tx][ty]*sinf(tmeani[tx][ty])+ustt[tx][ty]*sinf(tmeant[tx][ty]))*wetvi[tx][ty];
	vev = vv[i] - vsv;
  	uev = uv[i] - usv;
	vmagev[i]=sqrtf(uev*uev+vev*vev);
	vev_g[i]=vev;

	


}




__global__ void udepthmomcont(int nx,int ny,float dx,float eps,float ummn,int* wetu,float * zs,float * uu,float * hh,float *hum, float *hu,float * zb)
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


	float humi;
	float hui;
	
	__shared__ float  hhi[16][16];
	__shared__ float  hhip[16][16];
	//__shared__ float  hhjp[4][4];
	
	hui=hu[i];
	hhi[tx][ty]=hh[i];
	hhip[tx][ty]=hh[xplus+iy*nx];
	//hhjp[tx][ty]=hh[ix+(iy+1)*nx];
	__syncthreads;
	
	//Water depth in u-points do momentum equation: mean
    humi=0.5f*(hhi[tx][ty]+hhip[tx][ty]);
    // Water depth in u-points do continuity equation: upwind
    
	



	__syncthreads;
	if (hui>eps && humi>eps)
	{
		wetu[i]=1;
	}
	else
	{
		wetu[i]=0;
	}
	
	
	hum[i]=max(humi,eps);

}

__global__ void vdepthmomcont(int nx,int ny,float dx,float eps,float ummn,int* wetv,float * zs,float * vv,float * hh,float *hvm, float *hv,float * zb)
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

	float hvmi,hvi;
	
	__shared__ float  hhi[16][16];
	//__shared__ float  hhip[4][4];
	__shared__ float  hhjp[16][16];
	
	
	hhi[tx][ty]=hh[i];
	//hhip[tx][ty]=hh[ix+1+iy*nx];
	hhjp[tx][ty]=hh[ix+yplus*nx];
	__syncthreads;

	
	//Water depth in u-points do momentum equation: mean
    //hvmi=max(0.5f*(hh[i]+hh[ix+(min(iy,ny-2)+1)*nx]),eps);
    // Water depth in u-points do continuity equation: upwind

	//hvi=0.5f*(hhjp[tx][ty]-hhi[tx][ty])+hhi[tx][ty];
	hvmi=0.5f*(hhi[tx][ty]+hhjp[tx][ty]);
	//hvm(i,j)=max(.5d0*(hh(i,j)+hh(i,min(ny,j)+1)),par%eps)  
	
	
	hvi=hv[i];
	

	if (hvi>eps && hvmi>eps)
	{
		wetv[i]=1;
	}
	else
	{
		wetv[i]=0;
	}
    hvm[i]=max(hvmi,eps);
}

__global__ void depthhu(int nx,int ny,float dx,float ummn,float eps,float *hh,float * uu,float * hu,float *zs,float *zb)
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
	

	float hui=0.0f;
	
	__shared__ float  hhi[16][16];
	__shared__ float  hhip[16][16];
	hhi[tx][ty]=hh[i];
	hhip[tx][ty]=hh[xplus+iy*nx];

	if (uu[i]>ummn)
	{
		//hui=hhi[tx][ty];
		hui=zs[i]-max(-1.0f*zb[i],-1.0f*zb[xplus+iy*nx]);
	}
	else
	{
		if(uu[i]<-1.0f*ummn)
		{
			//hui=hhip[tx][ty];
			hui=zs[xplus+iy*nx]-max(-1.0f*zb[i],-1.0f*zb[xplus+iy*nx]);
		}
		else
		{
			hui=max(max(zs[i],zs[xplus+iy*nx])-max(-1.0f*zb[i],-1.0f*zb[xplus+iy*nx]),eps);
		}
				
	}
	//hui=0.5f*(hhip[tx][ty]+hhi[tx][ty]);
	hui=max(hui,0.0f);
	hu[i]=hui;

}

__global__ void depthhv(int nx,int ny,float dx,float ummn,float eps,float *hh,float * vv,float * hv,float *zs,float *zb)
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
	float hvi=0.0f;

	__shared__ float  hhi[16][16];
	//__shared__ float  hhip[4][4];
	__shared__ float  hhjp[16][16];
	
	
	hhi[tx][ty]=hh[i];
	//hhip[tx][ty]=hh[ix+1+iy*nx];
	hhjp[tx][ty]=hh[ix+yplus*nx];
	__syncthreads;
	if (vv[i]>ummn)
	{
		//hvi=hhi[tx][ty];
		hvi=zs[i]-max(-1.0f*zb[i],-1.0f*zb[ix+yplus*nx]);
	}
	else
	{
		if(vv[i]<-1*ummn)
		{
			//hvi=hhjp[tx][ty];
			hvi=zs[ix+yplus*nx]-max(-1.0f*zb[i],-1.0f*zb[ix+yplus*nx]);
		}
		else
		{
			hvi=max(max(zs[i],zs[ix+yplus*nx])-max(-1.0f*zb[i],-1.0f*zb[ix+yplus*nx]),eps);
			//hv[i]=hvm[i];
		}
	}
	hvi=max(hvi,0.0f);
	
	hv[i]=hvi;
}

__global__ void ududx_adv(int nx,int ny,float dx,float * hu,float * hum, float * uu, float * ududx)
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
	float qin,uududx;

	__shared__ float uui[16][16];
	__shared__ float uur[16][16];
	__shared__ float uul[16][16];


	__shared__ float hui[16][16];




	__shared__ float humi[16][16];


	uui[tx][ty]=uu[i];
	uur[tx][ty]=uu[xplus+iy*nx];
	uul[tx][ty]=uu[xminus+iy*nx];


	hui[tx][ty]=hu[i];
	humi[tx][ty]=hum[i];


	
	uududx=0.0f;
	qin=0.5f*(hui[tx][ty]*uui[tx][ty]+hu[xminus+iy*nx]*uul[tx][ty]);
	//ududx
	if (qin>0.0f)
	{
			uududx=uududx+qin/humi[tx][ty]*(uui[tx][ty]-uul[tx][ty])/dx;
	}
	qin=-0.5f*(hui[tx][ty]*uui[tx][ty]+hu[xplus+iy*nx]*uur[tx][ty]);
	if (qin>0.0f)
	{
			uududx=uududx+qin/humi[tx][ty]*(uui[tx][ty]-uur[tx][ty])/dx;
	}
	ududx[i]=uududx;
}


__global__ void ududx_adv2(int nx,int ny,float dx,float * hu,float * hum, float * uu, float * ududx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xminus2=mminus2(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int xplus2=pplus2(ix,nx);

	float qin,uududx;

	__shared__ float uui[16][16];
	__shared__ float uur[16][16];
	__shared__ float uur2[16][16];
	__shared__ float uul[16][16];
	__shared__ float uul2[16][16];

	__shared__ float hui[16][16];
	__shared__ float hur[16][16];
	__shared__ float hul[16][16];
	__shared__ float humi[16][16];


	uui[tx][ty]=uu[i];
	uur[tx][ty]=uu[xplus+iy*nx];
	uur2[tx][ty]=uu[xplus2+iy*nx];
	uul[tx][ty]=uu[xminus+iy*nx];
	uul2[tx][ty]=uu[xminus2+iy*nx];


	hui[tx][ty]=hu[i];
	hur[tx][ty]=hu[xplus+iy*nx];
	hul[tx][ty]=hu[xminus+iy*nx];
	humi[tx][ty]=hum[i];


	
	uududx=0.0f;
	qin=0.5f*(hui[tx][ty]*uui[tx][ty]+hul[tx][ty]*uul[tx][ty]);
	//ududx
	if (qin>0.0f)
	{
			uududx=uududx+qin/humi[tx][ty]*(3*uui[tx][ty]-4*uul[tx][ty]+uul2[tx][ty])/(2*dx);
	}
	qin=-0.5f*(hui[tx][ty]*uui[tx][ty]+hur[tx][ty]*uur[tx][ty]);
	if (qin>0.0f)
	{
			uududx=uududx+qin/humi[tx][ty]*(3*uui[tx][ty]-4*uur[tx][ty]+uur2[tx][ty])/(2*dx);
	}
	ududx[i]=uududx;
}



	__global__ void vdudy_adv(int nx,int ny,float dx,float * hv,float * hum, float * uu,float *vv, float * vdudy)
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
	float qin,vvdudy;

	__shared__ float uui[16][16];
	__shared__ float uut[16][16];
	__shared__ float uub[16][16];

	uui[tx][ty]=uu[i];
	uut[tx][ty]=uu[ix+yplus*nx];
	uub[tx][ty]=uu[ix+yminus*nx];

	vvdudy=0.0f;
	
	qin=0.5f*(vv[ix+yminus*nx]*hv[ix+yminus*nx]+vv[xplus+yminus*nx]*hv[xplus+yminus*nx]);
	if (qin>0.0f)
	{
			vvdudy=vvdudy+qin/hum[i]*(uui[tx][ty]-uub[tx][ty])/dx;
	}
	qin=-0.5f*(vv[i]*hv[i]+vv[xplus+iy*nx]*hv[xplus+iy*nx]);
	if (qin>0.0f)
	{
			vvdudy=vvdudy+qin/hum[i]*(uui[tx][ty]-uut[tx][ty])/dx;
	}
	vdudy[i]=vvdudy;

}

	__global__ void vdudy_adv2(int nx,int ny,float dx,float * hv,float * hum, float * uu,float *vv, float * vdudy)
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
	unsigned int yminus2=mminus2(iy,ny);
	unsigned int yplus2=pplus2(iy,ny);
	float qin,vvdudy;

	__shared__ float uui[16][16];
	__shared__ float uut[16][16];
	__shared__ float uut2[16][16];
	__shared__ float uub[16][16];
	__shared__ float uub2[16][16];

	__shared__ float vvi[16][16];
	__shared__ float vvr[16][16];
	__shared__ float vvb[16][16];
	__shared__ float vvbr[16][16];
	__shared__ float hvi[16][16];
	__shared__ float hvr[16][16];
	__shared__ float hvb[16][16];
	__shared__ float hvbr[16][16];
	__shared__ float humi[16][16];


	uui[tx][ty]=uu[i];
	uut[tx][ty]=uu[ix+yplus*nx];
	uub[tx][ty]=uu[ix+yminus*nx];
	uut2[tx][ty]=uu[ix+yplus2*nx];
	uub2[tx][ty]=uu[ix+yminus2*nx];

	vvi[tx][ty]=vv[i];
	vvr[tx][ty]=vv[xplus+iy*nx];
	vvb[tx][ty]=vv[ix+yminus*nx];
	vvbr[tx][ty]=vv[xplus+yminus*nx];
	hvi[tx][ty]=hv[i];
	hvr[tx][ty]=hv[xplus+iy*nx];
	hvb[tx][ty]=hv[ix+yminus*nx];
	hvbr[tx][ty]=hv[xplus+yminus*nx];
	humi[tx][ty]=hum[i];


	vvdudy=0.0f;
	
	qin=0.5f*(vvb[tx][ty]*hvb[tx][ty]+vvbr[tx][ty]*hvbr[tx][ty]);
	if (qin>0.0f)
	{
			vvdudy=vvdudy+qin/humi[tx][ty]*(3*uui[tx][ty]-4*uub[tx][ty]+uub2[tx][ty])/(2*dx);
	}
	qin=-0.5f*(vvi[tx][ty]*hvi[tx][ty]+vvr[tx][ty]*hvr[tx][ty]);
	if (qin>0.0f)
	{
			vvdudy=vvdudy+qin/humi[tx][ty]*(3*uui[tx][ty]-4*uut[tx][ty]+uut2[tx][ty])/(2*dx);
	}
	vdudy[i]=vvdudy;

}



__global__ void vdvdy_adv(int nx,int ny,float dx,float * hv,float * hvm, float * vv, float * vdvdy)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	float qin,vvdvdy;

	vvdvdy=0.0f;

	qin=0.5f*(vv[i]*hv[i]+vv[ix+yminus*nx]*hv[ix+yminus*nx]);
	if (qin>0.0f)
	{
		vvdvdy=vvdvdy+qin/hvm[i]*(vv[i]-vv[ix+(yminus)*nx])/(dx);
	}
	qin=-0.5f*(hv[i]*vv[i]+hv[ix+(yplus)*nx]*vv[ix+(yplus)*nx]);
	if (qin>0.0f)
	{
			vvdvdy=vvdvdy+qin/hvm[i]*(vv[i]-vv[ix+(yplus)*nx])/(dx);
	}
	vdvdy[i]=vvdvdy;
}

__global__ void vdvdy_adv2(int nx,int ny,float dx,float * hv,float * hvm, float * vv, float * vdvdy)
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
	unsigned int yminus2=mminus2(iy,ny);
	unsigned int yplus2=pplus2(iy,ny);

	__shared__ float vvi[16][16];
	__shared__ float vvb[16][16];
	__shared__ float vvb2[16][16];
	__shared__ float vvt[16][16];
	__shared__ float vvt2[16][16];
	__shared__ float hvi[16][16];
	__shared__ float hvb[16][16];
	__shared__ float hvt[16][16];
	__shared__ float hvmi[16][16];

	vvi[tx][ty]=vv[i];
	vvb[tx][ty]=vv[ix+yminus*nx];
	vvb2[tx][ty]=vv[ix+yminus2*nx];
	vvt[tx][ty]=vv[ix+yplus*nx];
	vvt2[tx][ty]=vv[ix+yplus2*nx];
	hvi[tx][ty]=hv[i];
	hvb[tx][ty]=hv[ix+yminus*nx];
	hvt[tx][ty]=hv[ix+yplus*nx];
	hvmi[tx][ty]=hvm[i];


	float qin,vvdvdy;

	vvdvdy=0.0f;

	qin=0.5*(vvi[tx][ty]*hvi[tx][ty]+vvb[tx][ty]*hvb[tx][ty]);
	if (qin>0.0f)
	{
		vvdvdy=vvdvdy+qin/hvmi[tx][ty]*(3.0f*vvi[tx][ty]-4.0f*vvb[tx][ty]+vvb2[tx][ty])/(2*dx);
	}
	qin=-0.5f*(hvi[tx][ty]*vvi[tx][ty]+hvt[tx][ty]*vvt[tx][ty]);
	if (qin>0.0f)
	{
			vvdvdy=vvdvdy+qin/hvmi[tx][ty]*(3.0f*vvi[tx][ty]-4.0f*vvt[tx][ty]+vvt2[tx][ty])/(2*dx);
	}
	vdvdy[i]=vvdvdy;
}

__global__ void udvdx_adv(int nx,int ny,float dx,float * hu,float * hvm,float * uu, float * vv, float * udvdx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	float qin,uudvdx;

	uudvdx=0.0f;
	qin=0.5*(uu[xminus+iy*nx]*hu[xminus+iy*nx]+uu[xminus+yplus*nx]*hu[xminus+yplus*nx]);
	if (qin>0.0f)
	{
			uudvdx=uudvdx+qin/hvm[i]*(vv[i]-vv[xminus+iy*nx])/(dx);
	}
	qin=-0.5*(uu[i]*hu[i]+uu[ix+yplus*nx]*hu[ix+yplus*nx]);
	if (qin>0.0f)
	{
			uudvdx=uudvdx+qin/hvm[i]*(vv[i]-vv[xplus+iy*nx])/(dx);
	}

	udvdx[i]=uudvdx;


}

__global__ void udvdx_adv2(int nx,int ny,float dx,float * hu,float * hvm,float * uu, float * vv, float * udvdx)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;

	unsigned int tx =threadIdx.x;
	unsigned int ty= threadIdx.y;

	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int xminus2=mminus2(ix,nx);
	unsigned int xplus2=pplus2(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	__shared__ float uui[16][16];
	__shared__ float uut[16][16];
	__shared__ float uul[16][16];
	__shared__ float uutl[16][16];
	__shared__ float vvi[16][16];
	__shared__ float vvl[16][16];
	__shared__ float vvl2[16][16];
	__shared__ float vvr[16][16];
	__shared__ float vvr2[16][16];
	__shared__ float hui[16][16];
	__shared__ float hut[16][16];
	__shared__ float hul[16][16];
	__shared__ float hutl[16][16];
	__shared__ float hvmi[16][16];

	uui[tx][ty]=uu[i];
	uut[tx][ty]=uu[ix+yplus*nx];
	uul[tx][ty]=uu[xminus+iy*nx];
	uutl[tx][ty]=uu[xminus+yplus*nx];
	vvi[tx][ty]=vv[i];
	vvl[tx][ty]=vv[xminus+iy*nx];
	vvl2[tx][ty]=vv[xminus2+iy*nx];
	vvr[tx][ty]=vv[xplus+iy*nx];
	vvr2[tx][ty]=vv[xplus2+iy*nx];
	hui[tx][ty]=hu[i];
	hut[tx][ty]=hu[ix+yplus*nx];
	hul[tx][ty]=hu[xminus+iy*nx];
	hutl[tx][ty]=hu[xminus+yplus*nx];
	hvmi[tx][ty]=hvm[i];


	float qin,uudvdx;

	uudvdx=0.0f;
	qin=0.5*(uul[tx][ty]*hul[tx][ty]+uutl[tx][ty]*hutl[tx][ty]);
	if (qin>0.0f)
	{
			uudvdx=uudvdx+qin/hvmi[tx][ty]*(3*vvi[tx][ty]-4*vvl[tx][ty]+vvl2[tx][ty])/(2*dx);
	}
	qin=-0.5*(uui[tx][ty]*hui[tx][ty]+uut[tx][ty]*hut[tx][ty]);
	if (qin>0.0f)
	{
			uudvdx=uudvdx+qin/hvmi[tx][ty]*(3*vvi[tx][ty]-4*vvr[tx][ty]+vvr2[tx][ty])/(2*dx);
	}

	udvdx[i]=uudvdx;


}


__global__ void smago(int nx,int ny,float dx,float * uu, float * vv,float nuh, float * nuhgrid,int usesmago)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	float dudx,dudy,dvdx,dvdy,tau;
	if (usesmago==1)
	{
		dudx=(uu[i]-uu[xminus+iy*nx])/dx;
		dudy=0.50f*(uu[ix+yplus*nx]-uu[ix+yminus*nx]+uu[xminus+yplus*nx]-uu[xminus+yminus*nx])/dx;
		dvdy=(vv[i]-vv[ix+yminus*nx])/dx;
		dvdx=0.50f*(vv[xplus+iy*nx]-vv[xminus+iy*nx]+vv[xplus+yminus*nx]-vv[xminus+yminus*nx])/dx;
		tau=sqrt(2.0f*dudx*dudx+2.0f*dvdy*dvdy+powf(dudy+dvdx,2.0f));
		nuhgrid[i]=nuh*nuh*dx*dx*tau;
	}
	else
	{
		nuhgrid[i]=nuh;
	}
	

}
	
__global__ void viscou(int nx,int ny,float dx,float rho,float eps,float nuhfac, float * nuhgrid,float *hh,float *hum,float *hvm,float * DR,float *uu,int * wetu,float * viscu)
{						
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	
	float nuh=nuhgrid[i];

	float nnuh,dudx1,dudx2,dudy1,dudy2;
	
	//if(ix>3)
	//{
		nnuh=max(nuh,nuhfac*hh[i]*powf(DR[i]/rho,1.0f/3.0f));
	//}
	//else
	//{
	//	nnuh=nuh*10;
	//}
	dudx1=hh[xplus+iy*nx]*(uu[xplus+iy*nx]-uu[i])/dx;
	dudx2=hh[i]*(uu[i]-uu[xminus+iy*nx])/dx;
	dudy1=0.5f*(hvm[i]+hvm[xplus+iy*nx])*(uu[ix+yplus*nx]-uu[i])/dx;
	dudy2=0.5f*(hvm[ix+yminus*nx]+hvm[xplus+yminus*nx])*(uu[i]-uu[ix+yminus*nx])/dx;
	viscu[i]=nnuh/hum[i]*((dudx1-dudx2)/(dx)*wetu[xplus+iy*nx]*wetu[xminus+iy*nx]+(dudy1-dudy2)/dx*wetu[ix+yplus*nx]*wetu[ix+yminus*nx]); 
	
//*wetu[xplus+iy*nx]*wetu[xplus+iy*nx]
}

__global__ void viscov(int nx,int ny,float dx,float rho,float eps,float nuhfac, float * nuhgrid,float *hh,float *hum,float *hvm,float * DR,float *vv,int * wetv,float * viscv)
{						
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	
	float nuh=nuhgrid[i];

	float nnuh,dvdx1,dvdx2,dvdy1,dvdy2;

	nnuh=max(nuh,nuhfac*hh[i]*powf(DR[i]/rho,1.0f/3.0f));
	dvdx1=0.5f*(hum[i]+hum[ix+yplus*nx])*(vv[xplus+iy*nx]-vv[i])/dx;
	dvdx2=0.5f*(hum[xminus+iy*nx]+hum[xminus+yplus*nx])*(vv[i]-vv[xminus+iy*nx])/dx;
	dvdy1=hh[ix+yplus*nx]*(vv[ix+yplus*nx]-vv[i])/dx;
	dvdy2=hh[i]*(vv[i]-vv[ix+yminus*nx])/dx;
	viscv[i]=nnuh/hvm[i]*((dvdx1-dvdx2)/(dx)*wetv[xplus+iy*nx]*wetv[xminus+iy*nx]+(dvdy1-dvdy2)/dx*wetv[ix+yplus*nx]*wetv[ix+yminus*nx]);

}

__global__ void viscovbnd(int nx,int ny,float * viscv )
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);
	
	if(iy==ny-1)
	{
		viscv[i]=viscv[ix+yminus*nx];
	}
	if(iy==0)
	{
		viscv[i]=viscv[ix+yplus*nx];
	}
}



__global__ void eulerustep(int nx,int ny,float dx,float dt,float g,float rho,float * zo,float fc,float windth,float windv,float Cd,float *uu,float * urms,float *ududx,float *vdudy,float *viscu,float *dzsdx,float *hu,float *hum,float *Fx,float *vu,float * ueu_g,float * vmageu,int *wetu)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i=ix+iy*nx;
	int tx =threadIdx.x;
	int ty= threadIdx.y;

	float ueu;
	float taubx;
	float hmin=1.0;
	unsigned int xminus=mminus(ix,nx);
	unsigned int xplus=pplus(ix,nx);
	unsigned int yminus=mminus(iy,ny);
	unsigned int yplus=pplus(iy,ny);

	

	__shared__ float  uui[16][16];
	
	uui[tx][ty]=uu[i];
	ueu=ueu_g[i];

	__syncthreads;

//&& ix>0
	if (wetu[i]==1 )
	{		
	taubx=zo[i]*rho*ueu*sqrtf(1.3456f*urms[i]*urms[i]+vmageu[i]*vmageu[i]);
		
	uui[tx][ty]=uui[tx][ty]-dt*(ududx[i]+vdudy[i]-viscu[i]+g*dzsdx[i]+taubx/(rho*hu[i])-Fx[i]/(rho*max(hum[i],hmin))-1.25f*Cd*cosf(windth)*windv*windv/(rho*hum[i])-fc*vu[i]);
	
	//viscu[i]=taubx;
	
	}
	else
	{
		uui[tx][ty]=0.0f;
		viscu[i]=0.0f;	

	}
	//if (ix>0)
	{
		uu[i]=uui[tx][ty];
		
	}


}

__global__ void eulervstep(int nx,int ny,float dx,float dt,float g,float rho,float * zo,float fc,float windth,float windv,float Cd,float *vv,float * urms,float *udvdx,float *vdvdy,float *viscv,float *dzsdy,float *hv,float *hvm,float *Fy,float *uv,float * vev_g,float * vmagev,int *wetv)
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

	__shared__ float  vvi[16][16];
	__shared__ float  urmsi[16][16];
	__shared__ float  vmagvi[16][16];
	__shared__ float  hvmi[16][16];

	float tauby,vev;

	float hmin=1.0;

	vvi[tx][ty]=vv[i];
	urmsi[tx][ty]=urms[i];
	vmagvi[tx][ty]=vmagev[i];
	hvmi[tx][ty]=hvm[i];
	
// && ix>0
	if (wetv[i]==1)
	{
	vev=vev_g[i];
	
	tauby=zo[i]*rho*vev*sqrtf(1.3456f*urmsi[tx][ty]*urmsi[tx][ty]+vmagvi[tx][ty]*vmagvi[tx][ty]);
	vvi[tx][ty]=vvi[tx][ty]-dt*(udvdx[i]+vdvdy[i]-viscv[i]+g*dzsdy[i]+tauby/(rho*hv[i])-Fy[i]/(rho*max(hvmi[tx][ty],hmin))+fc*uv[i]-1.25f*Cd*sinf(windth)*windv*windv/(rho*hvmi[tx][ty]));
	
	//viscv[i]=tauby;
	
	}
	else
	{
		vvi[tx][ty]=0.0f;
		viscv[i]=0.0f;
	}
	//if(ix>0 && iy>0 && iy<ny)
	{
		vv[i]=vvi[tx][ty];
	}//vdvdy[i]=tauby;

}




__global__ void continuity(int nx,int ny,float dx,float dt,float eps,float * uu,float* hu,float* vv,float* hv,float* zs,float *hh,float *zb,float * dzsdt)
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

	float qx,qy,qxm,qym,dzdt;
	float zz;
	
	__shared__ float uui[16][16];
	__shared__ float uul[16][16];
	__shared__ float vvi[16][16];
	__shared__ float vvb[16][16];
	__shared__ float hui[16][16];
	__shared__ float hul[16][16];
	__shared__ float hvi[16][16];
	__shared__ float hvb[16][16];
	
	uui[tx][ty]=uu[i];
	vvi[tx][ty]=vv[i];
	uul[tx][ty]=uu[xminus+iy*nx];
	vvb[tx][ty]=vv[ix+yminus*nx];
	hui[tx][ty]=hu[i];
	hul[tx][ty]=hu[xminus+iy*nx];
	hvi[tx][ty]=hv[i];
	hvb[tx][ty]=hv[ix+yminus*nx];
	
	
	
	
	zz=zs[i];

	qx=uui[tx][ty]*hui[tx][ty];
	qy=vvi[tx][ty]*hvi[tx][ty];
	
	qxm=uul[tx][ty]*hul[tx][ty];
	
	qym=vvb[tx][ty]*hvb[tx][ty];
	dzdt=(qxm-qx+qym-qy)/dx;
	
	
	__syncthreads;
	
	if(ix>0)
	{
	dzsdt[i]=dzdt;
	
	
	zs[i]=zz+dzdt*dt;
	
	//hh[i]=max(hh[i]+dzdt*dt,eps);
	}
	
}

__global__ void hsbnd(int nx,int ny,float eps,float * hh,float *zb,float *zs)
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
	
	
	__shared__ float Fi[16][16];
	__shared__ float Ft[16][16];
	__shared__ float Fb[16][16];
	__shared__ float Fr[16][16];

	
		Fi[tx][ty]=max(hh[i],eps);
		Ft[tx][ty]=max(hh[ix+yplus*nx],eps);
		Fb[tx][ty]=max(hh[ix+yminus*nx],eps);
		Fr[tx][ty]=max(hh[xplus+iy*nx],eps);

		//hh[i]=Fi[tx][ty];
		
		
		
		hh[i]=max(zb[i]+zs[i],eps);
		

}


__global__ void uuvvzslatbnd(int nx,int ny,float * uu,float * vv,float *zs)
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
	
	__shared__ float vvr[16][16];
	__shared__ float vvb[16][16];
	__shared__ float vvt[16][16];
	__shared__ float uut[16][16];
	__shared__ float uub[16][16];
	__shared__ float zst[16][16];
	__shared__ float zsb[16][16];
	__shared__ float zsl[16][16];
	
	

	
		
		uut[tx][ty]=uu[ix+yplus*nx];
		uub[tx][ty]=uu[ix+yminus*nx];
		vvr[tx][ty]=vv[xplus+iy*nx];
		vvt[tx][ty]=vv[ix+yplus*nx];
		vvb[tx][ty]=vv[ix+yminus*nx];
		zst[tx][ty]=zs[ix+yplus*nx];
		zsb[tx][ty]=zs[ix+yminus*nx];
		zsl[tx][ty]=zs[xminus+iy*nx];

		//F[i]=Fi[tx][ty]*wet;
		if (iy==0)
		{
			uu[i]=uut[tx][ty];
			vv[i]=vvt[tx][ty];
			zs[i]=zst[tx][ty];
		}
		if (iy==ny-1)
		{
			uu[i]=uub[tx][ty];
			vv[i]=vvb[tx][ty];
			zs[i]=zsb[tx][ty];
		}
		//if (iy==ny-2)
		//{
		//	vv[i]=0.0f;
		//}
		if (ix==0)
		{
			vv[i]=vvr[tx][ty];
		}
		if (ix==nx-1)
		{
			//zs[i]=zsl[tx][ty];
		}

}
