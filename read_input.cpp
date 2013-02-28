#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>

#define pi 3.14159265



extern "C" void readXbbndhead(char * wavebnd,float &thetamin,float &thetamax,float &dtheta,float &dtwavbnd,int &nwavbnd,int &nwavfile)
{
	FILE * fwav;
	fwav=fopen(wavebnd,"r");
	fscanf(fwav,"%f\t%f\t%f\t%f\t%d\t%d",&thetamin,&thetamax,&dtheta,&dtwavbnd,&nwavbnd,&nwavfile);
	fclose(fwav);
}
extern "C" void readbndhead(char * wavebnd,float &thetamin,float &thetamax,float &dtheta,float &dtwavbnd,int &nwavbnd)
{
	FILE * fwav;
	fwav=fopen(wavebnd,"r");
	fscanf(fwav,"%f\t%f\t%f\t%f\t%d",&thetamin,&thetamax,&dtheta,&dtwavbnd,&nwavbnd);
	printf("rtwavbnd=%f\tnwavbnd=%d\n",dtwavbnd,nwavbnd);
	
	fclose(fwav);
}

extern "C" void readXbbndstep(int nx, int ny,int ntheta,char * wavebnd,int step,float &Trep,double *&qfile,double *&Stfile )
{
	FILE * fwav;
	FILE * fXq,* fXE;
	char Xbqfile[256];
	char XbEfile[256];
	double dummy;
	size_t result;
	float thetamin,thetamax,dtheta,dtwavbnd;
	int nwavbnd,nwavfile;

	printf("Reading next bnd file... ");
	fwav=fopen(wavebnd,"r");
	fscanf(fwav,"%f\t%f\t%f\t%f\t%d\t%d",&thetamin,&thetamax,&dtheta,&dtwavbnd,&nwavbnd,&nwavfile);
	for (int n=0; n<step; n++)
	{
		fscanf(fwav,"%f\t%s\t%s\n",&Trep,&Xbqfile,&XbEfile);
	}
	fclose(fwav);

	//printf("Xbq: %s\n",Xbqfile);
		//printf("Xbe: %s\n",XbEfile);
		
		fXq=fopen(Xbqfile,"rb");
		if (!fXq)
		{
			printf("Unable to open file %s\t", Xbqfile);
			return;
		}
		else
		{
			
			
			for (int nn=0; nn<3*ny*nwavbnd; nn++)
			{
				result=fread (&dummy,sizeof(double),1,fXq);
				qfile[nn]=dummy;
			}
			
			fclose(fXq);
		}


		fXE=fopen(XbEfile,"rb");
		for (int nn=0; nn<ntheta*ny*nwavbnd; nn++)
		{
			fread (&dummy,sizeof(double),1,fXE);
			//printf("St=%f\n ",dummy);
			Stfile[nn]=dummy;
		}
		fclose(fXE);

		printf("done \n");
	
}

extern "C" void readStatbnd(int nx, int ny,int ntheta,float rho,float g,char * wavebnd,double *&Tpfile,double *&Stfile )
{
	FILE * fwav;
	
	float dumfloat;
	size_t result;
	float thetamin,thetamax,dtheta,dtwavbnd;
	float Trepdum;
	int nwavbnd;

	printf("Reading bnd file... ");
	fwav=fopen(wavebnd,"r");
	fscanf(fwav,"%f\t%f\t%f\t%f\t%d",&thetamin,&thetamax,&dtheta,&dtwavbnd,&nwavbnd);
	//printf("rtwavbnd=%f\n ",rtwavbnd);
	for (int ni=0; ni<nwavbnd; ni++)
		{
			fscanf(fwav,"%f\t",&Trepdum);
			//printf("Tp=%f\n ",Trepdum);
			Tpfile[ni]=Trepdum;
			for (int i=0; i<ntheta; i++)                             //! Fill St
			{
				fscanf(fwav,"%f\t",&dumfloat);
				//printf("St=%f\n ",dumfloat);
				for (int ii=0; ii<ny; ii++)
				{
					Stfile[ii+i*ny+ni*ny*ntheta]=dumfloat;
				}
			}
		}
		fclose(fwav);
		printf("done \n");

}




