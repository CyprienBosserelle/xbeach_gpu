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
#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>

#define pi 3.14159265
using DECNUM = float;


extern "C" void readXbbndhead(char * wavebnd, DECNUM &thetamin, DECNUM &thetamax, DECNUM &dtheta, DECNUM &dtwavbnd, int &nwavbnd, int &nwavfile)
{
	FILE * fwav;
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd, &nwavfile);
	fclose(fwav);
}
extern "C" void readbndhead(char * wavebnd, DECNUM &thetamin, DECNUM &thetamax, DECNUM &dtheta, DECNUM &dtwavbnd, int &nwavbnd)
{
	FILE * fwav;
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd);
	printf("rtwavbnd=%f\tnwavbnd=%d\n", dtwavbnd, nwavbnd);

	fclose(fwav);
}

extern "C" void readXbbndstep(int nx, int ny, int ntheta, char * wavebnd, int step, DECNUM &Trep, double *&qfile, double *&Stfile)
{
	FILE * fwav;
	FILE * fXq, *fXE;
	char Xbqfile[256];
	char XbEfile[256];
	double dummy;
	size_t result;
	DECNUM thetamin, thetamax, dtheta, dtwavbnd;
	int nwavbnd, nwavfile;

	printf("Reading next bnd file... ");
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd, &nwavfile);
	for (int n = 0; n < step; n++)
	{
		fscanf(fwav, "%f\t%s\t%s\n", &Trep, &Xbqfile, &XbEfile);
	}
	fclose(fwav);

	//printf("Xbq: %s\n",Xbqfile);
	//printf("Xbe: %s\n",XbEfile);

	fXq = fopen(Xbqfile, "rb");
	if (!fXq)
	{
		printf("Unable to open file %s\t", Xbqfile);
		return;
	}
	else
	{


		for (int nn = 0; nn < 4 * ny*nwavbnd; nn++)
		{
			result = fread(&dummy, sizeof(double), 1, fXq);
			qfile[nn] = (DECNUM)dummy;
		}

		fclose(fXq);
	}


	fXE = fopen(XbEfile, "rb");
	for (int nn = 0; nn < ntheta*ny*nwavbnd; nn++)
	{
		fread(&dummy, sizeof(double), 1, fXE);
		//printf("St=%f\n ",dummy);
		//Stfile[nn] = dummy;
		Stfile[nn] = (DECNUM)dummy;
	}
	fclose(fXE);

	printf("done \n");

}

extern "C" void readStatbnd(int nx, int ny, int ntheta, DECNUM rho, DECNUM g, char * wavebnd, double *&Tpfile, double *&Stfile)
{
	FILE * fwav;

	DECNUM dumDECNUM;
	size_t result;
	DECNUM thetamin, thetamax, dtheta, dtwavbnd;
	DECNUM Trepdum;
	int nwavbnd;

	printf("Reading bnd file... ");
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd);
	//printf("rtwavbnd=%f\n ",rtwavbnd);
	for (int ni = 0; ni < nwavbnd; ni++)
	{
		fscanf(fwav, "%f\t", &Trepdum);
		//printf("Tp=%f\n ",Trepdum);
		Tpfile[ni] = Trepdum;
		for (int i = 0; i < ntheta; i++)                             //! Fill St
		{
			fscanf(fwav, "%f\t", &dumDECNUM);
			//printf("St=%f\n ",dumDECNUM);
			for (int ii = 0; ii < ny; ii++)
			{
				Stfile[ii + i*ny + ni*ny*ntheta] = dumDECNUM;
			}
		}
	}
	fclose(fwav);
	printf("done \n");

}

/*
extern "C" void readbathy(int &nx, int &ny, float &dx, float &grdalpha, float *&zb)
{
	//read input data:
	printf("bathy: %s\n", filename);


	//read md file
	fid = fopen(filename, "r");
	fscanf(fid, "%u\t%u\t%f\t%*f\t%f", &nx, &ny, &dx, &grdalpha);
	grdalpha = grdalpha*pi / 180; // grid rotation

	int jread;
	//int jreadzs;
	for (int fnod = ny; fnod >= 1; fnod--)
	{

		fscanf(fid, "%u", &jread);
		
		for (int inod = 0; inod < nx; inod++)
		{
			fscanf(fid, "%f", &zb[inod + (jread - 1)*nx]);

		}
	}

	fclose(fid);
}
*/

