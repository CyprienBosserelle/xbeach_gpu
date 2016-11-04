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

#include "XBeachGPU.h"

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

XBGPUParam readparamstr(std::string line, XBGPUParam param)
{


	std::string parameterstr, parametervalue;

	///////////////////////////////////////////////////////
	// General parameters
	parameterstr = "bathy =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Bathymetryfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}
	
	//
	parameterstr = "depfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Bathymetryfile = parametervalue;
	}
		
	
	//
	parameterstr = "swave =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.swave = std::stoi(parametervalue);
	}

	//
	parameterstr = "flow =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.flow = std::stoi(parametervalue);
	}

	//
	parameterstr = "sedtrans =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.sedtrans = std::stoi(parametervalue);
	}

	//
	parameterstr = "morphology =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.morphology = std::stoi(parametervalue);
	}

	//
	parameterstr = "gpudevice =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.GPUDEVICE = std::stoi(parametervalue);
	}

	///////////////////////////////////////////////////////
	// Flow parameters
	//
	parameterstr = "eps =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.eps = std::stod(parametervalue);
	}
	
	parameterstr = "cf =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.cf = std::stod(parametervalue);
	}

	parameterstr = "cfsand =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.cfsand = std::stod(parametervalue);
	}

	parameterstr = "cfreef =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.cfreef = std::stod(parametervalue);
	}

	parameterstr = "nuh =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.nuh = std::stod(parametervalue);
	}

	parameterstr = "nuhfac =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.nuhfac = std::stod(parametervalue);
	}

	parameterstr = "usesmago =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.usesmago = std::stoi(parametervalue);
	}

	parameterstr = "smag =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.smag = std::stod(parametervalue);
	}

	parameterstr = "lat =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.lat = std::stod(parametervalue);
	}

	parameterstr = "Cd =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Cd = std::stod(parametervalue);
	}

	parameterstr = "wci =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wci = std::stod(parametervalue);
	}

	parameterstr = "hwci =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.hwci = std::stod(parametervalue);
	}

	///////////////////////////////////////////////////////
	// Wave parameters
	//
	parameterstr = "breakmodel =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.breakmodel = std::stoi(parametervalue);
	}

	parameterstr = "gamma =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.gammaa = std::stod(parametervalue);
	}

	parameterstr = "n =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.n = std::stod(parametervalue);
	}

	parameterstr = "alpha =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.alpha = std::stod(parametervalue);
	}

	parameterstr = "gammax =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.gammax = std::stod(parametervalue);
	}

	parameterstr = "beta =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.beta = std::stod(parametervalue);
	}

	parameterstr = "fw =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.fw = std::stod(parametervalue);
	}

	parameterstr = "fwsand =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.fwsand = std::stod(parametervalue);
	}

	parameterstr = "fwreef =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.fwreef = std::stod(parametervalue);
	}

	///////////////////////////////////////////////////////
	// Sediment parameters
	//
	parameterstr = "D50 =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.D50 = std::stod(parametervalue);
	}

	parameterstr = "D90 =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.D90 = std::stod(parametervalue);
	}

	parameterstr = "rhosed =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.rhosed = std::stod(parametervalue);
	}

	parameterstr = "wws =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wws = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "drydzmax =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.drydzmax = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "wetdzmax =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wetdzmax = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "maxslpchg =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.maxslpchg = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "por =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.por = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "morfac =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.morfac = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "sus =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.sus = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "bed =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.bed = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}


	parameterstr = "facsk =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.facsk = std::stod(parametervalue);
		// Value should be calculated in the sanity check
	}

	parameterstr = "facas =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.facas = std::stod(parametervalue);
		
	}

	///////////////////////////////////////////////////////
	// Timekeeping parameters
	//
	parameterstr = "dt =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.dt = std::stod(parametervalue);

	}

	parameterstr = "CFL =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.CFL = std::stod(parametervalue);

	}

	parameterstr = "sedstart =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.sedstart = std::stod(parametervalue);

	}

	parameterstr = "outputtimestep =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.outputtimestep = std::stod(parametervalue);

	}

	parameterstr = "endtime =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.endtime = std::stod(parametervalue);

	}

	///////////////////////////////////////////////////////
	// Input and output files
	//
	parameterstr = "SedThkfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.SedThkfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "wavebndfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wavebndfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "wavebndtype =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wavebndtype = std::stoi(parametervalue);
	}


	return param;
}

XBGPUParam checkparamsanity(XBGPUParam param)
{
	//First check that a bathy file was specified
	if (param.Bathymetryfile.empty())
	{
		std::cerr << "No bathymetry file specified. Please specify using 'bathy = Filename.bot'" << std::endl;
		exit(1);
	}

	return param;
}

std::string findparameter(std::string parameterstr, std::string line)
{
	std::size_t found, Numberstart, Numberend;
	std::string parameternumber;
	found = line.find(parameterstr);
	if (found != std::string::npos) // found a line that has Lonmin
	{
		//std::cout <<"found LonMin at : "<< found << std::endl;
		Numberstart = found + parameterstr.length();
		found = line.find(";");
		if (found != std::string::npos) // found a line that has Lonmin
		{
			Numberend = found;
		}
		else
		{
			Numberend = line.length();
		}
		parameternumber = line.substr(Numberstart, Numberend - Numberstart);
		//std::cout << parameternumber << std::endl;

	}
	return trim(parameternumber, " ");
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

std::string trim(const std::string& str, const std::string& whitespace)
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return ""; // no content

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
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

