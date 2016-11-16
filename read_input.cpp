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


extern "C" void readXbbndhead(const char * wavebnd, DECNUM &thetamin, DECNUM &thetamax, DECNUM &dtheta, DECNUM &dtwavbnd, int &nwavbnd, int &nwavfile)
{
	FILE * fwav;
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd, &nwavfile);
	fclose(fwav);
}
extern "C" void readbndhead(const char * wavebnd, DECNUM &thetamin, DECNUM &thetamax, DECNUM &dtheta, DECNUM &dtwavbnd, int &nwavbnd)
{
	FILE * fwav;
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd);
	printf("rtwavbnd=%f\tnwavbnd=%d\n", dtwavbnd, nwavbnd);

	fclose(fwav);
}

extern "C" void readXbbndstep(int nx, int ny, int ntheta, const char * wavebnd, int step, DECNUM &Trep, double *&qfile, double *&Stfile)
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

extern "C" void readStatbnd(int nx, int ny, int ntheta, DECNUM rho, DECNUM g, const char * wavebnd, double *&Tpfile, double *&Stfile)
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

std::vector<SLBnd> readWLfile(std::string WLfilename)
{
	std::vector<SLBnd> slbnd;

	std::ifstream fs(WLfilename);

	if (fs.fail()){
		std::cerr << WLfilename << " Water level bnd file could not be opened" << std::endl;
		exit(1);
	}

	std::string line;
	std::vector<std::string> lineelements;
	SLBnd slbndline;
	while (std::getline(fs, line))
	{
		//std::cout << line << std::endl;

		// skip empty lines
		if (!line.empty())
		{
			//Data should be in teh format :
			//BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata

			//by default we expect tab delimitation
			lineelements = split(line, '\t');
			slbndline.time = std::stod(lineelements[0]);
			slbndline.wlev = std::stod(lineelements[1]);
			
			//slbndline = readBSHline(line);
			slbnd.push_back(slbndline);
			//std::cout << line << std::endl;
		}

	}
	fs.close();

	//std::cout << slbnd[0].wlev << std::endl;


	return slbnd;
}

std::vector<WindBnd> readWNDfile(std::string WNDfilename, double grdalpha)
{
	std::vector<WindBnd> windbnd;

	std::ifstream fs(WNDfilename);

	if (fs.fail()){
		std::cerr << WNDfilename << " Wind bnd file could not be opened" << std::endl;
		exit(1);
	}

	std::string line;
	std::vector<std::string> lineelements;
	WindBnd wndbndline;
	while (std::getline(fs, line))
	{
		//std::cout << line << std::endl;

		// skip empty lines
		if (!line.empty())
		{
			//Data should be in teh format :
			//BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata

			//by default we expect tab delimitation
			lineelements = split(line, '\t');
			wndbndline.time = std::stod(lineelements[0]);
			wndbndline.spd = std::stod(lineelements[1]);
			wndbndline.dir = std::stod(lineelements[2]);
			
			//MAKE IT RELATIVE TO THE GRID X axis
			// warning this imnplies that grdalpha is in rad
			wndbndline.theta = (1.5*pi - grdalpha) - wndbndline.dir*pi / 180;

			wndbndline.U = wndbndline.spd*cos(wndbndline.theta);
			wndbndline.V = wndbndline.spd*sin(wndbndline.theta);

			
			windbnd.push_back(wndbndline);
			//std::cout << line << std::endl;
		}

	}
	fs.close();

	//std::cout << slbnd[0].wlev << std::endl;


	return windbnd;
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

	parameterstr = "fc =";
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
	
	parameterstr = "outfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.outfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

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

	parameterstr = "slbndfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.slbnd = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "windbndfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.windfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}
	

	parameterstr = "wavebndtype =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wavebndtype = std::stoi(parametervalue);
	}

	//Other parameters
	parameterstr = "GPUDEVICE =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.GPUDEVICE= std::stoi(parametervalue);
	}

	parameterstr = "nx =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.nx = std::stoi(parametervalue);
	}

	parameterstr = "ny =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ny = std::stoi(parametervalue);
	}

	parameterstr = "dx =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.dx = std::stod(parametervalue);
	}

	parameterstr = "grdalpha =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.grdalpha = std::stod(parametervalue);
	}

	parameterstr = "g =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.dx = std::stod(parametervalue);
	}

	parameterstr = "rho =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.dx = std::stod(parametervalue);
	}

	return param;
}

XBGPUParam checkparamsanity(XBGPUParam XParam, std::vector<SLBnd> slbnd, std::vector<WindBnd> wndbnd)
{
	XBGPUParam DefaultParams;

	double tiny = 0.00001;

	// Check endtime
	if (abs(XParam.endtime - DefaultParams.endtime) <= tiny)
	{
		//endtimne =0.0
		XParam.endtime = min(slbnd[slbnd.size()].time, wndbnd[wndbnd.size()].time);
	}




	return XParam;
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

double interptime(double next, double prev, double timenext, double time)
{
	return prev + (time) / (timenext)*(next - prev);
}

extern "C" void readbathyHead(std::string filename, int &nx, int &ny, double &dx, double &grdalpha )
{
	//read input data:
	//printf("bathy: %s\n", filename);
	FILE *fid;
	//int nx, ny;
	//double dx, grdalpha;
	//read md file
	fid = fopen(filename.c_str(), "r");
	fscanf(fid, "%u\t%u\t%lf\t%*f\t%lf", &nx, &ny, &dx, &grdalpha);
	grdalpha = grdalpha*pi / 180; // grid rotation
	fclose(fid);
}




extern "C" void readbathy(std::string filename, float *&zb)
{
	//read input data:
	//printf("bathy: %s\n", filename);
	FILE *fid;
	int nx, ny;
	double dx, grdalpha;
	//read md file
	fid = fopen(filename.c_str(), "r");
	fscanf(fid, "%u\t%u\t%lf\t%*f\t%lf", &nx, &ny, &dx, &grdalpha);
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


