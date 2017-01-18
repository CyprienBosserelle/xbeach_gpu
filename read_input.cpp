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


XBGPUParam readXbbndhead(XBGPUParam Param)
{
	std::ifstream fs(Param.wavebndfile);

	if (fs.fail()){
		std::cerr << Param.wavebndfile << " XBeach Reuse input bnd file could not be opened" << std::endl;
		write_text_to_log_file("ERROR: XBeach Reuse input bnd file could not be opened ");
		exit(1);
	}

	std::string line;
	std::vector<std::string> lineelements;
	
	std::getline(fs, line);
	//by default we expect tab delimitation
	lineelements = split(line, '\t');
	if (lineelements.size() < 5) // Expecting 5 parameters
	{
		lineelements.clear();
		lineelements = split(line, ' ');
	}
	
	Param.thetamin = std::stod(lineelements[0]);
	Param.thetamax = std::stod(lineelements[1]);
	Param.dtheta = std::stod(lineelements[2]);
	Param.dtbc = std::stod(lineelements[3]);
	Param.rtlength = std::stod(lineelements[4]);

	
	//FILE * fwav;
	//fwav = fopen(wavebnd, "r");
	//fscanf(fwav, "%f\t%f\t%f\t%f\t%d\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd, &nwavfile);
	//fclose(fwav);

	fs.close();

	return Param;

}
std::vector<Wavebndparam> readXbbndfile(XBGPUParam Param)
{
	std::vector<Wavebndparam> wavebndvec;
	std::ifstream fs(Param.wavebndfile);

	if (fs.fail()){
		std::cerr << Param.wavebndfile << " XBeach Reuse input bnd file could not be opened" << std::endl;
		write_text_to_log_file("ERROR: XBeach Reuse input bnd file could not be opened ");
		exit(1);
	}

	std::string line;
	std::vector<std::string> lineelements;
	Wavebndparam wavebndline;
	int linenumber = 0;
	while (std::getline(fs, line))
	{
		if (linenumber > 0)
		{
			//std::cout << line << std::endl;

			// skip empty lines
			if (!line.empty())
			{

				//by default we expect tab delimitation
				lineelements = split(line, '\t');
				if (lineelements.size() < 4)
				{
					lineelements.clear();
					lineelements = split(line, ' ');
				}
				wavebndline.time = std::stod(lineelements[0]);
				wavebndline.Trep = std::stod(lineelements[1]);
				wavebndline.qfile = lineelements[2];
				wavebndline.Efile = lineelements[3];

				//slbndline = readBSHline(line);
				wavebndvec.push_back(wavebndline);

			}
		}
		linenumber++;
	}
	return wavebndvec;
}

extern "C" void readbndhead(const char * wavebnd, DECNUM &thetamin, DECNUM &thetamax, DECNUM &dtheta, DECNUM &dtwavbnd, int &nwavbnd)
{
	FILE * fwav;
	fwav = fopen(wavebnd, "r");
	fscanf(fwav, "%f\t%f\t%f\t%f\t%d", &thetamin, &thetamax, &dtheta, &dtwavbnd, &nwavbnd);
	//printf("rtwavbnd=%f\tnwavbnd=%d\n", dtwavbnd, nwavbnd);

	fclose(fwav);
}

extern "C" void readXbbndstep(XBGPUParam Param, std::vector<Wavebndparam> wavebnd,int step, DECNUM &Trep, double *&qfile, double *&Stfile)
{

	int nx, ny, ntheta;
	nx = Param.nx;
	ny = Param.ny;
	ntheta = Param.ntheta;
	FILE * fXq, *fXE;
	
	double dummy;
	size_t result;
	DECNUM thetamin, thetamax, dtheta, dtwavbnd;
	int nwavbnd, nwavfile;
	nwavbnd = ceil(Param.rtlength / Param.dtbc); 


	printf("Reading next wave bnd file... ");
	write_text_to_log_file("Reading next bnd file... ");

	
	//printf("Xbq: %s\n",Xbqfile);
	//printf("Xbe: %s\n",XbEfile);

	fXq = fopen(wavebnd[step].qfile.c_str(), "rb");
	if (!fXq)
	{
		printf("Unable to open file %s\t", wavebnd[step].qfile.c_str());
		write_text_to_log_file("Unable to open file :" + wavebnd[step].qfile);
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


	fXE = fopen(wavebnd[step].Efile.c_str(), "rb");
	for (int nn = 0; nn < ntheta*ny*nwavbnd; nn++)
	{
		fread(&dummy, sizeof(double), 1, fXE);
		//printf("St=%f\n ",dummy);
		//Stfile[nn] = dummy;
		Stfile[nn] = (DECNUM)dummy;
	}
	fclose(fXE);

	Trep = wavebnd[step].Trep;

	printf("done \n");
	write_text_to_log_file("done");

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
	write_text_to_log_file("Reading bnd file... ");

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
	write_text_to_log_file("done");

}

std::vector<SLBnd> readWLfile(std::string WLfilename)
{
	std::vector<SLBnd> slbnd;

	std::ifstream fs(WLfilename);

	if (fs.fail()){
		std::cerr << WLfilename << " Water level bnd file could not be opened" << std::endl;
		write_text_to_log_file("ERROR: Water level bnd file could not be opened ");
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
			if (lineelements.size() < 2)
			{
				lineelements.clear();
				lineelements = split(line, ' ');
			}
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
		write_text_to_log_file("ERROR: Wind bnd file could not be opened ");
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



std::vector<Wavebndparam> ReadCstBnd(XBGPUParam XParam)
{
	std::vector<Wavebndparam> wavebnd;


	std::ifstream fs(XParam.wavebndfile);

	if (fs.fail()){
		std::cerr << XParam.wavebndfile << " Wave bnd file could not be opened" << std::endl;
		write_text_to_log_file("ERROR: Wave bnd file could not be opened ");
		exit(1);
	}

	std::string line;
	std::vector<std::string> lineelements;
	Wavebndparam waveline;

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
			if (lineelements.size() < 5) // If we cant find all the elements it must be space delimited
			{
				lineelements.clear();
				lineelements = split(line, ' ');
			}

			waveline.time = std::stod(lineelements[0]);
			waveline.Hs = std::stod(lineelements[1]);
			waveline.Tp = std::stod(lineelements[2]);
			// make bnd normal wave direction
			waveline.Dp = (1.5*pi - XParam.grdalpha) - std::stod(lineelements[3])*pi / 180; // Why make it in degree?
			waveline.s = std::stod(lineelements[4]);
			wavebnd.push_back(waveline);
		}
	}
	fs.close();

	return wavebnd;
}
std::vector<Wavebndparam> ReadJSWPBnd(XBGPUParam XParam)
{
	std::vector<Wavebndparam> wavebnd;


	std::ifstream fs(XParam.wavebndfile);

	if (fs.fail()){
		std::cerr << XParam.wavebndfile << " Wave bnd file could not be opened" << std::endl;
		write_text_to_log_file("ERROR: Wave bnd file could not be opened ");
		exit(1);
	}

	std::string line;
	std::vector<std::string> lineelements;
	Wavebndparam waveline;

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
			if (lineelements.size() < 6) // If we cant find all the elements it must be space delimited
			{
				lineelements.clear();
				lineelements = split(line, ' ');
			}

			waveline.time = std::stod(lineelements[0]);
			waveline.Hs = std::stod(lineelements[1]);
			waveline.Tp = std::stod(lineelements[2]);
			// make bnd normal wave direction
			waveline.Dp = (1.5*pi - XParam.grdalpha) - std::stod(lineelements[3])*pi / 180; // Why make it in degree?
			waveline.s = std::stod(lineelements[4]);
			waveline.gamma = std::stod(lineelements[5]);
			wavebnd.push_back(waveline);
		}
	}
	fs.close();

	return wavebnd;
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

	parameterstr = "roller =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.roller = std::stoi(parametervalue);
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

	

	// Below is a bit more complex than usual because more than 1 node can be outputed as a timeseries
	parameterstr = "TSOfile =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.TSoutfile.push_back(parametervalue);
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "TSnode =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		std::vector<std::string> nodes = split(parametervalue, ',');
		TSnode node;
		node.i = std::stoi(nodes[0]);
		node.j = std::stoi(nodes[1]);
		param.TSnodesout.push_back(node);

		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	//outvars
	parameterstr = "outvars =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		std::vector<std::string> vars = split(parametervalue, ',');
		for (int nv = 0; nv < vars.size(); nv++)
		{
			//Verify that the variable name makes sense?
			//Need to add more here
			std::vector<std::string> SupportedVarNames = { "hh", "uu", "vv", "wci", "zs", "zb", "cfm", "dzb", "stdep", "Fx", "Fy", "cgx", "cgy", "cx", "cy", "ctheta", "D", "E", "H", "urms", "ueu", "vev", "thetamean", "Hmean", "uumean", "vvmean", "hhmean", "zsmean", "Cmean", "sigm", "k", "c", "kh", "cg", "sinh2kh", "dhdx", "dhdy", "dudx", "dudy", "dvdx", "dvdy", "C", "R", "DR", "ee", "vmageu", "vmagev", "dzsdx", "dzsdy", "dzsdt", "fwm", "hu", "hum", "hv", "hvm", "uv", "vu", "ududx", "vdvdy", "udvdx", "vdudy", "ust", "rr", "kturb", "rolthick", "ceqsg" };

			for (int isup = 0; isup < SupportedVarNames.size(); isup++)
			{
				std::string vvar = trim(vars[nv]," ");
				//std::cout << "..." << vvar << "..." << std::endl;
				if (vvar.compare(SupportedVarNames[isup]) == 0)
				{
					param.outvars.push_back(vvar);
					break;
				}
			}
			
		}
		

		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "wavebndtype =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.wavebndtype = std::stoi(parametervalue);
	}

	parameterstr = "thetamin =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.thetamin = std::stod(parametervalue);
	}

	parameterstr = "thetamax =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.thetamax = std::stod(parametervalue);
	}

	parameterstr = "dtheta =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.dtheta = std::stod(parametervalue);
	}

	parameterstr = "dtbc =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.dtbc = std::stod(parametervalue);
	}

	parameterstr = "rtlength =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.rtlength = std::stod(parametervalue);
	}

	parameterstr = "sprdthr =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.sprdthr = std::stod(parametervalue);
	}

	parameterstr = "random =";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.random = std::stoi(parametervalue);
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

	double tiny = 0.0000001;

	//Check Bathy input type
	std::string bathyext;
	std::vector<std::string> extvec = split(XParam.Bathymetryfile, '.');
	bathyext = extvec.back();
	if (bathyext.compare("nc") == 0)
	{
		if (abs(XParam.grdalpha - DefaultParams.grdalpha) < tiny)
		{
			
			write_text_to_log_file("For nc of bathy file please specify grdalpha in the XBG_param.txt (if different then default [0])");
		}
	}
	if (bathyext.compare("dep") == 0 || bathyext.compare("bot") == 0)
	{
		if (XParam.nx <= 0 || XParam.ny <= 0 || XParam.dx < tiny)
		{
			std::cerr << "FATAL ERROR: nx or ny or dx were not specified. These parameters are required when using ." << bathyext << " file" << std::endl;
			write_text_to_log_file("FATAL ERROR: nx or ny or dx were not specified. These parameters are required when using ." + bathyext + " file");
			exit(1);
		}
	}

	if (XParam.nx <= 0 || XParam.ny <= 0 || XParam.dx < tiny)
	{
		std::cerr << "FATAL ERROR: nx or ny or dx could not specified." << std::endl;
		write_text_to_log_file("FATAL ERROR: nx or ny or dx could not specified.");
		exit(1);
	}


	// Check whether endtime was specified by the user
	if (abs(XParam.endtime - DefaultParams.endtime) <= tiny)
	{
		//No; i.e. endtimne =0.0
		if (slbnd.back().time>0.0 && wndbnd.back().time > 0.0)
		{
			XParam.endtime = min(slbnd.back().time, wndbnd.back().time);
		}

		
	}
	else
	{
		//Check that endtime is no longer than the shortest boundary
		double endbnd = min(slbnd.back().time, wndbnd.back().time);

		XParam.endtime = min(XParam.endtime, endbnd);
	}


	// Check that a wave bnd file was specified otherwise kill the app
	// This is temporary until the wave boundary scheme is improved
	if (XParam.wavebndfile.empty())
	{
		std::cerr << "Fatal error: No wave boundary file specified. Please specify using 'wavebndfile = wave_boundary.bnd;'" << std::endl;
		write_text_to_log_file("ERROR: No wave boundary file specified. Please specify using 'wavebndfile = wave_boundary.bnd;");
		exit(1);
	}

	// Check that outputtimestep is not zero, so at least the first and final time step are saved
	// If only the model stepup is needed than just run with endtime=0.0
	if (abs(XParam.outputtimestep - DefaultParams.outputtimestep) <= tiny)
	{
		XParam.outputtimestep = XParam.endtime;
		//otherwise there is really no point running the model
	}

	//Check that sand and reef friction is not zero (Default) if it is the case then use the default cf
	if (abs(XParam.cfsand - DefaultParams.cfsand) < tiny)
	{
		if (XParam.cf < tiny)
		{
			XParam.cf = tiny;
		}
		XParam.cfsand = XParam.cf;
	}

	if (abs(XParam.cfreef - DefaultParams.cfreef) < tiny)
	{
		
		XParam.cfreef = XParam.cfsand;
	}

	//Check that sand and reef friction is not zero (Default) if it is the case then use the default cf
	if (abs(XParam.fwsand - DefaultParams.fwsand) < tiny)
	{
		if (XParam.fw < tiny)
		{
			XParam.fw = tiny;
		}
		XParam.fwsand = XParam.fw;
	}

	if (abs(XParam.fwreef - DefaultParams.fwreef) < tiny)
	{

		XParam.fwreef = XParam.fwsand;
	}


	// Check that if smagorinsky formulation is used then nuh == samgo otherwise use the specified value for smago
	if (XParam.usesmago == 1)
	{
		//print a warning message if nuh was user specified
		if (abs(XParam.nuh - DefaultParams.nuh) > tiny)
		{
			std::cout << "WARNING: Using Smagorinsky formulation. User specified value for 'nuh' will be ignored. Use 'smag = 0.3' to control the smagorinsky parameter" << std::endl;
			write_text_to_log_file("WARNING: Using Smagorinsky formulation. User specified value for 'nuh' will be ignored. Use 'smag = 0.3' to control the smagorinsky parameter");
		}

		//force nuh to be == to smago
		XParam.nuh = XParam.smag;
	}

	//Check that there are as many file specified for Time series output as there are vectors of nodes
	if (XParam.TSoutfile.size() != XParam.TSnodesout.size())
	{
		// Issue a Warning
		std::cout << "WARNING: the number of timeseries output files is not equal to the number of nodes specified" << std::endl;
		std::cout << "for each location where timeseries output file is required, the XBG_param.txt file shoud contain 2 lines see example felow to extract in 2 locations:" << std::endl;
		std::cout << "TSOfile = Reef_Crest.txt" << std::endl;
		std::cout << "TSnode = 124,239;" << std::endl;
		std::cout << "TSOfile = shore.txt" << std::endl;
		std::cout << "TSnode = 233,256;" << std::endl;

		write_text_to_log_file("WARNING: the number of timeseries output files is not equal to the number of nodes specified");
		write_text_to_log_file("for each location where timeseries output file is required, the XBG_param.txt file shoud contain 2 lines see example felow to extract in 2 locations:");
		write_text_to_log_file("TSOfile = Reef_Crest.txt");
		write_text_to_log_file("TSnode = 124,239;");
		write_text_to_log_file("TSOfile = Shore.txt");
		write_text_to_log_file("TSnode = 233,256;");
		//min not defined for const so use this convoluted statement below
		int minsize;
		if (XParam.TSoutfile.size() > XParam.TSnodesout.size())
		{
			minsize = XParam.TSnodesout.size();
		}

		if (XParam.TSoutfile.size() < XParam.TSnodesout.size())
		{
			minsize = XParam.TSoutfile.size();
		}

				
		XParam.TSoutfile.resize(minsize);
		XParam.TSnodesout.resize(minsize);
	}

	//Chaeck that if timeseries output nodes are specified that they are within nx and ny
	if (XParam.TSnodesout.size() > 0)
	{
		for (int o = 0; o < XParam.TSnodesout.size(); o++)
		{
			if (XParam.TSnodesout[o].i < 0)
			{
				//An idiot is in charge
				XParam.TSnodesout[o].i = 0;
			}

			if (XParam.TSnodesout[o].i > XParam.nx-1)
			{
				XParam.TSnodesout[o].i = XParam.nx - 1;
			}

			if (XParam.TSnodesout[o].j < 0)
			{
				//An idiot is in charge
				XParam.TSnodesout[o].j = 0;
			}

			if (XParam.TSnodesout[o].j > XParam.ny - 1)
			{
				XParam.TSnodesout[o].j = XParam.ny - 1;
			}
		}

	}

	if (XParam.outvars.empty() && XParam.outputtimestep > 0)
	{
		//a nc file was specified but no output variable were specified
		std::vector<std::string> SupportedVarNames = { "zb", "zs", "uu", "vv", "H", "thetamean", "D", "urms", "ueu", "vev", "C", "dzb", "Fx", "Fy", "hh", "Hmean", "uumean", "vvmean", "hhmean", "zsmean", "Cmean" };
		for (int isup = 0; isup < SupportedVarNames.size(); isup++)
		{
			XParam.outvars.push_back(SupportedVarNames[isup]);
				
		}

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
		if (!item.empty())//skip empty tokens
		{
			elems.push_back(item);
		}
		
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
	//grdalpha = grdalpha*pi / 180; // grid rotation
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

extern "C" void readXBbathy(std::string filename, int nx,int ny, float *&zb)
{
	//read input data:
	//printf("bathy: %s\n", filename);
	FILE *fid;
	
	//read md file
	fid = fopen(filename.c_str(), "r");
	
	

	
	//int jreadzs;
	for (int jnod = 0; jnod < ny; jnod++)
	{

		

		for (int inod = 0; inod < nx; inod++)
		{
			fscanf(fid, "%f", &zb[inod + (jnod)*nx]);

		}
	}

	fclose(fid);
}


void write_text_to_log_file(std::string text)
{
	std::ofstream log_file(
		"XBG_log.txt", std::ios_base::out | std::ios_base::app);
	log_file << text << std::endl;
	log_file.close(); //destructoir implicitly does it
}

void SaveParamtolog(XBGPUParam XParam)
{
	write_text_to_log_file("#################################");
	write_text_to_log_file("# Bathymetry file");
	write_text_to_log_file("bathy = " + XParam.Bathymetryfile + ";");
	write_text_to_log_file("nx = " + std::to_string(XParam.nx) + ";");
	write_text_to_log_file("ny = " + std::to_string(XParam.ny) + ";");
	write_text_to_log_file("dx = " + std::to_string(XParam.dx) + ";");
	write_text_to_log_file("grdalpha = " + std::to_string(XParam.grdalpha*180.0/pi) + ";");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Model controls");
	write_text_to_log_file("swave = " + std::to_string(XParam.swave) + ";");
	write_text_to_log_file("flow = " + std::to_string(XParam.flow) + ";");
	write_text_to_log_file("sedtrans = " + std::to_string(XParam.sedtrans) + ";");
	write_text_to_log_file("morphology = " + std::to_string(XParam.morphology) + ";");
	write_text_to_log_file("gpudevice = " + std::to_string(XParam.GPUDEVICE) + ";");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Flow parameters");
	write_text_to_log_file("eps = " + std::to_string(XParam.eps) + ";");
	write_text_to_log_file("cfsand = " + std::to_string(XParam.cfsand) + ";");
	write_text_to_log_file("cfreef = " + std::to_string(XParam.cfreef) + ";");
	write_text_to_log_file("usesmago = " + std::to_string(XParam.usesmago) + ";");
	if (XParam.usesmago == 0)
	{
		write_text_to_log_file("nuh = " + std::to_string(XParam.sedtrans) + ";");
	}
	else
	{
		write_text_to_log_file("smag = " + std::to_string(XParam.smag) + ";");
	}
	write_text_to_log_file("nuhfac = " + std::to_string(XParam.nuhfac) + ";");
	write_text_to_log_file("lat = " + std::to_string(XParam.lat) + ";");
	write_text_to_log_file("fc = " + std::to_string(XParam.fc) + ";");
	write_text_to_log_file("Cd = " + std::to_string(XParam.Cd) + ";");
	write_text_to_log_file("wci = " + std::to_string(XParam.wci) + ";");
	write_text_to_log_file("hwci = " + std::to_string(XParam.hwci) + ";");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Waves parameters");
	write_text_to_log_file("breakmodel = " + std::to_string(XParam.breakmodel) + ";");
	write_text_to_log_file("gamma = " + std::to_string(XParam.gammaa) + ";");
	write_text_to_log_file("n = " + std::to_string(XParam.n) + ";");
	write_text_to_log_file("alpha = " + std::to_string(XParam.alpha) + ";");
	write_text_to_log_file("gammax = " + std::to_string(XParam.gammax) + ";");
	write_text_to_log_file("beta = " + std::to_string(XParam.beta) + ";");
	write_text_to_log_file("fwsand = " + std::to_string(XParam.fwsand) + ";");
	write_text_to_log_file("fwreef = " + std::to_string(XParam.fwreef) + ";");
	write_text_to_log_file("roller = " + std::to_string(XParam.roller) + ";");
	write_text_to_log_file("thetamin = " + std::to_string(XParam.thetamin) + ";");
	write_text_to_log_file("thetamax = " + std::to_string(XParam.thetamax) + ";");
	write_text_to_log_file("dtheta = " + std::to_string(XParam.dtheta) + "; ");
	write_text_to_log_file("ntheta = " + std::to_string(XParam.ntheta) + "; ");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Wave boundary parameters");
	write_text_to_log_file("dtbc = " + std::to_string(XParam.dtbc) + "; ");
	write_text_to_log_file("rtlength = " + std::to_string(XParam.rtlength) + "; ");
	write_text_to_log_file("sprdthr = " + std::to_string(XParam.sprdthr) + "; ");
	write_text_to_log_file("random = " + std::to_string(XParam.random) + "; ");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Sediment parameters");
	write_text_to_log_file("D50 = " + std::to_string(XParam.D50) + ";");
	write_text_to_log_file("D90 = " + std::to_string(XParam.D90) + ";");
	write_text_to_log_file("rhosed = " + std::to_string(XParam.rhosed) + ";");
	write_text_to_log_file("wws = " + std::to_string(XParam.wws) + ";");
	write_text_to_log_file("drydzmax = " + std::to_string(XParam.drydzmax) + ";");
	write_text_to_log_file("wetdzmax = " + std::to_string(XParam.wetdzmax) + ";");
	write_text_to_log_file("maxslpchg = " + std::to_string(XParam.maxslpchg) + ";");
	write_text_to_log_file("por = " + std::to_string(XParam.por) + ";");
	write_text_to_log_file("morfac = " + std::to_string(XParam.morfac) + ";");
	write_text_to_log_file("sus = " + std::to_string(XParam.sus) + ";");
	write_text_to_log_file("bed = " + std::to_string(XParam.bed) + ";");
	write_text_to_log_file("facsk = " + std::to_string(XParam.facsk) + ";");
	write_text_to_log_file("facas = " + std::to_string(XParam.facas) + ";");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Timekeeping parameters");
	write_text_to_log_file("CFL = " + std::to_string(XParam.CFL) + ";");
	write_text_to_log_file("sedstart = " + std::to_string(XParam.sedstart) + ";");
	write_text_to_log_file("outputtimestep = " + std::to_string(XParam.outputtimestep) + ";");
	std::string alloutvars= "";
	for (int nvar = 0; nvar < XParam.outvars.size(); nvar++)
	{
		if (nvar > 0)
		{
			alloutvars = alloutvars + ", ";
		}
		alloutvars = alloutvars + XParam.outvars[nvar];
	}
	write_text_to_log_file("outvars = " + alloutvars + ";");


	write_text_to_log_file("endtime = " + std::to_string(XParam.endtime) + ";");
	write_text_to_log_file("\n");
	write_text_to_log_file("# Files");
	write_text_to_log_file("outfile = " + XParam.outfile + ";");
	write_text_to_log_file("SedThkfile = " + XParam.SedThkfile + ";");
	write_text_to_log_file("wavebndfile = " + XParam.wavebndfile + ";");
	write_text_to_log_file("wavebndtype = " + std::to_string(XParam.wavebndtype) + ";");
	write_text_to_log_file("slbndfile = " + XParam.slbnd + ";");
	write_text_to_log_file("windbndfile = " + XParam.windfile + ";");
	if (!XParam.TSoutfile.empty())
	{
		for (int o = 0; o < XParam.TSoutfile.size(); o++)
		{
			write_text_to_log_file("TSOfile = " + XParam.TSoutfile[o] + ";");
			write_text_to_log_file("TSnode = " + std::to_string(XParam.TSnodesout[o].i) + "," + std::to_string(XParam.TSnodesout[o].j) + ";");
		}
	}
	write_text_to_log_file("\n");
	write_text_to_log_file("# Others");
	write_text_to_log_file("g = " + std::to_string(XParam.g) + ";");
	write_text_to_log_file("rho = " + std::to_string(XParam.rho) + ";");
}



