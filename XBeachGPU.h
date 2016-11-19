#ifndef XBeachGPU_H
#define XBeachGPU_H
 
#define pi 3.14159265

#include <stdio.h>
#include <math.h>


#include <iostream>
#include <fstream>
#include <iomanip>


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cmath>
#include <sstream>
#include <iterator>
#include <netcdf.h>
#include <algorithm>
#include <vector>
#include <ctime>

//using namespace std;



using DECNUM = float;


//class template to create vecotr ofg output locations fro the model
class TSnode{
public:
	int i, j;
};


class XBGPUParam{
public:
	
	//General parameters 
	//int modeltype;// Type of model: 1: wave only; 2: currents only 3: waves+currents 4:waves+currents+sediment(+ morphology if morfac>0) //Obsolete Need to remove
	int swave=1, flow=1, sedtrans=0, morphology=0;
	int GPUDEVICE=0;// What GPU device to use default is the firt one availabe (i.e. 0) CPU only should be -1 (not operational yet) and 1 for second GPU etc...
	int nx, ny; // grid size
	double dx; // grid resolution
	double grdalpha=0.0; // grid rotation Y axis from the North input in degrees but later converted to rad
	
	//Flow parameters
	double g = 9.81;
	double rho = 1025.0;
	double eps=0.01;//drying height in m
	double cf = 0.01; // bottom friction for flow model cf 
	double cfsand = 0.0, cfreef=0.0;// bottom friction for sand and for reef area (Reef and sand discrimination is done based on sediment thickness file if none is present cf2 cannot be used )
	double nuh = 1.0;// Viscosity coeff ,
	double nuhfac = 1.0;//nuhfac=1.0f;//0.001f; //viscosity coefficient for roller induced turbulent horizontal viscosity// it should be small contrary to what XBeach recommend as default
	int usesmago = 1;// Uses smagorynsky formulation to calculate viscosity 0: No 1: Yes
	double smag = 0.3; // Smagorinsky coeff only used if usesmago = 1
	double lat = 0.0; // Latitude of the grid use negative for south hemisphere (this implies the grid is small on earth scale)
	double fc = 0.0;
	double Cd = 0.002; // Wind drag coeff
	double wci = 0;// Wave current interaction switch (can also be used as a number between 0 and 1 to reduce the interaction if unstable) 
	double hwci = 0.1; // hwci=0.010f;//min depth for wci

	//Wave parameters
	int breakmodel = 1;// Wave dissipation model 1: roelvink 2: Baldock. use 1 for unsteady runs (i.e. with wave group) and use 2 for steady runs
	double gammaa = 0.6; // Wave breaking gamma param 
	double n = 8.0; // exponential; in Roelving breaking model
	double alpha = 1.0; // calibration for wave dissipation (should be 1)
	double gammax = 2.0; //gammax=2.0f; //maximum ratio Hrms/hh
	double beta = 0.15; // Roller slope dissipation param
	double fw = 0.001;//Wave bottom dissipation parameters fw 
	double fwsand = 0.0, fwreef = 0.0; //Wave bottom dissipation parameters fw is for sand fw2 is for reefs.see cf comments
	
	//Sediment parameters
	double D50=0.00038, D90=0.00053; // sand grain size in m
	double rhosed = 2650.0; // sand density in kg/m3
	double wws=0.0509; //// sand fall velocity (should be calculated) m/s
	double drydzmax=1.0, wetdzmax=2.0; // max slope in avalannching model
	double maxslpchg=0.01; // max change within a step to avoid avalanching tsunami
	double por=0.4; // sand porosity (should not be constant)
	double morfac=1.0; // morphological factor 0 no changes in morphology 1 normal changes in morpho >1 accelerated morphological changes (beware this doesn't accelerate the bnd you have to do this manually)
	double sus=1.0, bed=1.0; // calibration coeff for suspended load and bed load
	double facsk=0.2, facas=0.2;// calibration factor for wave skewness and Asymetry

	// File
	std::string Bathymetryfile;// bathymetry file name
	std::string SedThkfile; // Structure file 
	std::string wavebndfile;// wave bnd file
	int wavebndtype = 2; // 1 is quasistationary wave spectrum; 2 is for infrgravity and long bound waves Xbeach type
	std::string slbnd; // tide/surge bnd file
	std::string windfile; // Wind forcing file
	std::string outfile = "XBGPU_output.nc"; //outputfile
	

	//Timekeeping
	double dt=0.0, CFL=0.7;// Model time step in s. either one should be defined if both then dt is constant
	double sedstart=3600.0;// time to start sediment transport and morpho
	double outputtimestep=0.0; //number of seconds between output 0.0 for none
	double endtime=0.0; // Total runtime in s will be calculated based on bnd input as min(length of the shortest time series, user defined)
	
	//Timeseries output
	std::vector<std::string> TSoutfile; //filename of output time series (Only save time, H U,V and zs)
	std::vector<TSnode> TSnodesout; // vector containing i and j of each variables

};


class SLBnd {
public:
	double time, wlev;
};


class WindBnd{
public:
	double time; // time in s
	double U, V; // U and V wind in m/s relative to the grid
	double spd, dir; //speed in m/s and dir in deg relative to true North
	double theta; // wind direction in rad relative to the grid y axis
};




// additional functions
void makjonswap(DECNUM hm0gew,DECNUM fp,DECNUM mainang,DECNUM rt,DECNUM scoeff,DECNUM gam,DECNUM * theta,int ntheta,DECNUM& TTrep,DECNUM * &Stt);
extern "C" void creatncfile(XBGPUParam XParam, DECNUM totaltime, DECNUM *zb, DECNUM *zs, DECNUM * uu, DECNUM * vv, DECNUM * H, DECNUM * Tp, DECNUM * Dp, DECNUM * D, DECNUM * Urms, DECNUM * ueu, DECNUM * vev, DECNUM * C, DECNUM *Fx, DECNUM *Fy, DECNUM * hh, DECNUM *Hmean, DECNUM *uumean, DECNUM *vvmean, DECNUM *hhmean, DECNUM *zsmean, DECNUM *Cmean);
extern "C" void writestep2nc(XBGPUParam XParam, DECNUM totaltime, DECNUM *zb, DECNUM *zs, DECNUM * uu, DECNUM * vv, DECNUM * H, DECNUM * Tp, DECNUM * Dp, DECNUM * D, DECNUM * Urms, DECNUM *ueu, DECNUM * vev, DECNUM * C, DECNUM *dzb, DECNUM *Fx, DECNUM *Fy, DECNUM *hh, DECNUM *Hmean, DECNUM *uumean, DECNUM *vvmean, DECNUM *hhmean, DECNUM *zsmean, DECNUM *Cmean);

extern "C" void create3dnc(int nx,int ny,int nt,DECNUM dx,DECNUM totaltime,DECNUM *theta,DECNUM * var);
extern "C" void write3dvarnc(int nx,int ny,int nt,DECNUM totaltime,DECNUM * var);

extern "C" void read3Dnc(int nx, int ny,int ntheta,char ncfile[],DECNUM * &ee);
extern "C" void read2Dnc(int nx, int ny,char ncfile[],DECNUM * &hh);

extern "C" void readXbbndhead(const char * wavebndfile,DECNUM &thetamin,DECNUM &thetamax,DECNUM &dtheta,DECNUM &dtwavbnd,int &nwavbnd,int &nwavfile);
extern "C" void readXbbndstep(int nx, int ny,int ntheta,const char * wavebndfile,int step,DECNUM &Trep,double *&qfile,double *&Stfile );
extern "C" void readStatbnd(int nx, int ny,int ntheta,DECNUM rho,DECNUM g,const char * wavebndfile,double *&Tpfile,double *&Stfile );
extern "C" void readbndhead(const char * wavebndfile,DECNUM &thetamin,DECNUM &thetamax,DECNUM &dtheta,DECNUM &dtwavbnd,int &nwavbnd);

//Below is for the new CPU routine
extern "C" int mminusC(int ix,int nx);
extern "C" int pplusC(int ix, int nx);
extern "C" int mminus2C(int ix,int nx);
extern "C" int pplus2C(int ix, int nx);
extern "C" int signC(DECNUM x);
extern "C" void set_bndCPU(int nx, int ny,DECNUM Trep,int ntheta,DECNUM * theta,DECNUM *&sigm);
extern "C" void dispersion_initCPU(int nx,int ny,DECNUM twopi,DECNUM g,DECNUM aphi,DECNUM bphi,DECNUM * sigm,DECNUM * hh,DECNUM * &cg);
extern "C" void updatezomCPU(int nx, int ny,DECNUM cf,DECNUM cf2,DECNUM fw,DECNUM fw2,DECNUM * structdepth, DECNUM * &cfm,DECNUM * &fwm);
extern "C" void ubndCPU(int nx, int ny, DECNUM dx, DECNUM dt,DECNUM g, DECNUM rho,DECNUM totaltime,DECNUM wavbndtime,DECNUM rt,DECNUM zsbnd,DECNUM Trep,DECNUM * qbndold, DECNUM * qbndnew,DECNUM *&zs, DECNUM * &uu,DECNUM * &vv, DECNUM *vu, DECNUM * umean, DECNUM * vmean,DECNUM * zb,DECNUM * cg,DECNUM * hum, DECNUM * zo, DECNUM *Fx,DECNUM *&hh);
extern "C" void offshorebndWavCPU(int nx, int ny,int ntheta,DECNUM totaltime,DECNUM Trep,DECNUM *St,DECNUM *&sigm, DECNUM *&ee);
extern "C" void sanityCPU(int nx, int ny,DECNUM eps,DECNUM * hh,DECNUM * sigm, int ntheta,DECNUM * &ee);
extern "C" void dispersionCPU(int nx,int ny,DECNUM twopi,DECNUM g,DECNUM aphi,DECNUM bphi,DECNUM * sigm,DECNUM * hh,DECNUM * &k,DECNUM * &c,DECNUM * &kh,DECNUM * &sinh2kh,DECNUM * &cg);
extern "C" void calcwciCPU(int nx, int ny, DECNUM wci, DECNUM hwci, DECNUM * hh, DECNUM * &wcig);
extern "C" void slopesCPU(int nx,int ny,DECNUM dx,DECNUM * hh,DECNUM * uu,DECNUM * vv,DECNUM * &dhdx,DECNUM * &dhdy,DECNUM * &dudx,DECNUM * &dudy,DECNUM * &dvdx,DECNUM * &dvdy);
extern "C" DECNUM slopes2DxCPU(int nx,DECNUM dx,int i,int ix, int iy,DECNUM * hh);
extern "C" DECNUM slopes2DyCPU(int nx,int ny,DECNUM dx,int i,int ix, int iy,DECNUM * hh);
extern "C" void propagthetaCPU(int nx,int ny,int ntheta,DECNUM * wci,DECNUM * &ctheta,DECNUM *cxsth,DECNUM *sxnth,DECNUM *dhdx,DECNUM *dhdy,DECNUM *dudx,DECNUM *dudy,DECNUM *dvdx,DECNUM *dvdy,DECNUM *sigm,DECNUM *kh);
extern "C" void actionCPU(int ntheta,int nx,int ny,DECNUM * &ee,DECNUM * sigm);
extern "C" void xadvecupwindCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM * wci,DECNUM *ee,DECNUM *cg,DECNUM *cxsth,DECNUM *uu,DECNUM * &xadvec);
extern "C" void xadvecupwind2CPU(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *cxsth, DECNUM *uu, DECNUM * xadvec);

extern "C" void yadvecupwindCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM * wci,DECNUM *ee,DECNUM *cg,DECNUM *sxnth,DECNUM *vv,DECNUM * &yadvec);
extern "C" void yadvecupwind2CPU(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM * wci, DECNUM *ee, DECNUM *cg, DECNUM *sxnth, DECNUM *vv, DECNUM * yadvec);

extern "C" void thetaadvecuw1hoCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM wci,DECNUM *ee,DECNUM *ctheta,DECNUM * &thetaadvec);
extern "C" void thetaadvecuw2hoCPU(int nx, int ny, int ntheta, DECNUM dtheta, DECNUM dx, DECNUM dt, DECNUM wci, DECNUM *ee, DECNUM *ctheta, DECNUM * thetaadvec);

extern "C" void eulerupwindCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM dx,DECNUM dt,DECNUM wci,DECNUM *&ee,DECNUM *xadvec,DECNUM *yadvec,DECNUM * thetaadvec);
extern "C" void rollerlatbndCPU(int nx,int ny,int ntheta,DECNUM eps,DECNUM *hh,DECNUM *&rr);
extern "C" void energyCPU(int nx,int ny,int ntheta,DECNUM * &ee,DECNUM * sigm);
extern "C" void energintCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM rho,DECNUM g,DECNUM gammax,DECNUM * &E,DECNUM * &H,DECNUM * hh,DECNUM * &ee);
extern "C" void roelvinkCPU(int nx, int ny,DECNUM rho,DECNUM g,DECNUM gamma,DECNUM alpha,DECNUM n,DECNUM Trep,DECNUM * fwm,DECNUM * cfm,DECNUM *hh,DECNUM *H,DECNUM *E,DECNUM *&D, DECNUM *k);
extern "C" void baldockCPU(int nx, int ny,DECNUM rho,DECNUM g,DECNUM gamma,DECNUM alpha,DECNUM n,DECNUM Trep,DECNUM *fwm,DECNUM * cfm,DECNUM *hh,DECNUM *H,DECNUM *E,DECNUM *&D, DECNUM *k);
extern "C" void dissipationCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM eps,DECNUM dt,DECNUM g,DECNUM beta,DECNUM * wci,DECNUM *hh,DECNUM *&ee,DECNUM *D,DECNUM *E,DECNUM *&rr,DECNUM *c,DECNUM *cxsth,DECNUM *sxnth,DECNUM * uu,DECNUM * vv,DECNUM *&DR,DECNUM *&R);
extern "C" void meandirCPU(int nx,int ny,int ntheta,DECNUM rho,DECNUM g,DECNUM dtheta,DECNUM * ee,DECNUM * thet, DECNUM * &thetamean, DECNUM * &E,DECNUM * &H);
extern "C" void radstressCPU(int nx,int ny, int ntheta,DECNUM dx,DECNUM dtheta,DECNUM * ee,DECNUM *rr,DECNUM * cxsth,DECNUM * sxnth,DECNUM * cg,DECNUM * c,DECNUM * &Sxx,DECNUM * &Sxy,DECNUM * &Syy);
extern "C" void wavforceCPU(int nx,int ny, int ntheta,DECNUM dx,DECNUM dtheta,DECNUM * Sxx,DECNUM * Sxy,DECNUM * Syy,DECNUM * &Fx,DECNUM * &Fy,DECNUM * hh);
extern "C" void twodimbndnoixCPU(int nx,int ny,DECNUM eps,DECNUM * hh,DECNUM * &F);
extern "C" void twodimbndCPU(int nx,int ny,DECNUM eps,DECNUM * hh,DECNUM * &F);
extern "C" void breakerdelayCPU(int nx,int ny,int ntheta,DECNUM dtheta,DECNUM g, DECNUM rho,DECNUM Trep, DECNUM eps,DECNUM * &urms,DECNUM *&ust,DECNUM *H,DECNUM *E,DECNUM *c,DECNUM *k, DECNUM *hh,DECNUM *rr);
extern "C" void addavg_varCPU(int nx, int ny,DECNUM * &Varmean,DECNUM * Var);
extern "C" void divavg_varCPU(int nx, int ny,DECNUM ntdiv,DECNUM * &Varmean);
extern "C" void resetavg_varCPU(int nx, int ny,DECNUM * &Varmean);

XBGPUParam readparamstr(std::string line, XBGPUParam param);
std::string findparameter(std::string parameterstr, std::string line);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
std::string trim(const std::string& str, const std::string& whitespace);
XBGPUParam checkparamsanity(XBGPUParam XParam, std::vector<SLBnd> slbnd, std::vector<WindBnd> wndbnd);
std::vector<SLBnd> readWLfile(std::string WLfilename);
std::vector<WindBnd> readWNDfile(std::string WNDfilename, double grdalpha);
double interptime(double next, double prev, double timenext, double time);


extern "C" void readbathyHead(std::string filename, int &nx, int &ny, double &dx, double &grdalpha);
extern "C" void readbathy(std::string filename, float *&zb);
// General functions
//void CUDA_CHECK(cudaError CUDerr);

template <class T> const T& min (const T& a, const T& b);
template <class T> const T& max (const T& a, const T& b);



extern "C"
void waveinitGPU(XBGPUParam Param);
void wavebnd(XBGPUParam Param);
void flowbnd(XBGPUParam Param, std::vector<SLBnd> slbnd, std::vector<WindBnd> wndbnd);
void wavestep(XBGPUParam Param);
void flowstep(XBGPUParam Param);
void sedimentstep(XBGPUParam Param);

void write_text_to_log_file( std::string text);
void SaveParamtolog(XBGPUParam XParam);

void readgridncsize(std::string ncfile, int &nx, int &ny, double &dx);
extern "C" void readnczb(int nx, int ny, std::string ncfile, DECNUM * &zb);



#endif