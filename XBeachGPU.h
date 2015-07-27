#ifndef XBeachGPU_H
#define XBeachGPU_H
 
#define pi 3.14159265

using DECNUM = float;



// additional functions
void makjonswap(DECNUM hm0gew,DECNUM fp,DECNUM mainang,DECNUM rt,DECNUM scoeff,DECNUM gam,DECNUM * theta,int ntheta,DECNUM& TTrep,DECNUM * &Stt);
extern "C" void creatncfile(char outfile[], int nx,int ny,/*int npart,*/DECNUM dx,DECNUM totaltime,int imodel,/*DECNUM * xxp,DECNUM * yyp,*/DECNUM *zb,DECNUM *zs,DECNUM * uu, DECNUM * vv, DECNUM * H,DECNUM * Tp,DECNUM * Dp,DECNUM * D,DECNUM * Urms,DECNUM * ueu,DECNUM * vev,DECNUM * C,DECNUM *Fx,DECNUM *Fy,DECNUM *hh,DECNUM *Hmean,DECNUM *uumean,DECNUM *vvmean,DECNUM *hhmean,DECNUM *zsmean,DECNUM *Cmean);
extern "C" void writestep2nc(char outfile[], int nx,int ny,/*int npart,*/DECNUM totaltime,int imodel/*,DECNUM *xxp,DECNUM *yyp*/,DECNUM *zb,DECNUM *zs,DECNUM * uu, DECNUM * vv, DECNUM * H, DECNUM * Tp, DECNUM *Dp,DECNUM *D,DECNUM *Urms,DECNUM *ueu,DECNUM *vev,DECNUM *C,DECNUM *dzb,DECNUM *Fx,DECNUM *Fy,DECNUM *hh,DECNUM *Hmean,DECNUM *uumean,DECNUM *vvmean,DECNUM *hhmean,DECNUM *zsmean,DECNUM *Cmean);

extern "C" void create3dnc(int nx,int ny,int nt,DECNUM dx,DECNUM totaltime,DECNUM *theta,DECNUM * var);
extern "C" void write3dvarnc(int nx,int ny,int nt,DECNUM totaltime,DECNUM * var);

extern "C" void read3Dnc(int nx, int ny,int ntheta,char ncfile[],DECNUM * &ee);
extern "C" void read2Dnc(int nx, int ny,char ncfile[],DECNUM * &hh);

extern "C" void readXbbndhead(char * wavebndfile,DECNUM &thetamin,DECNUM &thetamax,DECNUM &dtheta,DECNUM &dtwavbnd,int &nwavbnd,int &nwavfile);
extern "C" void readXbbndstep(int nx, int ny,int ntheta,char * wavebndfile,int step,DECNUM &Trep,double *&qfile,double *&Stfile );
extern "C" void readStatbnd(int nx, int ny,int ntheta,DECNUM rho,DECNUM g,char * wavebndfile,double *&Tpfile,double *&Stfile );
extern "C" void readbndhead(char * wavebndfile,DECNUM &thetamin,DECNUM &thetamax,DECNUM &dtheta,DECNUM &dtwavbnd,int &nwavbnd);

//Below is for teh new CPU routine
extern "C" int mminusC(int ix,int nx);
extern "C" int pplusC(int ix, int nx);
extern "C" int mminus2C(int ix,int nx);
extern "C" int pplus2C(int ix, int nx);
extern "C" int signC(DECNUM x);
extern "C" void set_bndCPU(int nx, int ny,DECNUM Trep,int ntheta,DECNUM * theta,DECNUM *&sigm);
extern "C" void dispersion_initCPU(int nx,int ny,DECNUM twopi,DECNUM g,DECNUM aphi,DECNUM bphi,DECNUM * sigm,DECNUM * hh,DECNUM * &cg);
extern "C" void updatezomCPU(int nx, int ny,DECNUM cf,DECNUM cf2,DECNUM fw,DECNUM fw2,DECNUM * structdepth, DECNUM * &cfm,DECNUM * &fwm);
extern "C" void ubndCPU(int nx, int ny, DECNUM dx, DECNUM dt,DECNUM g, DECNUM rho,DECNUM totaltime,DECNUM wavbndtime,DECNUM rt,DECNUM slbndtime, DECNUM rtsl,DECNUM zsbndold,DECNUM zsbndnew,DECNUM Trep,DECNUM * qbndold, DECNUM * qbndnew,DECNUM *&zs, DECNUM * &uu,DECNUM * &vv, DECNUM *vu, DECNUM * umean, DECNUM * vmean,DECNUM * zb,DECNUM * cg,DECNUM * hum, DECNUM * zo, DECNUM *Fx,DECNUM *&hh);
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


// General functions
//void CUDA_CHECK(cudaError CUDerr);

template <class T> const T& min (const T& a, const T& b);
template <class T> const T& max (const T& a, const T& b);



extern "C"
void waveinitGPU(void);
void wavebnd(void);
void flowbnd(void);
void wavestep(void);
void flowstep(void);
void sedimentstep(void);







#endif