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


extern "C" void creatncfile(XBGPUParam XParam, DECNUM totaltime, DECNUM *zb, DECNUM *zs, DECNUM * uu, DECNUM * vv, DECNUM * H, DECNUM * Tp, DECNUM * Dp, DECNUM * D, DECNUM * Urms, DECNUM * ueu, DECNUM * vev, DECNUM * C, DECNUM *Fx, DECNUM *Fy, DECNUM * hh, DECNUM *Hmean, DECNUM *uumean, DECNUM *vvmean, DECNUM *hhmean, DECNUM *zsmean, DECNUM *Cmean)
{               
	int status;
	int nx = XParam.nx;
	int ny = XParam.ny;
	double dx = XParam.dx;

   	int ncid,xx_dim,yy_dim,time_dim,p_dim;
	size_t nxx,nyy,nnpart;
	int  var_dimids[3], var_dimzb[2];
	
	int dzb_id,zb_id,zs_id,uu_id,vv_id,H_id,Tp_id,Dp_id,D_id,Urms_id,ueu_id,vev_id,time_id,xx_id,yy_id,xxp_id,yyp_id,Ceq_id;
	int Fx_id, Fy_id, hh_id,Hmean_id,uumean_id,vvmean_id,hhmean_id,zsmean_id,Cmean_id;
	nxx=nx;
	nyy=ny;
	//nnpart=npart;

	static size_t start[] = {0, 0, 0}; // start at first value 
    	static size_t count[] = {1, ny, nx};
	static size_t zbstart[] = {0, 0}; // start at first value 
    	static size_t zbcount[] = {ny, nx};
	//static size_t pstart[] = {0, 0}; // start at first value 
   // 	static size_t pcount[] = {1, npart};
	static size_t tst[]={0};
	static size_t xstart[] = {0}; // start at first value 
    	static size_t xcount[] = {nx};
	DECNUM *xval;
	static size_t ystart[] = {0}; // start at first value 
    	static size_t ycount[] = {ny};
	DECNUM *yval;

	xval=(DECNUM *)malloc(nx*sizeof(DECNUM));
	yval=(DECNUM *)malloc(ny*sizeof(DECNUM));

	for (int i=0; i<nx; i++)
	{
		xval[i]=i*dx;
	}
	for (int i=0; i<ny; i++)
	{
		yval[i]=i*dx;
	}



	//create the netcdf dataset
	status = nc_create(XParam.outfile.c_str(), NC_NOCLOBBER, &ncid);
	
	//Define dimensions: Name and length
	
	status = nc_def_dim(ncid, "xx", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	//status = nc_def_dim(ncid, "npart",nnpart,&p_dim);
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	int tdim[]={time_dim};
	int xdim[]={xx_dim};
	int ydim[]={yy_dim};
	//int pdim[2];
	//pdim[0]=time_dim;
	//pdim[1]=p_dim;
	//define variables: Name, Type,...
	var_dimids[0] = time_dim;
    	var_dimids[1] = yy_dim;
    	var_dimids[2] = xx_dim;
	var_dimzb[0]= yy_dim;
	var_dimzb[1]= xx_dim;
	
    	status = nc_def_var (ncid, "time", NC_FLOAT,1,tdim, &time_id);
	status = nc_def_var (ncid, "x", NC_FLOAT,1,xdim, &xx_id);
	status = nc_def_var (ncid, "y", NC_FLOAT,1,ydim, &yy_id);
	//status = nc_def_var (ncid, "xxp", NC_FLOAT,2,pdim, &xxp_id);
	//status = nc_def_var (ncid, "yyp", NC_FLOAT,2,pdim, &yyp_id);
	if(XParam.morphology == 1)
	{
		status = nc_def_var (ncid, "zb", NC_FLOAT, 3, var_dimids, &zb_id);
		status = nc_def_var (ncid, "dzb", NC_FLOAT, 3, var_dimids, &dzb_id);
	}
	else
	{
		status = nc_def_var (ncid, "zb", NC_FLOAT, 2, var_dimzb, &zb_id);
		
	}
		
	
    	status = nc_def_var (ncid, "zs", NC_FLOAT, 3, var_dimids, &zs_id);
    	status = nc_def_var (ncid, "uu", NC_FLOAT, 3, var_dimids, &uu_id);
    	status = nc_def_var (ncid, "vv", NC_FLOAT, 3, var_dimids, &vv_id);
    	status = nc_def_var (ncid, "H", NC_FLOAT, 3, var_dimids, &H_id);
		status = nc_def_var (ncid, "Fx", NC_FLOAT, 3, var_dimids, &Fx_id);
		status = nc_def_var (ncid, "Fy", NC_FLOAT, 3, var_dimids, &Fy_id);
	status = nc_def_var (ncid, "Tp", NC_FLOAT, 3, var_dimids, &Tp_id);
	status = nc_def_var (ncid, "Dp", NC_FLOAT, 3, var_dimids, &Dp_id);
	status = nc_def_var (ncid, "D", NC_FLOAT, 3, var_dimids, &D_id);
	status = nc_def_var (ncid, "Urms", NC_FLOAT, 3, var_dimids, &Urms_id);
	status = nc_def_var (ncid, "ueu", NC_FLOAT, 3, var_dimids, &ueu_id);
	status = nc_def_var (ncid, "vev", NC_FLOAT, 3, var_dimids, &vev_id);
	status = nc_def_var (ncid, "C", NC_FLOAT, 3, var_dimids, &Ceq_id);
	status = nc_def_var (ncid, "hh", NC_FLOAT, 3, var_dimids, &hh_id);
	status = nc_def_var (ncid, "Hmean", NC_FLOAT, 3, var_dimids, &Hmean_id);
	status = nc_def_var (ncid, "uumean", NC_FLOAT, 3, var_dimids, &uumean_id);
	status = nc_def_var (ncid, "vvmean", NC_FLOAT, 3, var_dimids, &vvmean_id);
	status = nc_def_var (ncid, "hhmean", NC_FLOAT, 3, var_dimids, &hhmean_id);
	status = nc_def_var (ncid, "zsmean", NC_FLOAT, 3, var_dimids, &zsmean_id);
	status = nc_def_var (ncid, "Cmean", NC_FLOAT, 3, var_dimids, &Cmean_id);
	//put attriute: assign attibute values
	//nc_put_att

	//End definitions: leave define mode
	status = nc_enddef(ncid);  

	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id,tst,&totaltime);
	status = nc_put_vara_float(ncid, xx_id,xstart,xcount, xval);
	status = nc_put_vara_float(ncid, yy_id,ystart,ycount, yval);
	//status = nc_put_vara_float(ncid, xxp_id,pstart,pcount, xxp);
	//status = nc_put_vara_float(ncid, yyp_id,pstart,pcount, yyp);
	if (XParam.morphology == 1)
	{
		status = nc_put_vara_float(ncid, zb_id, start, count, zb);
		status = nc_put_vara_float(ncid, dzb_id, start, count, H);
	}
	else {
		status = nc_put_vara_float(ncid, zb_id, zbstart, zbcount, zb);
	}

		
	
	status = nc_put_vara_float(ncid, zs_id, start, count, zs);
	status = nc_put_vara_float(ncid, uu_id, start, count, uu);
	status = nc_put_vara_float(ncid, vv_id, start, count, vv);
	status = nc_put_vara_float(ncid, H_id, start, count, H);
	status = nc_put_vara_float(ncid, Fx_id, start, count, Fx);
	status = nc_put_vara_float(ncid, Fy_id, start, count, Fy);
	status = nc_put_vara_float(ncid, Tp_id, start, count, Tp);
	status = nc_put_vara_float(ncid, Dp_id, start, count, Dp);
	status = nc_put_vara_float(ncid, D_id, start, count, D);
	status = nc_put_vara_float(ncid, Urms_id, start, count, Urms);
	status = nc_put_vara_float(ncid, ueu_id, start, count, ueu);
	status = nc_put_vara_float(ncid, vev_id, start, count, vev);
	status = nc_put_vara_float(ncid, Ceq_id, start, count, C);
	status = nc_put_vara_float(ncid, hh_id, start, count, hh);

	status = nc_put_vara_float(ncid, Hmean_id, start, count, Hmean);
	status = nc_put_vara_float(ncid, uumean_id, start, count, uumean);
	status = nc_put_vara_float(ncid, vvmean_id, start, count, vvmean);
	status = nc_put_vara_float(ncid, hhmean_id, start, count, hhmean);
	status = nc_put_vara_float(ncid, zsmean_id, start, count, zsmean);	
	status = nc_put_vara_float(ncid, Cmean_id, start, count, Cmean);

	//close and save new file
	status = nc_close(ncid);  
}
extern "C" void writestep2nc(XBGPUParam XParam, DECNUM totaltime,DECNUM *zb,DECNUM *zs,DECNUM * uu, DECNUM * vv, DECNUM * H,DECNUM * Tp,DECNUM * Dp,DECNUM * D,DECNUM * Urms,DECNUM *ueu,DECNUM * vev,DECNUM * C,DECNUM *dzb, DECNUM *Fx, DECNUM *Fy,DECNUM *hh,DECNUM *Hmean,DECNUM *uumean,DECNUM *vvmean,DECNUM *hhmean,DECNUM *zsmean,DECNUM *Cmean)
{
	int status;

	int nx = XParam.nx;
	int ny = XParam.ny;
	double dx = XParam.dx;

   	int ncid,time_dim,recid;
	size_t nxx,nyy;
	int dzb_id,zb_id,zs_id,uu_id,vv_id,H_id,Tp_id,Dp_id,D_id,Urms_id,ueu_id,vev_id,time_id,xxp_id,yyp_id,C_id,Fx_id,Fy_id,hh_id,Hmean_id,uumean_id,vvmean_id,hhmean_id,zsmean_id,Cmean_id;
	static size_t start[] = {0, 0, 0}; // start at first value 
    	static size_t count[] = {1, ny, nx};
 	//static size_t pstart[] = {0, 0}; // start at first value 
    //	static size_t pcount[] = {1, npart};
	static size_t tst[]={0};

	nxx=nx;
	nyy=ny;

	
	static size_t nrec;
	status = nc_open(XParam.outfile.c_str(), NC_WRITE, &ncid);
	
	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen  (ncid, recid, &nrec);
	printf("nrec=%d\n",nrec);

	//read file for variable ids
	status = nc_inq_varid(ncid, "time", &time_id);
	if(XParam.morphology == 1)
	{
		status = nc_inq_varid(ncid, "zb", &zb_id);
		status = nc_inq_varid(ncid, "dzb", &dzb_id);
	}
	status = nc_inq_varid(ncid, "zs", &zs_id);
	status = nc_inq_varid(ncid, "uu", &uu_id);
	status = nc_inq_varid(ncid, "vv", &vv_id);
	status = nc_inq_varid(ncid, "H", &H_id);
	status = nc_inq_varid(ncid, "Fx", &Fx_id);
	status = nc_inq_varid(ncid, "Fy", &Fy_id);
	status = nc_inq_varid(ncid, "Tp", &Tp_id);
	status = nc_inq_varid(ncid, "Dp", &Dp_id);
	status = nc_inq_varid(ncid, "D", &D_id);
	status = nc_inq_varid(ncid, "Urms", &Urms_id);
	status = nc_inq_varid(ncid, "ueu", &ueu_id);
	status = nc_inq_varid(ncid, "vev", &vev_id);
	status = nc_inq_varid(ncid, "C", &C_id);
	status = nc_inq_varid(ncid, "hh", &hh_id);

	status = nc_inq_varid(ncid, "Hmean", &Hmean_id);
	status = nc_inq_varid(ncid, "uumean", &uumean_id);
	status = nc_inq_varid(ncid, "vvmean", &vvmean_id);
	status = nc_inq_varid(ncid, "hhmean", &hhmean_id);
	status = nc_inq_varid(ncid, "zsmean", &zsmean_id);
	status = nc_inq_varid(ncid, "Cmean", &Cmean_id);
	//status = nc_inq_varid(ncid, "xxp", &xxp_id);
	//status = nc_inq_varid(ncid, "yyp", &yyp_id);
	
	
	start[0] = nrec;
	//pstart[0] = nrec;    
	tst[0]=nrec;

	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id,tst,&totaltime);
	if (XParam.morphology == 1)
	{
		status = nc_put_vara_float(ncid, zb_id, start, count, zb);
		status = nc_put_vara_float(ncid, dzb_id, start, count, dzb);
	}
	status = nc_put_vara_float(ncid, zs_id, start, count, zs);
	status = nc_put_vara_float(ncid, uu_id, start, count, uu);
	status = nc_put_vara_float(ncid, vv_id, start, count, vv);
	status = nc_put_vara_float(ncid, H_id, start, count, H);
	status = nc_put_vara_float(ncid, Fx_id, start, count, Fx);
	status = nc_put_vara_float(ncid, Fy_id, start, count, Fy);
	status = nc_put_vara_float(ncid, Tp_id, start, count, Tp);
	status = nc_put_vara_float(ncid, Dp_id, start, count, Dp);
	status = nc_put_vara_float(ncid, D_id, start, count, D);
	status = nc_put_vara_float(ncid, Urms_id, start, count, Urms);
	status = nc_put_vara_float(ncid, ueu_id, start, count, ueu);
	status = nc_put_vara_float(ncid, vev_id, start, count, vev);
	status = nc_put_vara_float(ncid, C_id, start, count, C);
	status = nc_put_vara_float(ncid, hh_id, start, count, hh);

	status = nc_put_vara_float(ncid, Hmean_id, start, count, Hmean);
	status = nc_put_vara_float(ncid, uumean_id, start, count, uumean);
	status = nc_put_vara_float(ncid, vvmean_id, start, count, vvmean);
	status = nc_put_vara_float(ncid, hhmean_id, start, count, hhmean);
	status = nc_put_vara_float(ncid, zsmean_id, start, count, zsmean);	
	status = nc_put_vara_float(ncid, Cmean_id, start, count, Cmean);
	//status = nc_put_vara_float(ncid, xxp_id, pstart, pcount, xxp);
	//status = nc_put_vara_float(ncid, yyp_id, pstart, pcount, yyp);


	//close and save
	status = nc_close(ncid);


}


extern "C" void create3dnc(int nx,int ny,int nt,DECNUM dx,DECNUM totaltime,DECNUM *theta,DECNUM * var)
{
	int status;
   	int ncid,xx_dim,yy_dim,time_dim,p_dim,tvar_id;
	
	size_t nxx,nyy,ntt;
	static size_t start[] = {0, 0, 0, 0}; // start at first value 
    static size_t count[] = {1, nt, ny, nx};
	int time_id,xx_id,yy_id,tt_id;	//
	nxx=nx;
	nyy=ny;
	ntt=nt;

	//create the netcdf dataset
	status = nc_create("3Dvar.nc", NC_NOCLOBBER, &ncid);
	
	//Define dimensions: Name and length
	
	status = nc_def_dim(ncid, "xx", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	status = nc_def_dim(ncid, "ntheta",ntt,&p_dim);
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	int tdim[]={time_dim};
	int xdim[]={xx_dim};
	int ydim[]={yy_dim};
	int pdim[]={p_dim};

	//define variables: Name, Type,...
	int  var_dimids[4];
	var_dimids[0] = time_dim;
    var_dimids[1] = p_dim;
	var_dimids[2] = yy_dim;
    var_dimids[3] = xx_dim;
	
	
    status = nc_def_var (ncid, "time", NC_FLOAT,1,tdim, &time_id);
	status = nc_def_var (ncid, "x", NC_FLOAT,1,xdim, &xx_id);
	status = nc_def_var (ncid, "y", NC_FLOAT,1,ydim, &yy_id);
	status = nc_def_var (ncid, "theta", NC_FLOAT,1,pdim, &tt_id);


	status = nc_def_var (ncid, "3Dvar", NC_FLOAT, 4, var_dimids, &tvar_id);


	status = nc_enddef(ncid); 


	static size_t tst[]={0};
	static size_t xstart[] = {0}; // start at first value 
    static size_t xcount[] = {nx};
	DECNUM *xval;
	static size_t ystart[] = {0}; // start at first value 
    static size_t ycount[] = {ny};
	DECNUM *yval;
	static size_t tstart[] = {0}; // start at first value 
    static size_t tcount[] = {nt};
	

	xval=(DECNUM *)malloc(nx*sizeof(DECNUM));
	yval=(DECNUM *)malloc(ny*sizeof(DECNUM));

	for (int i=0; i<nx; i++)
	{
		xval[i]=i*dx;
	}
	for (int i=0; i<ny; i++)
	{
		yval[i]=i*dx;
	}


	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id,tst,&totaltime);
	status = nc_put_vara_float(ncid, xx_id,xstart,xcount, xval);
	status = nc_put_vara_float(ncid, yy_id,ystart,ycount, yval);
	status = nc_put_vara_float(ncid, tt_id,tstart,tcount, theta);

	status = nc_put_vara_float(ncid, tvar_id, start, count, var);
	status = nc_close(ncid);

}
extern "C" void write3dvarnc(int nx,int ny,int nt,DECNUM totaltime,DECNUM * var)
{
	int status;
   	int ncid,time_dim,recid;
	size_t nxx,nyy;
	static size_t start[] = {0, 0, 0, 0}; // start at first value 
    static size_t count[] = {1, nt, ny, nx};
	static size_t tst[]={0};
	int time_id,var_id;


	nxx=nx;
	nyy=ny;

	
	static size_t nrec;
	status = nc_open("3Dvar.nc", NC_WRITE, &ncid);
	
	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen  (ncid, recid, &nrec);
	//printf("nrec=%d\n",nrec);

	//read file for variable ids
	status = nc_inq_varid(ncid, "time", &time_id);
	status = nc_inq_varid(ncid, "3Dvar", &var_id);
	
	start[0] = nrec;//
	tst[0]=nrec;

	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id,tst,&totaltime);
	status = nc_put_vara_float(ncid, var_id, start, count, var);
	status = nc_close(ncid);

}

extern "C" void read3Dnc(int nx, int ny,int ntheta,char ncfile[],DECNUM * &ee)
{
	int status;
	int ncid,ee_id;
	static size_t count[] = {nx, ny,ntheta};
	
	status = nc_open(ncfile, NC_NOWRITE, &ncid);
	status = nc_inq_varid(ncid, "z", &ee_id);
	status = nc_get_var_float(ncid,ee_id,ee);
	status = nc_close(ncid);
	
	
}

extern "C" void read2Dnc(int nx, int ny,char ncfile[],DECNUM * &hh)
{
	int status;
	int ncid,hh_id;
	static size_t count[] = {nx, ny};
	
	status = nc_open(ncfile, NC_NOWRITE, &ncid);
	status = nc_inq_varid(ncid, "hh", &hh_id);
	status = nc_get_var_float(ncid,hh_id,hh);
	status = nc_close(ncid);
	
	
}
