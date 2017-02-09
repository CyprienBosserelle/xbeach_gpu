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

void handle_error(int status) {
	if (status != NC_NOERR) {
		fprintf(stderr, "Netcdf %s\n", nc_strerror(status));
		write_text_to_log_file("Netcdf " + std::string(nc_strerror(status)));
		//fprintf(logfile, "Netcdf: %s\n", nc_strerror(status));
		exit(-1);
	}
}


extern "C" void creatncfileUD(XBGPUParam XParam, double totaltime,int ntheta,float dtheta,float  thetamin,float thetamax)
{
	int status;
	int nx = XParam.nx;
	int ny = XParam.ny;
	double dx = XParam.dx;
	size_t nxx, nyy, nth;
	int ncid, xx_dim, yy_dim, time_dim, theta_dim;
	float * xval, *yval, *thetaval;
	static size_t xcount[] = { nx };
	static size_t ycount[] = { ny };
	static size_t thcount[] = { ntheta };
	int time_id, xx_id, yy_id, th_id;

	static size_t tst[] = { 0 };
	static size_t xstart[] = { 0 }; // start at first value
	static size_t ystart[] = { 0 }; // start at first value 
	static size_t thstart[] = { 0 }; // start at first value
	nxx = nx;
	nyy = ny;
	nth = ntheta;

	//Recreat the x, y and theta array
	xval = (DECNUM *)malloc(nx*sizeof(DECNUM));
	yval = (DECNUM *)malloc(ny*sizeof(DECNUM));
	thetaval = (float *)malloc(ntheta*sizeof(float));

	for (int i = 0; i<nx; i++)
	{
		xval[i] = i*dx;
	}
	for (int i = 0; i<ny; i++)
	{
		yval[i] = i*dx;
	}
	for (int i = 0; i<ntheta; i++)
	{
		thetaval[i] = i*(dtheta)+thetamin + 0.5f*dtheta;
	}

	//create the netcdf datasetXParam.outfile.c_str()
	status = nc_create(XParam.outfile.c_str(), NC_NOCLOBBER, &ncid);

	//Define dimensions: Name and length

	status = nc_def_dim(ncid, "xx", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	status = nc_def_dim(ncid, "theta", nth, &theta_dim);
	//status = nc_def_dim(ncid, "npart",nnpart,&p_dim);
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);

	int tdim[] = { time_dim };
	int xdim[] = { xx_dim };
	int ydim[] = { yy_dim };
	int thdim[] = {theta_dim};
	
	status = nc_def_var(ncid, "time", NC_FLOAT, 1, tdim, &time_id);
	status = nc_def_var(ncid, "xx", NC_FLOAT, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "yy", NC_FLOAT, 1, ydim, &yy_id);
	status = nc_def_var(ncid, "theta", NC_FLOAT, 1, thdim, &th_id);
	//End definitions: leave define mode
	status = nc_enddef(ncid);

	//Provide values for variables
	status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_float(ncid, xx_id, xstart, xcount, xval);
	status = nc_put_vara_float(ncid, yy_id, ystart, ycount, yval);
	status = nc_put_vara_float(ncid, th_id, thstart, thcount, thetaval);
	//close and save new file
	status = nc_close(ncid);

	free(xval);
	free(yval);
	free(thetaval);
}

extern "C" void defncvar(XBGPUParam XParam, std::string varst, int vdim, float * var)
{
	int status;
	int ncid, var_id;
	int  var_dimid2D[2];
	int  var_dimid3D[3];
	int  var_dimid4D[4];

	int recid, xid, yid, thid;
	size_t ntheta;// nx and ny are stored in XParam not yet for ntheta

	status = nc_open(XParam.outfile.c_str(), NC_WRITE, &ncid);
	status = nc_redef(ncid);

	//Inquire dimensions ids
	status = nc_inq_unlimdim(ncid, &recid);//time
	status = nc_inq_dimid(ncid, "xx", &xid);
	status = nc_inq_dimid(ncid, "yy", &yid);
	status = nc_inq_dimid(ncid, "theta", &thid);
	status = nc_inq_dimlen(ncid, thid, &ntheta);

	var_dimid2D[0] = yid;
	var_dimid2D[1] = xid;

	var_dimid3D[0] = recid;
	var_dimid3D[1] = yid;
	var_dimid3D[2] = xid;

	var_dimid4D[0] = recid;
	var_dimid4D[1] = thid;
	var_dimid4D[2] = yid;
	var_dimid4D[3] = xid;

	static size_t start2D[] = { 0, 0 }; // start at first value 
	static size_t count2D[] = { XParam.ny, XParam.nx };

	static size_t start3D[] = { 0, 0, 0 }; // start at first value 
	static size_t count3D[] = { 1, XParam.ny, XParam.nx };

	static size_t start4D[] = { 0, 0, 0, 0 }; // start at first value 
	static size_t count4D[] = { 1, ntheta, XParam.ny, XParam.nx };


	if (vdim == 2)
	{
		status = nc_def_var(ncid, varst.c_str(), NC_FLOAT, vdim, var_dimid2D, &var_id);
		status = nc_enddef(ncid);
		status = nc_put_vara_float(ncid, var_id, start2D, count2D, var);
	}
	if (vdim == 3)
	{
		status = nc_def_var(ncid, varst.c_str(), NC_FLOAT, vdim, var_dimid3D, &var_id);
		status = nc_enddef(ncid);
		status = nc_put_vara_float(ncid, var_id, start3D, count3D, var);
	}
	if (vdim == 4)
	{
		status = nc_def_var(ncid, varst.c_str(), NC_FLOAT, vdim, var_dimid4D, &var_id);
		status = nc_enddef(ncid);
		status = nc_put_vara_float(ncid, var_id, start4D, count4D, var);
	}

	//close and save new file
	status = nc_close(ncid);

}


extern "C" void writenctimestep(XBGPUParam XParam, double totaltime)
{
	int status, ncid, recid, time_id;
	status = nc_open(XParam.outfile.c_str(), NC_WRITE, &ncid);
	static size_t nrec;
	static size_t tst[] = { 0 };
	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);
	status = nc_inq_varid(ncid, "time", &time_id);
	tst[0] = nrec;
	status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	//close and save
	status = nc_close(ncid);

}

extern "C" void writencvarstep(XBGPUParam XParam, std::string varst, float * var)
{
	int status, ncid,time_dim, recid, var_id,ndims;
	static size_t nrec;
	int dimids[NC_MAX_VAR_DIMS];
	size_t  *ddim, *start,*count;
//XParam.outfile.c_str()
	status = nc_open(XParam.outfile.c_str(), NC_WRITE, &ncid);

	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);

	status = nc_inq_varid(ncid, varst.c_str(), &var_id);

	status = nc_inq_varndims(ncid, var_id, &ndims);
	if (status != NC_NOERR) handle_error(status);
	//printf("hhVar:%d dims\n", ndimshh);

	status = nc_inq_vardimid(ncid, var_id, dimids);
	if (status != NC_NOERR) handle_error(status);

	ddim = (size_t *)malloc(ndims*sizeof(size_t));
	start = (size_t *)malloc(ndims*sizeof(size_t));
	count = (size_t *)malloc(ndims*sizeof(size_t));

	//Read dimensions nx_u ny_u 
	for (int iddim = 0; iddim < ndims; iddim++)
	{
		status = nc_inq_dimlen(ncid, dimids[iddim], &ddim[iddim]);
		if (status != NC_NOERR) handle_error(status);
		start[iddim] = 0;
		count[iddim] = ddim[iddim];
		//printf("dim:%d=%d\n", iddim, ddimhh[iddim]);
	}

	start[0] = nrec-1;
	count[0] = 1;

	status = nc_put_vara_float(ncid, var_id, start, count, var);
	//close and save
	status = nc_close(ncid);

	free(ddim);
	free(start);
	free(count);
}


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
	status = nc_def_var (ncid, "xx", NC_FLOAT,1,xdim, &xx_id);
	status = nc_def_var (ncid, "yy", NC_FLOAT,1,ydim, &yy_id);
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


	free(xval);
	free(yval);
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


extern "C" void create2dnc(int nx, int ny, double dx, double dy, double totaltime, double *xx, double *yy, double * var)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim, tvar_id;

	size_t nxx, nyy, ntt;
	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { 1, ny, nx };
	int time_id, xx_id, yy_id, tt_id;	//
	nxx = nx;
	nyy = ny;
	

	//create the netcdf dataset
	status = nc_create("2Dvar.nc", NC_NOCLOBBER, &ncid);

	//Define dimensions: Name and length

	status = nc_def_dim(ncid, "xx", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	int tdim[] = { time_dim };
	int xdim[] = { xx_dim };
	int ydim[] = { yy_dim };
	

	//define variables: Name, Type,...
	int  var_dimids[3];
	var_dimids[0] = time_dim;
	
	var_dimids[1] = yy_dim;
	var_dimids[2] = xx_dim;


	status = nc_def_var(ncid, "time", NC_DOUBLE, 1, tdim, &time_id);
	status = nc_def_var(ncid, "xx", NC_DOUBLE, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "yy", NC_DOUBLE, 1, ydim, &yy_id);



	status = nc_def_var(ncid, "2Dvar", NC_DOUBLE, 3, var_dimids, &tvar_id);


	status = nc_enddef(ncid);


	static size_t tst[] = { 0 };
	static size_t xstart[] = { 0 }; // start at first value 
	static size_t xcount[] = { nx };
	;
	static size_t ystart[] = { 0 }; // start at first value 
	static size_t ycount[] = { ny };






	//Provide values for variables
	status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_double(ncid, xx_id, xstart, xcount, xx);
	status = nc_put_vara_double(ncid, yy_id, ystart, ycount, yy);


	status = nc_put_vara_double(ncid, tvar_id, start, count, var);
	status = nc_close(ncid);

}

extern "C" void create3dnc(int nx, int ny, int nt, double dx, double dy, double dtheta,double totaltime, double *xx, double *yy,double *theta, double * var)
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
	
	
    status = nc_def_var (ncid, "time", NC_DOUBLE,1,tdim, &time_id);
	status = nc_def_var(ncid, "xx", NC_DOUBLE, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "yy", NC_DOUBLE, 1, ydim, &yy_id);
	status = nc_def_var(ncid, "theta", NC_DOUBLE, 1, pdim, &tt_id);


	status = nc_def_var(ncid, "3Dvar", NC_DOUBLE, 4, var_dimids, &tvar_id);


	status = nc_enddef(ncid); 


	static size_t tst[]={0};
	static size_t xstart[] = {0}; // start at first value 
    static size_t xcount[] = {nx};
	
	static size_t ystart[] = {0}; // start at first value 
    static size_t ycount[] = {ny};

	static size_t tstart[] = {0}; // start at first value 
    static size_t tcount[] = {nt};
	

	//Provide values for variables
	status = nc_put_var1_double(ncid, time_id,tst,&totaltime);
	status = nc_put_vara_double(ncid, xx_id,xstart,xcount, xx);
	status = nc_put_vara_double(ncid, yy_id,ystart,ycount, yy);
	status = nc_put_vara_double(ncid, tt_id,tstart,tcount, theta);

	status = nc_put_vara_double(ncid, tvar_id, start, count, var);
	status = nc_close(ncid);

}

extern "C" void write3dvarnc(int nx,int ny,int nt,double totaltime,double * var)
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
	status = nc_put_var1_double(ncid, time_id,tst,&totaltime);
	status = nc_put_vara_double(ncid, var_id, start, count, var);
	status = nc_close(ncid);

}

extern "C" void write2varnc(int nx, int ny,  double totaltime, double * var)
{
	int status;
	int ncid, time_dim, recid;
	size_t nxx, nyy;
	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { 1, ny, nx };
	static size_t tst[] = { 0 };
	int time_id, var_id;


	nxx = nx;
	nyy = ny;


	static size_t nrec;
	status = nc_open("3Dvar.nc", NC_WRITE, &ncid);

	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);
	//printf("nrec=%d\n",nrec);

	//read file for variable ids
	status = nc_inq_varid(ncid, "time", &time_id);
	status = nc_inq_varid(ncid, "3Dvar", &var_id);

	start[0] = nrec;//
	tst[0] = nrec;

	//Provide values for variables
	status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_double(ncid, var_id, start, count, var);
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

extern "C" void readnczb(int nx, int ny, std::string ncfile, DECNUM * &zb)
{
	int status;
	int ncid, hh_id;
	static size_t count[] = { nx, ny };

	status = nc_open(ncfile.c_str(), NC_NOWRITE, &ncid);
	status = nc_inq_varid(ncid, "zb", &hh_id);
	status = nc_get_var_float(ncid, hh_id, zb);
	status = nc_close(ncid);


}




void readgridncsize(std::string ncfile, int &nx, int &ny, double &dx)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, ndimshh, ndims;
	double *xcoord, *ycoord;
	int varid;


	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char coordname[NC_MAX_NAME + 1];
	size_t  *ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	
	//Open NC file
	printf("Open file\n");
	status = nc_open(ncfile.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//printf(" %s...\n", hhvar);
	status = nc_inq_varid(ncid, "zb", &varid);
	if (status != NC_NOERR)	handle_error(status);



	status = nc_inq_varndims(ncid, varid, &ndimshh);
	if (status != NC_NOERR) handle_error(status);
	//printf("hhVar:%d dims\n", ndimshh);

	status = nc_inq_vardimid(ncid, varid, dimids);
	if (status != NC_NOERR) handle_error(status);

	ddimhh = (size_t *)malloc(ndimshh*sizeof(size_t));

	//Read dimensions nx_u ny_u 
	for (int iddim = 0; iddim < ndimshh; iddim++)
	{
		status = nc_inq_dimlen(ncid, dimids[iddim], &ddimhh[iddim]);
		if (status != NC_NOERR) handle_error(status);

		//printf("dim:%d=%d\n", iddim, ddimhh[iddim]);
	}

	if (ndimshh > 2)
	{
		ny = ddimhh[1];
		nx = ddimhh[2];
	}
	else
	{
		ny = ddimhh[0];
		nx = ddimhh[1];
	}

	//allocate
	xcoord = (double *)malloc(nx*ny*sizeof(double));
	ycoord = (double *)malloc(nx*ny*sizeof(double));

	//inquire variable name for x dimension
	//aka x dim of hh
	int ycovar, xcovar;

	if (ndimshh > 2)
	{
		ycovar = dimids[1];
		xcovar = dimids[2];
	}
	else
	{
		ycovar = dimids[0];
		xcovar = dimids[1];
	}

	//ycoord
	status = nc_inq_dimname(ncid, ycovar, coordname);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, coordname, &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varndims(ncid, varid, &ndims);
	if (status != NC_NOERR) handle_error(status);

	if (ndims < 2)
	{
		double * ytempvar;
		ytempvar = (double *)malloc(ny*sizeof(double)); 
		size_t start[] = { 0 };
		size_t count[] = { ny };
		status = nc_get_vara_double(ncid, varid, start, count, ytempvar);
		if (status != NC_NOERR) handle_error(status);

		for (int i = 0; i<nx; i++)
		{
			for (int j = 0; j<ny; j++)
			{

				ycoord[i + j*nx] = ytempvar[j];

			}
		}
		free(ytempvar);
	}
	else
	{
		size_t start[] = { 0, 0 };
		size_t count[] = { ny, nx };
		status = nc_get_vara_double(ncid, varid, start, count, ycoord);
		if (status != NC_NOERR) handle_error(status);

	}
	//xcoord
	status = nc_inq_dimname(ncid, xcovar, coordname);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, coordname, &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varndims(ncid, varid, &ndims);
	if (status != NC_NOERR) handle_error(status);

	if (ndims < 2)
	{
		double * xtempvar;
		xtempvar = (double *)malloc(nx*sizeof(double));
		size_t start[] = { 0 };
		size_t count[] = { nx };
		status = nc_get_vara_double(ncid, varid, start, count, xtempvar);
		if (status != NC_NOERR) handle_error(status);

		for (int i = 0; i<nx; i++)
		{
			for (int j = 0; j<ny; j++)
			{

				xcoord[i + j*nx] = xtempvar[i];

			}
		}
		free(xtempvar);
	}
	else
	{
		size_t start[] = { 0, 0 };
		size_t count[] = { ny, nx };
		status = nc_get_vara_double(ncid, varid, start, count, xcoord);
		if (status != NC_NOERR) handle_error(status);

	}

	float dxx, dyy;
	//check dx
	dxx = abs(xcoord[0] - xcoord[nx - 1]) / (nx-1);
	dyy = abs(ycoord[0] - ycoord[(ny - 1)*nx]) / (ny-1);


	dx = dxx;


	status = nc_close(ncid);

	free(ddimhh);
	free(xcoord);
	free(ycoord);


}

void createbndnc(int tslen, int ny, int ntheta, double dy, double dtheta, double totaltime, double Hs, double Trep, double Tp, double Dp, double * timevec, double *yy, double *theta, double * ee, double * qx, double * qy)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim, ee_var_id, qx_var_id, qy_var_id;
	
	size_t nxx,nyy,ntt;
	static size_t start_ee[] = {0, 0, 0, 0}; // start at first value 
    static size_t count_ee[] = {1, ntheta, ny, tslen};
	static size_t start_q[] = { 0, 0, 0 }; // start at first value 
	static size_t count_q[] = { 1, ny, tslen };
	int time_id,xx_id,yy_id,tt_id;	//
	int Hs_id, Trep_id, Tp_id, Dp_id;
	nxx=tslen;
	nyy=ny;
	ntt=ntheta;

	//create the netcdf dataset
	status = nc_create("XBG_bnd_reuse.nc", NC_NOCLOBBER, &ncid);
	
	//Define dimensions: Name and length
	
	status = nc_def_dim(ncid, "tt", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	status = nc_def_dim(ncid, "theta",ntt,&p_dim);
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	int tdim[]={time_dim};
	int xdim[]={xx_dim};
	int ydim[]={yy_dim};
	int pdim[]={p_dim};

	//define variables: Name, Type,...
	int  var_dimids_ee[4];
	var_dimids_ee[0] = time_dim;
    var_dimids_ee[1] = p_dim;
	var_dimids_ee[2] = yy_dim;
    var_dimids_ee[3] = xx_dim;
	

	int var_dimids_q[3];
	var_dimids_q[0] = time_dim;
	var_dimids_q[1] = yy_dim;
	var_dimids_q[2] = xx_dim;
	
    status = nc_def_var (ncid, "time", NC_DOUBLE,1,tdim, &time_id);
	status = nc_def_var(ncid, "tt", NC_DOUBLE, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "yy", NC_DOUBLE, 1, ydim, &yy_id);
	status = nc_def_var(ncid, "theta", NC_DOUBLE, 1, pdim, &tt_id);
	status = nc_def_var(ncid, "Trep", NC_DOUBLE, 1, tdim, &Trep_id);
	status = nc_def_var(ncid, "Hs", NC_DOUBLE, 1, tdim, &Hs_id);
	status = nc_def_var(ncid, "Tp", NC_DOUBLE, 1, tdim, &Tp_id);
	status = nc_def_var(ncid, "Dp", NC_DOUBLE, 1, tdim, &Dp_id);


	status = nc_def_var(ncid, "ee_bnd", NC_DOUBLE, 4, var_dimids_ee, &ee_var_id);
	status = nc_def_var(ncid, "qx_bnd", NC_DOUBLE, 3, var_dimids_q, &qx_var_id);
	status = nc_def_var(ncid, "qy_bnd", NC_DOUBLE, 3, var_dimids_q, &qy_var_id);



	status = nc_enddef(ncid); 


	static size_t tst[]={0};
	static size_t xstart[] = {0}; // start at first value 
    static size_t xcount[] = {tslen};
	
	static size_t ystart[] = {0}; // start at first value 
    static size_t ycount[] = {ny};

	static size_t tstart[] = {0}; // start at first value 
    static size_t tcount[] = {ntheta};
	

	//Provide values for variables
	status = nc_put_var1_double(ncid, time_id,tst,&totaltime);
	status = nc_put_vara_double(ncid, xx_id, xstart, xcount, timevec);
	status = nc_put_vara_double(ncid, yy_id,ystart,ycount, yy);
	status = nc_put_vara_double(ncid, tt_id,tstart,tcount, theta);

	status = nc_put_var1_double(ncid, Trep_id, tst, &Trep);
	status = nc_put_var1_double(ncid, Tp_id, tst, &Tp);
	status = nc_put_var1_double(ncid, Hs_id, tst, &Hs);
	status = nc_put_var1_double(ncid, Dp_id, tst, &Dp);

	status = nc_put_vara_double(ncid, ee_var_id, start_ee, count_ee, ee);
	status = nc_put_vara_double(ncid, qx_var_id, start_q, count_q, qx);
	status = nc_put_vara_double(ncid, qy_var_id, start_q, count_q, qy);


	//close file
	status = nc_close(ncid);



}

void writebndnc(int tslen, int ny, int ntheta, double dy, double dtheta, double totaltime, double Hs, double Trep, double Tp, double Dp, double * timevec, double *yy, double *theta, double * ee, double * qx, double * qy)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim, ee_var_id, qx_var_id, qy_var_id, recid;

	size_t nxx, nyy, ntt;
	static size_t nrec;
	static size_t start_ee[] = { 0, 0, 0, 0 }; // start at first value 
	static size_t count_ee[] = { 1, ntheta, ny, tslen };
	static size_t start_q[] = { 0, 0, 0 }; // start at first value 
	static size_t count_q[] = { 1, ny, tslen };
	int time_id, xx_id, yy_id, tt_id;	//
	int Hs_id, Trep_id, Tp_id, Dp_id;
	nxx = tslen;
	nyy = ny;
	ntt = ntheta;

	//create the netcdf dataset
	status = nc_open("XBG_bnd_reuse.nc", NC_WRITE, &ncid);

	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);
	//printf("nrec=%d\n",nrec);

	//read file for variable ids
	status = nc_inq_varid(ncid, "time", &time_id);

	status = nc_inq_varid(ncid, "Hs", &Hs_id);
	status = nc_inq_varid(ncid, "Trep", &Trep_id);
	status = nc_inq_varid(ncid, "Tp", &Tp_id);
	status = nc_inq_varid(ncid, "Dp", &Dp_id);

	status = nc_inq_varid(ncid, "ee_bnd", &ee_var_id);
	status = nc_inq_varid(ncid, "qx_bnd", &qx_var_id);
	status = nc_inq_varid(ncid, "qy_bnd", &qy_var_id);

	start_ee[0] = nrec;//
	start_q[0] = nrec;//
	
	
	static size_t tst[] = { nrec };
	


	//Provide values for variables
	status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	
	status = nc_put_var1_double(ncid, Hs_id, tst, &Hs);
	status = nc_put_var1_double(ncid, Trep_id, tst, &Trep);
	status = nc_put_var1_double(ncid, Tp_id, tst, &Tp);
	status = nc_put_var1_double(ncid, Dp_id, tst, &Dp);

	status = nc_put_vara_double(ncid, ee_var_id, start_ee, count_ee, ee);
	status = nc_put_vara_double(ncid, qx_var_id, start_q, count_q, qx);
	status = nc_put_vara_double(ncid, qy_var_id, start_q, count_q, qy);


	//close file
	status = nc_close(ncid);



}

XBGPUParam read_reuse_bndnc_head(XBGPUParam Param)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, recid, theta_dimid, theta_varid, tt_dimid, tt_varid;
	
	static size_t nrec, ntheta, ntt;
	status = nc_open(Param.wavebndfile.c_str(), NC_NOWRITE, &ncid);

	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimid(ncid, "theta", &theta_dimid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimlen(ncid, theta_dimid, &ntheta);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_varid(ncid, "theta", &theta_varid);
	if (status != NC_NOERR)	handle_error(status);

	double thetamin, thetamax, dtheta;


	size_t start[] = { 0 };
	size_t count[] = { 1 };
	status = nc_get_vara_double(ncid, theta_varid, start, count, &thetamin);
	start[0] = ntheta-1;
	status = nc_get_vara_double(ncid, theta_varid, start, count, &thetamax);

	Param.ntheta = ntheta;
	Param.thetamin = thetamin;
	Param.thetamax = thetamax;
	
	Param.dtheta = (Param.thetamax - Param.thetamin) / Param.ntheta;

	status = nc_inq_dimid(ncid, "tt", &tt_dimid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimlen(ncid, tt_dimid, &ntt);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_varid(ncid, "tt", &tt_varid);
	if (status != NC_NOERR)	handle_error(status);

	double rtl,rtlp;

	start[0] = ntt - 1;
	status = nc_get_vara_double(ncid, tt_varid, start, count, &rtl);
	start[0] = ntt - 2;
	status = nc_get_vara_double(ncid, tt_varid, start, count, &rtlp);
	
	
	//close file
	status = nc_close(ncid);

	Param.rtlength = rtl;//ntt*(rtl - rtlp); // Not totally sure why this is bigger than rtl...
	Param.dtbc = rtl - rtlp;
	
	return Param;
}

std::vector<Wavebndparam> read_reuse_bndnc_vec(XBGPUParam Param)
{
	std::vector<Wavebndparam> wavebndvec;

	Wavebndparam wavebndline;

	//read the dimentions of grid, levels and time 
	int status;
	int ncid, recid, theta_dimid, theta_varid, tt_dimid, tt_varid;

	static size_t nrec, ntheta, ntt;
	status = nc_open(Param.wavebndfile.c_str(), NC_NOWRITE, &ncid);

	//read id from time dimension
	status = nc_inq_unlimdim(ncid, &recid);
	status = nc_inq_dimlen(ncid, recid, &nrec);
	if (status != NC_NOERR)	handle_error(status);
	
	//close file
	status = nc_close(ncid);

	//rtlength has just been calculated by teh previous function

	for (int i = 0; i < (nrec+1); i++)
	{
		wavebndline.time = (Param.rtlength)*i;
		
		//slbndline = readBSHline(line);
		wavebndvec.push_back(wavebndline);
	}


	return wavebndvec;
}


void read_reuse_bndnc(XBGPUParam Param, int rec, float &Trep, double * &qfile, double * &Stfile)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, theta_dimid, tt_dimid, yy_dimid, ee_varid,qx_varid, qy_varid;

	static size_t nrec, ntheta, ntt,nyy; 
	status = nc_open(Param.wavebndfile.c_str(), NC_NOWRITE, &ncid);


	status = nc_inq_dimid(ncid, "theta", &theta_dimid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimlen(ncid, theta_dimid, &ntheta);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimid(ncid, "tt", &tt_dimid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimlen(ncid, tt_dimid, &ntt);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimid(ncid, "yy", &yy_dimid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_dimlen(ncid, yy_dimid, &nyy);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_varid(ncid, "ee_bnd", &ee_varid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_varid(ncid, "qx_bnd", &qx_varid);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_inq_varid(ncid, "qy_bnd", &qy_varid);
	if (status != NC_NOERR)	handle_error(status);

	int tslenbc = ceil(Param.rtlength / Param.dtbc)+1; // should be equal to ntt!

	//Sanity check
	if (nyy > Param.ny)
	{
		// Error
		//exit
	}
	if (ntt > tslenbc)
	{
		// Error
		//exit
	}
	
	double * qxtemp, *qytemp, *eetemp;
	//Allocate the temporary array
	qxtemp = (double *)malloc(nyy*ntt*sizeof(double));
	qytemp = (double *)malloc(nyy*ntt*sizeof(double));
	eetemp = (double *)malloc(nyy*ntheta*ntt*sizeof(double));
	
	static size_t start_ee[] = { rec, 0, 0, 0 }; // start at first value 
	static size_t count_ee[] = { 1, ntheta, nyy, ntt };
	static size_t start_q[] = { rec, 0, 0 }; // start at first value 
	static size_t count_q[] = { 1, nyy, ntt };

	status = nc_get_vara_double(ncid, ee_varid, start_ee, count_ee, eetemp);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_get_vara_double(ncid, qx_varid, start_q, count_q, qxtemp);
	if (status != NC_NOERR)	handle_error(status);

	status = nc_get_vara_double(ncid, qy_varid, start_q, count_q, qytemp);
	if (status != NC_NOERR)	handle_error(status);

	
	//close file
	status = nc_close(ncid);

	for (int j = 0; j < nyy; j++)
	{
		for (int m = 0; m < ntt; m++)
		{
			 qfile[j + 0 * nyy + m*nyy * 4]=qxtemp[m + j*ntt] ;
			 qfile[j + 1 * nyy + m*nyy * 4]=qytemp[m + j*ntt] ;

			for (int itheta = 0; itheta < Param.ntheta; itheta++)
			{
				  Stfile[j + itheta*nyy + m*nyy*Param.ntheta] = eetemp[m + j*ntt + itheta*nyy*ntt];
			}
		}
	}

	Trep = 15.0; /// Temporary debug

	free(qxtemp);
	free(qytemp);
	free(eetemp);

}