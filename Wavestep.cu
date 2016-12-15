#include "XBeachGPU.h"


void CUDA_CHECK(cudaError CUDerr)
{


	if (cudaSuccess != CUDerr) {

		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \

			__FILE__, __LINE__, cudaGetErrorString(CUDerr));

		exit(EXIT_FAILURE);

	}
}





XBGPUParam waveinitGPU(XBGPUParam Param, std::vector<Wavebndparam> wavebnd)
{
	// Initialize wave model
	int nx, ny;
	nx = Param.nx;
	ny = Param.ny;


	
	if (Param.dtheta > 0.0)
	{
		Param.ntheta = round((Param.thetamax - Param.thetamin) / Param.dtheta);
	}
	else
	{
		if (Param.ntheta == 0)
		{
			Param.ntheta = 1;
		}
		Param.dtheta = (Param.thetamax - Param.thetamin) / Param.ntheta;
	}

	
	ntheta = Param.ntheta;
	dtheta = Param.dtheta;


	printf("ntheta=%d\tdtheta=%f\n", ntheta, dtheta);
	write_text_to_log_file("ntheta=" + std::to_string(ntheta) + "\t dtheta=" + std::to_string(dtheta));
	//printf("nwavbnd=%d\n", nwavbnd);

	if (Param.wavebndtype == 1)
	{
		nwavbnd = wavebnd.size(); // one Stfile/qfile will be used throughout the simulation
	}
	if (Param.wavebndtype >= 2)
	{
		nwavbnd = ceil(Param.rtlength / Param.dtbc);
	}
		
	theta = (DECNUM *)malloc(ntheta*sizeof(DECNUM));

	Stfile = (double *)malloc(ntheta*ny*nwavbnd*sizeof(double));
	qfile = (double *)malloc(4 * ny*nwavbnd*sizeof(double));
	Tpfile = (double *)malloc(nwavbnd*sizeof(double));
	//dummy=(double *)malloc(1000*sizeof(double));

	qbndnew = (DECNUM *)malloc(4 * ny*sizeof(DECNUM));
	qbndold = (DECNUM *)malloc(4 * ny*sizeof(DECNUM));
	St = (DECNUM *)malloc(ntheta*ny*sizeof(DECNUM));
	Stold = (DECNUM *)malloc(ntheta*ny*sizeof(DECNUM));
	Stnew = (DECNUM *)malloc(ntheta*ny*sizeof(DECNUM));
	cxsth = (DECNUM *)malloc(ntheta*sizeof(DECNUM));
	sxnth = (DECNUM *)malloc(ntheta*sizeof(DECNUM));

	for (int i = 0; i < ntheta; i++)
	{
		theta[i] = i*(Param.dtheta)+Param.thetamin + 0.5f*Param.dtheta;
		cxsth[i] = cos(theta[i]);
		sxnth[i] = sin(theta[i]);

		//printf("theta=%f\tcxsth=%f\tsxnth=%f\n", theta[i], cxsth[i], sxnth[i]);
	}

	dang = theta[1] - theta[0];
	//dtheta=dang;


	ee = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	dd = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	wete = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));

	rr = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));

	cgx = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cgy = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cx = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	cy = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	ctheta = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	thet = (DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));
	//drr=(DECNUM *)malloc(nx*ny*ntheta*sizeof(DECNUM));


	printf("Generating initial wave bnd data\n");
	write_text_to_log_file("Generating initial wave bnd data");
	if (Param.wavebndtype == 1)
	{
		//GenCstWave(Param, wavebnd, Stfile, qfile, Tpfile);
		GenCstWave(Param, wavebnd, theta, Stfile, qfile, Tpfile);
		Trep = Tpfile[0];
		//readStatbnd(nx, ny, ntheta, Param.rho, Param.g, Param.wavebndfile.c_str(), Tpfile, Stfile);
		//Trepold = Tpfile[0];
		//Trepnew = Tpfile[1];
		//rt = dtwavbnd;

	}

	if (Param.wavebndtype == 2)
	{
		//readXbbndstep(nx, ny, ntheta, Param.wavebndfile.c_str(), 1, Trepold, qfile, Stfile);
		readXbbndstep(Param, wavebnd, 0, Trep, qfile, Stfile);
				
	}
	if (Param.wavebndtype == 3)
	{
		//JONSWAP
		//First generate a Highres 2D spec
		double * HRfreq;
		double * HRdir;
		double * HRSpec;

		int nfHR, ndHR;

		makjonswap(Param, wavebnd, 0, nfHR, ndHR, HRfreq, HRdir, HRSpec);
		//Then generate wave group timeseries based on that spectra
		//GenWGnLBW(Param,HRfreq,HRdir,HRSpec,Trep,qfile,Stfile);



	}
	nwbndstep = 0;
	for (int ni = 0; ni < ny; ni++)
	{
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			Stold[ni + itheta*ny] = Stfile[ni + itheta*ny + nwbndstep*ny*ntheta];
			Stnew[ni + itheta*ny] = Stfile[ni + itheta*ny + (nwbndstep + 1)*ny*ntheta];


		}
		for (int xi = 0; xi < 4; xi++)
		{
			qbndold[ni + xi*ny] = qfile[ni + xi*ny + nwbndstep*ny * 4];
			qbndnew[ni + xi*ny] = qfile[ni + xi*ny + (nwbndstep + 1)*ny * 4];
		}
	}




	//fscanf(fwav,"%f\t%f\t%f\t%f\t%f\t%f",&hm0gew,&fp,&mainang,&scoeff,&gam,&rt);
	//mainang=(1.5*pi-grdalpha)-mainang*pi/180;
	//fp=1/fp;
	//printf("init rt=%f\n",rt);


	//makjonswap(hm0gew,fp,mainang,rt,scoeff,gam,theta,ntheta,Trepnew, Stnew);


	//Clac Stat

	
	for (int i = 0; i < ntheta; i++)                             //! Fill St
	{
		//St[i]=Stold[i];
		//printf("St[%d]=%f\n",i,St[i]);
		for (int ii = 0; ii < ny; ii++)
		{
			St[ii + i*ny] = Stold[ii + i*ny];
		}
	}









	//printf("hh=%f\n",hh[0]);



	// Apply bnd on CPU side

	for (int ii = 0; ii < nx; ii++)
	{

		for (int jj = 0; jj < ny; jj++)
		{

			for (int nt = 0; nt < ntheta; nt++)
			{
				if (ii == 0)
				{
					ee[0 + jj*nx + nt*nx*ny] = St[jj + nt*ny];// not on gpu since it is a bank conflicting problem
				}
				else
				{
					ee[ii + jj*nx + nt*nx*ny] = 0.0f;
				}
				rr[ii + jj*nx + nt*nx*ny] = 0.0f;

			}
		}
	}


	//run dispersion relation	

	return Param;

}



void wavebndOLD(XBGPUParam Param)
{
	int nx, ny;
	nx = Param.nx;
	ny = Param.ny;

	if (totaltime >= dtwavbnd*(nwavbnd*wxstep - 1))//The -1 here is so that we read the next file before the last step of the previous file runs out
	{
		if (Param.wavebndtype == 2)
		{
			//readXbbndstep(nx, ny, ntheta, Param.wavebndfile.c_str(), wxstep, Trep, qfile, Stfile);

		}
		nwbndstep = 0;



		for (int ni = 0; ni < ny; ni++)
		{

			for (int itheta = 0; itheta < ntheta; itheta++)
			{
				Stold[ni + itheta*ny] = Stfile[ni + itheta*ny + nwbndstep*ny*ntheta];
				Stnew[ni + itheta*ny] = Stfile[ni + itheta*ny + (nwbndstep + 1)*ny*ntheta];
			}
			if (Param.wavebndtype == 2)
			{
				for (int xi = 0; xi < 4; xi++)
				{
					qbndold[ni + xi*ny] = qfile[ni + xi*ny + nwbndstep*ny * 4];
					qbndnew[ni + xi*ny] = qfile[ni + xi*ny + (nwbndstep + 1)*ny * 4];
				}
			}
		}


		wxstep = wxstep + 1;

	}






	//if ((nstep==1 || nstep==nwstp) && (imodel==1 || imodel>=3))
	//{
	//update wave bnd

	if (totaltime >= wavbndtime /*&& wavebndtype==2*/)
	{
		for (int i = 0; i < ntheta; i++)                             //! Fill Stold
		{
			for (int ni = 0; ni < ny; ni++)
			{
				Stold[ni + i*ny] = Stfile[ni + i*ny + nwbndstep*ntheta*ny];

			}
		}
		//fscanf(fwav,"%f\t%f\t%f\t%f\t%f\t%f",&hm0gew,&fp,&mainang,&scoeff,&gam,&rt);
		//mainang=(1.5*pi-grdalpha)-mainang*pi/180;
		//printf("rt=%f\n",rt);

		//fp=1/fp;

		//fscanf(fwav,"%f\t%f",&rt,&Trepnew);
		nwbndstep = nwbndstep + 1;
		for (int i = 0; i < ntheta; i++)                             //! Fill St
		{
			for (int ni = 0; ni < ny; ni++)
			{


				Stnew[ni + i*ny] = Stfile[ni + i*ny + nwbndstep*ntheta*ny];

				if (Param.wavebndtype == 1)
				{

					Trep = Tpfile[nwbndstep];
				}
			}
		}
		if (Param.wavebndtype == 2)
		{
			for (int ni = 0; ni < ny; ni++)
			{
				for (int xi = 0; xi < 4; xi++)
				{
					qbndold[ni + xi*ny] = qbndnew[ni + ny*xi];
					qbndnew[ni + xi*ny] = qfile[ni + xi*ny + nwbndstep*ny * 4];
				}
			}
			if (Param.GPUDEVICE >= 0)
			{
				CUDA_CHECK(cudaMemcpy(qbndold_g, qbndold, 4 * ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
				CUDA_CHECK(cudaMemcpy(qbndnew_g, qbndnew, 4 * ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
			}

			//printf("qbndold[300]=%f\n",qbndold[300]);
			//printf("qbndnew[300]=%f\n",qbndnew[300]);
		}
		//printf("Stfile[0]=%f\n",Stfile[0]);


		//makjonswap(hm0gew,fp,mainang,rt,scoeff,gam,theta,ntheta,Trepnew, Stnew);
		//wavbndtime=wavbndtime+dtwavbnd;
		wavbndtime = nwbndstep*dtwavbnd + (wxstep - 1)*nwavbnd*dtwavbnd; //should be better than above as it will not accumulate the rounbd off error
	}

	for (int i = 0; i < ntheta; i++)                             //! Fill St
	{
		for (int ni = 0; ni < ny; ni++)
		{
			St[ni + i*ny] = Stold[ni + i*ny] + (totaltime - wavbndtime + dtwavbnd)*(Stnew[ni + i*ny] - Stold[ni + i*ny]) / dtwavbnd;
		}
		//printf("St[%d]=%f\n",i,St[i*ny]);
	}
	//printf("Wave timestep:%f\n",wdt);
	//Wave model step
	//wavestep();
	//nwstp = nstep + nstpw;
	//wdt = dt;
	//}

}

void wavebnd(XBGPUParam Param, std::vector<Wavebndparam> wavebndvec)
{
	int nx, ny;
	double timenext, timesincelast;
	nx = Param.nx;
	ny = Param.ny;
	//update sl bnd

	// find next timestep

	double difft = wavebndvec[WAVstepinbnd].time - totaltime;

	
	if (difft < 0.0)
	{
		WAVstepinbnd++;

		
		if (Param.wavebndtype == 2)
		{
			//Read new STfile and qfile
			readXbbndstep(Param, wavebndvec, WAVstepinbnd - 1, Trep, qfile, Stfile);
			
		}
	}

	if (Param.wavebndtype == 1)
	{
		nwbndstep = WAVstepinbnd - 1;
		timenext = wavebndvec[WAVstepinbnd].time - wavebndvec[WAVstepinbnd - 1].time;
		timesincelast = (totaltime - wavebndvec[WAVstepinbnd - 1].time);

	}
	if (Param.wavebndtype == 2)
	{
		nwbndstep = floor((totaltime - wavebndvec[WAVstepinbnd - 1].time) / Param.dtbc);
		timenext = Param.dtbc;
		timesincelast = totaltime - (nwbndstep*Param.dtbc + wavebndvec[WAVstepinbnd - 1].time);
	}



	for (int ni = 0; ni < ny; ni++)
	{
		for (int itheta = 0; itheta < ntheta; itheta++)
		{
			Stold[ni + itheta*ny] = Stfile[ni + itheta*ny + nwbndstep*ny*ntheta];
			Stnew[ni + itheta*ny] = Stfile[ni + itheta*ny + (nwbndstep + 1)*ny*ntheta];


		}
		for (int xi = 0; xi < 4; xi++)
		{
			qbndold[ni + xi*ny] = qfile[ni + xi*ny + nwbndstep*ny * 4];
			qbndnew[ni + xi*ny] = qfile[ni + xi*ny + (nwbndstep + 1)*ny * 4];
		}
	}

	for (int i = 0; i < ntheta; i++)                             //! Fill St
	{
		for (int ni = 0; ni < ny; ni++)
		{
			St[ni + i*ny] = interptime(Stnew[ni + i*ny], Stold[ni + i*ny], timenext, timesincelast);
		}
		
	}

	if (Param.flow == 1)
	{
		CUDA_CHECK(cudaMemcpy(qbndold_g, qbndold, 4 * ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(qbndnew_g, qbndnew, 4 * ny*sizeof(DECNUM), cudaMemcpyHostToDevice));
	}

}

void wavestep(XBGPUParam Param)
{
	int nx, ny;
	nx = Param.nx;
	ny = Param.ny;
	double dt = Param.dt;
	ntheta = Param.ntheta;
	//Subroutine runs the wave model

	dim3 blockDim(16, 16, 1);
	dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);



	dim3 blockDim4(4, 4, 1);
	dim3 gridDim4(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);


	CUDA_CHECK(cudaMemcpy(St_g, St, ny*ntheta*sizeof(DECNUM), cudaMemcpyHostToDevice));
	//offshorebndWav(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g)
	offshorebndWav << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, totaltime, Trep, St_g, sigm_g, ee_g);
	//CUT_CHECK_ERROR("Offshore Wave bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//Sanity check
	sanity << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, hh_g, sigm_g, ntheta, ee_g);

	//CUT_CHECK_ERROR("sanity execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//CUDA_CHECK( cudaMalloc((void **)&cg_g, nx*ny*sizeof(DECNUM )) );
	//CUDA_CHECK( cudaMalloc((void **)&cx_g, nx*ny*ntheta*sizeof(DECNUM )) );
	//	CUDA_CHECK( cudaMalloc((void **)&c_g, nx*ny*sizeof(DECNUM )) );
	//CUDA_CHECK( cudaMalloc((void **)&cy_g, nx*ny*ntheta*sizeof(DECNUM )) );
	//CUDA_CHECK( cudaMalloc((void **)&k_g, nx*ny*sizeof(DECNUM )) );
	CUDA_CHECK(cudaMalloc((void **)&kh_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&sinh2kh_g, nx*ny*sizeof(DECNUM)));



	//dispersion
	dispersion << <gridDim, blockDim, 0 >> >(nx, ny, twopi, Param.g, aphi, bphi, sigm_g, hh_g, k_g, c_g, kh_g, sinh2kh_g, cg_g);
	//CUT_CHECK_ERROR("dispersion execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//CUDA_CHECK( cudaMemcpy(C,kh_g,  ny*nx*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );

	CUDA_CHECK(cudaMalloc((void **)&dhdx_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&dhdy_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&dudx_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&dudy_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&dvdx_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&dvdy_g, nx*ny*sizeof(DECNUM)));


	// Wave current interaction	(i.e remove wci in shallow water)
	calcwci << <gridDim, blockDim, 0 >> >(nx, ny, Param.wci, Param.hwci, hh_g, wci_g);
	//CUT_CHECK_ERROR("calcwci execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	// // Slopes of water depth and velocities
	slopes << <gridDim, blockDim, 0 >> >(nx, ny, Param.dx, hh_g, uu_g, vv_g, dhdx_g, dhdy_g, dudx_g, dudy_g, dvdx_g, dvdy_g);//
	//CUT_CHECK_ERROR("slopes execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//CUDA_CHECK( cudaMalloc((void **)&cgx_g, nx*ny*ntheta*sizeof(DECNUM )) );
	//CUDA_CHECK( cudaMalloc((void **)&cgy_g, nx*ny*ntheta*sizeof(DECNUM )) );
	//CUDA_CHECK( cudaMalloc((void **)&ctheta_g, nx*ny*ntheta*sizeof(DECNUM )) );
	//CUDA_CHECK( cudaMemcpy(C,kh_g,  ny*nx*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
	//Propagation speed in theta space
	propagtheta << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, wci_g, ctheta_g,/*c_g,cx_g,cy_g,*/cxsth_g, sxnth_g,/*uu_g,vv_g,*/dhdx_g, dhdy_g, dudx_g, dudy_g, dvdx_g, dvdy_g, sigm_g, kh_g);//
	//CUT_CHECK_ERROR("propagtheta execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//////////
	//CUDA_CHECK( cudaMemcpy(ctheta,ctheta_g,  ny*nx*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );
	//////////

	CUDA_CHECK(cudaFree(dhdx_g));
	CUDA_CHECK(cudaFree(dhdy_g));
	CUDA_CHECK(cudaFree(dudx_g));
	CUDA_CHECK(cudaFree(dudy_g));
	CUDA_CHECK(cudaFree(dvdx_g));
	CUDA_CHECK(cudaFree(dvdy_g));


	//


	//read3Dnc(nx,ny,ntheta,"eeX.nc",ee);
	//CUDA_CHECK( cudaMemcpy(ee_g, ee, nx*ny*ntheta*sizeof(DECNUM ), cudaMemcpyHostToDevice) );




	//
	// transform to wave action
	//
	action << <gridDim, blockDim, 0 >> >(ntheta, nx, ny, ee_g, sigm_g);
	//CUT_CHECK_ERROR("action execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());




	//
	// Upwind Euler timestep propagation
	//
	CUDA_CHECK(cudaMalloc((void **)&xadvec_g, nx*ny*ntheta*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&yadvec_g, nx*ny*ntheta*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&thetaadvec_g, nx*ny*ntheta*sizeof(DECNUM)));

	xadvecupwind2 << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, wci_g, ee_g, cg_g, cxsth_g, uu_g, xadvec_g);
	//CUT_CHECK_ERROR("eulerupwind xadvec execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	yadvecupwind2 << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, wci_g, ee_g, cg_g, sxnth_g, vv_g, yadvec_g);
	//CUT_CHECK_ERROR("eulerupwind yadvec execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	//CUDA_CHECK( cudaMalloc((void **)&eect_g, nx*ny*ntheta*sizeof(DECNUM )) );

	//eectheta<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,ee_g,ctheta_g,eect_g);
	////CUT_CHECK_ERROR("eulerupwind eectheta execution failed\n");
	//CUDA_CHECK( cudaThreadSynchronize() );

	//thetaadvecuw<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,eect_g,thetaadvec_g);
	////CUT_CHECK_ERROR("eulerupwind thetaadvecuw execution failed\n");
	//CUDA_CHECK( cudaThreadSynchronize() );

	thetaadvecuw2ho << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, Param.wci, ee_g, ctheta_g, thetaadvec_g);
	//CUT_CHECK_ERROR("eulerupwind thetaadvec execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());
	//CUDA_CHECK( cudaMemcpy(ctheta,yadvec_g,  ny*nx*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );

	//CUDA_CHECK( cudaMemcpy(ctheta,thetaadvec_g,  ny*nx*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );


	//read3Dnc(nx,ny,ntheta,"xadvecX.nc",ee);
	//CUDA_CHECK( cudaMemcpy(xadvec_g, ee, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read3Dnc(nx,ny,ntheta,"yadvecX.nc",ee);
	//CUDA_CHECK( cudaMemcpy(yadvec_g, ee, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );

	//read3Dnc(nx,ny,ntheta,"thetaadvecX.nc",ee);
	//CUDA_CHECK( cudaMemcpy(thetaadvec_g, ee, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );



	eulerupwind << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, Param.wci, ee_g, xadvec_g, yadvec_g, thetaadvec_g);
	//CUT_CHECK_ERROR("eulerupwind  execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());






	//CUDA_CHECK( cudaFree(cgx_g));
	//CUDA_CHECK( cudaFree(cgy_g));
	//Fix lateraL BND
	rollerlatbnd << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.eps, hh_g, ee_g);
	//CUT_CHECK_ERROR("energy latbnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//
	// transform back to wave energy
	//
	energy << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, ee_g, sigm_g);
	//CUT_CHECK_ERROR("energy execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//CUDA_CHECK( cudaMemcpy(ctheta,ee_g,  ny*nx*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );

	//CUDA_CHECK( cudaMalloc((void **)&H_g, nx*ny*sizeof(DECNUM )) );
	CUDA_CHECK(cudaMalloc((void **)&E_g, nx*ny*sizeof(DECNUM)));
	//CUDA_CHECK( cudaMalloc((void **)&D_g, nx*ny*sizeof(DECNUM )) );

	//CUDA_CHECK( cudaMemcpy(ctheta,ee_g,  ny*nx*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );

	//
	// Energy integrated over wave directions,Hrms
	//
	energint << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.rho, Param.g, Param.gammax, E_g, H_g, hh_g, ee_g);
	//CUT_CHECK_ERROR("energint execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());




	//
	// calculate change in intrinsic frequency // removed because it is super slow and doesn't do much
	//
	// tm is thetamean and it is calculated in the mean dir scheme
	//	CUDA_CHECK( cudaMalloc((void **)&tm_g, nx*ny*sizeof(DECNUM )) );
	//	calctm<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,tm_g,theta_g,ee_g);
	//	//CUT_CHECK_ERROR("energint execution failed\n");
	//    CUDA_CHECK( cudaThreadSynchronize() );
	/*
	//Change of intrinsec frequency
	*/





	// 
	//  Total dissipation from breaking  and bottom friction
	//

	if (Param.breakmodel == 1)
	{
		roelvink << <gridDim, blockDim, 0 >> >(nx, ny, Param.rho, Param.g, Param.gammaa, Param.alpha, Param.n, Trep, fwm_g, cfm_g, hh_g, H_g, E_g, D_g, k_g);
		//CUT_CHECK_ERROR("roelvink execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());
	}
	else
	{
		baldock << <gridDim, blockDim, 0 >> > (nx, ny, Param.rho, Param.g, Param.gammaa, Param.alpha, Param.n, Trep, fwm_g, cfm_g, hh_g, H_g, E_g, D_g, k_g);//Baldock more appropriate for pseudo stationary cases
		//CUT_CHECK_ERROR("baldoc execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());
	}
	//
	//  Calculate roller energy balance
	//
	//CUDA_CHECK( cudaMemcpy(hhmean,E_g, nx*ny*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );

	if (Param.roller == 1)
	{
		xadvecupwind2 << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, wci_g, rr_g, c_g, cxsth_g, uu_g, xadvec_g);
		//CUT_CHECK_ERROR("eulerupwind xadvec execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		yadvecupwind2 << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, wci_g, rr_g, c_g, sxnth_g, vv_g, yadvec_g);
		//CUT_CHECK_ERROR("eulerupwind yadvec execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		//eectheta<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,rr_g,ctheta_g,eect_g);
		////CUT_CHECK_ERROR("eulerupwind eectheta execution failed\n");
		//CUDA_CHECK( cudaThreadSynchronize() );

		//thetaadvecuw<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,dtheta,eect_g,thetaadvec_g);
		////CUT_CHECK_ERROR("eulerupwind thetaadvecuw execution failed\n");
		//CUDA_CHECK( cudaThreadSynchronize() );	

		thetaadvecuw2ho << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, Param.wci, rr_g, ctheta_g, thetaadvec_g);
		//CUT_CHECK_ERROR("eulerupwind thetaadvec execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		eulerupwind << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.dx, dt, Param.wci, rr_g, xadvec_g, yadvec_g, thetaadvec_g);
		//CUT_CHECK_ERROR("eulerupwind  execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());

		//
		//  Adjust lateral bnds
		//

		rollerlatbnd << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.eps, hh_g, rr_g);
		//CUT_CHECK_ERROR("rollerlatbnd execution failed\n");
		CUDA_CHECK(cudaThreadSynchronize());
	}
	//CUDA_CHECK( cudaFree(eect_g));
	CUDA_CHECK(cudaFree(xadvec_g));
	CUDA_CHECK(cudaFree(yadvec_g));
	CUDA_CHECK(cudaFree(thetaadvec_g));

	//read2Dnc(nx,ny,"D.nc",uu);
	//CUDA_CHECK( cudaMemcpy(D_g, uu, nx*ny*sizeof(DECNUM ), cudaMemcpyHostToDevice) );


	// 
	//  Distribution of dissipation over directions and frequencies
	//                               
	dissipation << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.eps, dt, Param.g, Param.beta, wci_g, hh_g, ee_g, D_g, E_g, rr_g, c_g, cxsth_g, sxnth_g, uu_g, vv_g, DR_g, R_g);
	//CUT_CHECK_ERROR("dissipation execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//
	//Fix lateraL BND
	//
	rollerlatbnd << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.eps, hh_g, ee_g);
	//CUT_CHECK_ERROR("energy latbnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	// 
	//  Compute mean wave direction
	// 

	meandir << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.rho, Param.g, dtheta, ee_g, theta_g, thetamean_g, E_g, H_g);
	//CUT_CHECK_ERROR("meandir execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());



	//
	//	Constant warm start // WARNING ONLY TO BE USED FOR DEBUGGING
	//
	//	read3Dnc(nx,ny,ntheta,"eeX.nc",ee);
	//	CUDA_CHECK( cudaMemcpy(ee_g, ee, nx*ny*ntheta*sizeof(DECNUM ), cudaMemcpyHostToDevice) );





	//
	// Radiation stresses and forcing terms
	//

	CUDA_CHECK(cudaMalloc((void **)&Sxx_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Sxy_g, nx*ny*sizeof(DECNUM)));
	CUDA_CHECK(cudaMalloc((void **)&Syy_g, nx*ny*sizeof(DECNUM)));

	radstress << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.dx, dtheta, ee_g, rr_g, cxsth_g, sxnth_g, cg_g, c_g, Sxx_g, Sxy_g, Syy_g);

	//CUT_CHECK_ERROR("radstress execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//	
	// Wave forces
	//
	wavforce << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, Param.dx, dtheta, Sxx_g, Sxy_g, Syy_g, Fx_g, Fy_g, hh_g);
	//CUT_CHECK_ERROR("wavforce execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	twodimbndnoix << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, hh_g, Fx_g);
	//CUT_CHECK_ERROR("wave force X bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	twodimbnd << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, hh_g, Fy_g);
	//CUT_CHECK_ERROR("wave force Y bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//CUDA_CHECK( cudaMemcpy(ctheta,ctheta_g,  ny*nx*ntheta*sizeof(DECNUM ), cudaMemcpyDeviceToHost) );


	//
	// CAlculate stokes velocity and breaker delay //Breaker delay removed because it is slow and kinda useless
	//
	breakerdelay << <gridDim, blockDim, 0 >> >(nx, ny, ntheta, dtheta, Param.g, Param.rho, Trep, Param.eps, urms_g, ust_g, H_g, E_g, c_g, k_g, hh_g, R_g);
	//CUT_CHECK_ERROR("breakerdelay execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());

	//twodimbnd<<<gridDim, blockDim, 0>>>(nx,ny,eps,hh_g,urms_g);
	//CUT_CHECK_ERROR("wave force Y bnd execution failed\n");
	//CUDA_CHECK( cudaThreadSynchronize() );

	twodimbnd << <gridDim, blockDim, 0 >> >(nx, ny, Param.eps, hh_g, ust_g);
	//CUT_CHECK_ERROR("wave force Y bnd execution failed\n");
	CUDA_CHECK(cudaThreadSynchronize());


	CUDA_CHECK(cudaFree(Sxy_g));
	CUDA_CHECK(cudaFree(Sxx_g));
	CUDA_CHECK(cudaFree(Syy_g));
	//CUDA_CHECK( cudaFree(cg_g));
	//CUDA_CHECK( cudaFree(c_g));
	CUDA_CHECK(cudaFree(tm_g));




	//
	// Adjust Offshore Bnd
	//

	//CUDA_CHECK( cudaMemcpy(St_g, St, ny*ntheta*sizeof(DECNUM ), cudaMemcpyHostToDevice) );
	//offshorebndWav(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g)
	//offshorebndWav<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g);
	////CUT_CHECK_ERROR("Offshore Wave bnd execution failed\n");
	//CUDA_CHECK( cudaThreadSynchronize() );



	CUDA_CHECK(cudaFree(E_g));
	//CUDA_CHECK( cudaFree(H_g));
	//CUDA_CHECK( cudaFree(D_g));

	//CUDA_CHECK( cudaFree(k_g));
	CUDA_CHECK(cudaFree(kh_g));
	CUDA_CHECK(cudaFree(sinh2kh_g));








}

void wavestepCPU(XBGPUParam Param)
{
	int nx, ny;
	nx = Param.nx;
	ny = Param.ny;
	//Subroutine runs the wave model

	//printf("%2.2f\t",ee_g[16+8*nx+1*ntheta]);
	//offshorebndWav(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g)
	offshorebndWavCPU(nx, ny, ntheta, totaltime, Trep, St, sigm_g, ee_g);


	//Sanity check
	sanityCPU(nx, ny, Param.eps, hh_g, sigm_g, ntheta, ee_g);
	//printf("%2.2f\t", ee_g[16 + 8 * nx + 1 * ntheta]);

	//dispersion
	//printf("%2.2f\t", ee_g[16 + 8 * nx + 1 * ntheta]);
	dispersionCPU(nx, ny, twopi, Param.g, aphi, bphi, sigm_g, hh_g, k_g, c_g, kh_g, sinh2kh_g, cg_g);


	// Wave current interaction	(i.e remove wci in shallow water)
	calcwciCPU(nx, ny, Param.wci, Param.hwci, hh_g, wci_g);
	//printf("%f\t",ee_g[0+16*nx+6*nx*ny]);

	// Slopes of water depth and velocities
	slopesCPU(nx, ny, Param.dx, hh_g, uu_g, vv_g, dhdx_g, dhdy_g, dudx_g, dudy_g, dvdx_g, dvdy_g);//
	//printf("%f\t",ee_g[0+16*nx+6*nx*ny]);

	//Propagation speed in theta space
	propagthetaCPU(nx, ny, ntheta, wci_g, ctheta_g, cxsth_g, sxnth_g, dhdx_g, dhdy_g, dudx_g, dudy_g, dvdx_g, dvdy_g, sigm_g, kh_g);//
	//printf("%f\n",ee_g[200+16*nx]);


	//read3Dnc(nx, ny, ntheta, "eeX.nc", ee_g);
	//printf("%f\n", ee_g[200 + 16 * nx]);
	// transform to wave action
	actionCPU(ntheta, nx, ny, ee_g, sigm_g);



	// Upwind Euler timestep propagation
	xadvecupwind2CPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, wci_g, ee_g, cg_g, cxsth_g, uu_g, xadvec_g);

	yadvecupwind2CPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, wci_g, ee_g, cg_g, sxnth_g, vv_g, yadvec_g);

	thetaadvecuw2hoCPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, Param.wci, ee_g, ctheta_g, thetaadvec_g);

	//Apply advection
	eulerupwindCPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, Param.wci, ee_g, xadvec_g, yadvec_g, thetaadvec_g);

	//Fix lateraL BND
	rollerlatbndCPU(nx, ny, ntheta, Param.eps, hh_g, ee_g);

	// transform back to wave energy
	energyCPU(nx, ny, ntheta, ee_g, sigm_g);

	// Energy integrated over wave directions,Hrms
	energintCPU(nx, ny, ntheta, dtheta, Param.rho, Param.g, Param.gammax, E_g, H_g, hh_g, ee_g);

	//
	// calculate change in intrinsic frequency // removed because it is super slow and doesn't do much
	//
	// tm is thetamean and it is calculated in the mean dir scheme
	//	CUDA_CHECK( cudaMalloc((void **)&tm_g, nx*ny*sizeof(DECNUM )) );
	//	calctm<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,tm_g,theta_g,ee_g);
	//	//CUT_CHECK_ERROR("energint execution failed\n");
	//    CUDA_CHECK( cudaThreadSynchronize() );

	//Change of intrinsec frequency


	//  Total dissipation from breaking  and bottom friction
	if (Param.breakmodel == 1)
	{
		roelvinkCPU(nx, ny, Param.rho, Param.g, Param.gammaa, Param.alpha, Param.n, Trep, fwm_g, cfm_g, hh_g, H_g, E_g, D_g, k_g);

	}
	else
	{
		baldockCPU(nx, ny, Param.rho, Param.g, Param.gammaa, Param.alpha, Param.n, Trep, fwm_g, cfm_g, hh_g, H_g, E_g, D_g, k_g);//Baldock more appropriate for pseudo stationary cases

	}

	//  Calculate roller energy balance
	if (roller == 1)
	{
		xadvecupwind2CPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, wci_g, rr_g, c_g, cxsth_g, uu_g, xadvec_g);

		yadvecupwind2CPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, wci_g, rr_g, c_g, sxnth_g, vv_g, yadvec_g);

		thetaadvecuw2hoCPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, Param.wci, rr_g, ctheta_g, thetaadvec_g);

		eulerupwindCPU(nx, ny, ntheta, dtheta, Param.dx, Param.dt, Param.wci, rr_g, xadvec_g, yadvec_g, thetaadvec_g);

		rollerlatbndCPU(nx, ny, ntheta, Param.eps, hh_g, rr_g);

	}

	//  Distribution of dissipation over directions and frequencies
	dissipationCPU(nx, ny, ntheta, dtheta, Param.eps, Param.dt, Param.g, Param.beta, wci_g, hh_g, ee_g, D_g, E_g, rr_g, c_g, cxsth_g, sxnth_g, uu_g, vv_g, DR_g, R_g);

	//Fix lateraL BND
	rollerlatbndCPU(nx, ny, ntheta, Param.eps, hh_g, ee_g);

	//  Compute mean wave direction
	meandirCPU(nx, ny, ntheta, Param.rho, Param.g, dtheta, ee_g, theta_g, thetamean_g, E_g, H_g);

	// Radiation stresses and forcing terms
	radstressCPU(nx, ny, ntheta, Param.dx, dtheta, ee_g, rr_g, cxsth_g, sxnth_g, cg_g, c_g, Sxx_g, Sxy_g, Syy_g);

	// Wave forces
	wavforceCPU(nx, ny, ntheta, Param.dx, dtheta, Sxx_g, Sxy_g, Syy_g, Fx_g, Fy_g, hh_g);

	//Lat Bnd
	twodimbndnoixCPU(nx, ny, Param.eps, hh_g, Fx_g);
	twodimbndCPU(nx, ny, Param.eps, hh_g, Fy_g);

	// CAlculate stokes velocity and breaker delay //Breaker delay removed because it is slow and kinda useless
	breakerdelayCPU(nx, ny, ntheta, dtheta, Param.g, Param.rho, Trep, Param.eps, urms_g, ust_g, H_g, E_g, c_g, k_g, hh_g, R_g);
	twodimbndCPU(nx, ny, Param.eps, hh_g, ust_g);

	// Adjust Offshore Bnd
	//CUDA_CHECK( cudaMemcpy(St_g, St, ny*ntheta*sizeof(DECNUM ), cudaMemcpyHostToDevice) );
	//offshorebndWav(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g)
	//offshorebndWav<<<gridDim, blockDim, 0>>>(nx,ny,ntheta,totaltime,Trep,St_g,sigm_g,ee_g);
	////CUT_CHECK_ERROR("Offshore Wave bnd execution failed\n");
	//CUDA_CHECK( cudaThreadSynchronize() );





}
