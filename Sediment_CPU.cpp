#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>

#define pi 3.14159265
using DECNUM = float
;



extern "C" void updatezomCPU(int nx, int ny,DECNUM cf,DECNUM cf2,DECNUM fw,DECNUM fw2,DECNUM * structdepth, DECNUM * &cfm,DECNUM * &fwm)
{
	for (int ix=0; ix<nx; ix++)
		{
			for (int iy=0; iy<ny; iy++)
				{
					int i=ix+iy*nx;
	
	
	
	
	
					if (structdepth[i]<0.05f)
					{
						cfm[i]=cf2;
						fwm[i]=fw2;
					}
					else
					{
						cfm[i]=cf;
						fwm[i]=fw;
					}
			}
	}
	
	
	
}

