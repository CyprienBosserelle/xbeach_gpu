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
double interptime(double next, double prev, double timenext, double time)
{
	return prev + (time) / (timenext)*(next - prev);
}


double interp1D(int nx, double *x, double *y, double xx)
{
	double yy;
	double prevx = x[0];
	double nextx = x[1];
	double prevy = y[0];
	double nexty = y[1];
	int indx = 0;
	int indxp;

	double diffx = 0;



	for (int i = 0; i < nx; i++)
	{
		diffx = xx - x[i];
		if (diffx <= 0.0)
		{
			indx = (int)max(i*1.0 - 1, 0.0);
			break;
		}
	}
	indxp = (int)min(indx*1.0 + 1, nx*1.0 - 1);
	prevx = x[indx];
	nextx = x[indxp];
	prevy = y[indx];
	nexty = y[indxp];


	yy = prevy + (xx - prevx) / (nextx - prevx)*(nexty - prevy);
	return yy;
}
double Interp2(int nx, int ny, double *x, double *y, double *z, double xx, double yy)
//double BilinearInterpolation(double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2, double x, double y)
{
	double x2x1, y2y1, x2x, y2y, yy1, xx1;
	double q11, q12, q21, q22, x1, x2, y1, y2;


	// find x1 and x2
	int indx = 0;
	int indxp;

	double diffx;



	for (int i = 0; i < nx; i++)
	{
		diffx = xx - x[i];
		if (diffx <= 0.0)
		{
			indx = (int)max(i*1.0 - 1, 0.0);
			break;
		}
	}

	x1 = x[indx];
	indxp = (int)min(indx*1.0 + 1, nx*1.0 - 1);
	x2 = x[indxp];

	int indy = 0;
	int indyp;

	double diffy;



	for (int i = 0; i < ny; i++)
	{
		diffy = yy - y[i];
		if (diffy <= 0.0)
		{
			indy = (int)max(i*1.0 - 1, 0.0);
			break;
		}
	}

	y1 = y[indy];
	indyp = (int)min(indy*1.0 + 1, ny*1.0 - 1);
	y2 = y[indyp];

	q11 = z[indx + indy*nx];
	q12 = z[indx + indyp*nx];
	q21 = z[indxp + indy*nx];
	q22 = z[indxp + indyp*nx];

	x2x1 = x2 - x1;
	y2y1 = y2 - y1;
	x2x = x2 - xx;
	y2y = y2 - yy;
	yy1 = yy - y1;
	xx1 = xx - x1;
	return 1.0 / (x2x1 * y2y1) * (
		q11 * x2x * y2y +
		q21 * xx1 * y2y +
		q12 * x2x * yy1 +
		q22 * xx1 * yy1
		);
}
// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive

void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * pi * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
/*void fft(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}
*/
// inverse fft (in-place)
void ifft(CArray& x)
{
	// conjugate the complex numbers
	for (int j = 0; j < x.size(); j++)
	{
		x[j] = std::conj(x[j]);
	}

	// forward fft
	fft(x);

	// conjugate the complex numbers again
	for (int j = 0; j < x.size(); j++)
	{
		x[j] = std::conj(x[j]);
	}

	// scale the numbers
	x /= x.size();
}

void hilbert(CArray& xi)
{
	int n,m;
	double p;

	m = xi.size();
	n =ceil( log(m*1.0) / log(2.0));
	n = pow(2, n);

	
	CArray x(n);
	x = 0.0;

	for (int i = 0; i < xi.size(); i++)
	{
		x[i] = std::real(xi[i]);
	}
	
	fft(x);

	for (int j = 0; j < x.size(); j++)
	{
		x[j] = x[j] * sqrt(n);
	}

	std::valarray<double> h(0.0, n);

	h[0] = 1.0;
	h[std::slice(1, n / 2 - 1, 1)] = 2.0;
	h[n / 2] = 1.0;
	//Below is a redundant operation since h should be initialised with zeros
	//h[std::slice(n / 2 + 1, n - (n / 2 + 1), 1)] = 0.0;

	for (int j = 0; j < x.size(); j++)
	{
		x[j] = x[j]*h[j];
	}

	ifft(x);

	//scale?
	for (int j = 0; j < x.size(); j++)
	{
		x[j] = x[j] / sqrt(n);
	}

	for (int i = 0; i < xi.size(); i++)
	{
		xi[i] = x[i];
	}


}

void flipiv(CArray &x)
{
	//flip of a vall array of complex upside down

	CArray xf = x;
	int m = x.size();
	for (int j = 0; j < m; j++)
	{
		x[j] = xf[m - (j + 1)];
	}
}
