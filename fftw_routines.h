/* fftw_routines.h */
/*
* Copyright 2012 IPOL Image Processing On Line http://www.ipol.im/
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


void forward_fftw_simple(double* in,int width,int height, double* reOut,
                         double* imOut);
//Compute the DFT


void backward_fftw_simple(double* reIn,
                          double* imIn,
                          double* out,
                          unsigned int width,
                          unsigned int height);
//Compute the DFT inverse


void fftw_multiplication(double* reTabImage,
                         double* imTabImage,
                         double* reTabFilter,
                         double* imTabFilter,
                         unsigned int width,
                         unsigned int height);
//Used for easy convolution





void fftw_division(double* reTabImage,
                   double* imTabImage,
                   double* reTabFilter,
                   double* imTabFilter,
                   unsigned int width,
                   unsigned int height);

//Used for easy deconvolution


