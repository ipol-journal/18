/* standard_routines.h */
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


void image_difference(float *groundtruth,float *restored, int width,
                      int height,int channel_number, float *difference);
//Compute thdifference.


void RMSE(float *difference, int width, int height,int channel_number,
          int code_length, int flag_rmse_ci);
//Compute the Root Mean Squarred Error (RMSE)

void dynamic_renormalization(float *image, int width, int height,
                             int channel_number);
//Affine contrast change => [0,255]

float abs_hat_alpha(const float* code, int code_length, float xi, float deltat);
//Compute the MODULUS of the Fourier transform of the
//flutter shutter function with code at \xi
