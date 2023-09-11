/* flutter.h */
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
extern  const int num_code; //number of code
extern const int code_length; //Number of coefficients of the code
//extern  const float code[(const)code_length][(const)code_length];
void flutter_kernel(double [], double [], int,
                    int, int, float, float);

void flutter_restore(double [], int,int, int,
                     float, double [],
                     int, float);

float integral_code(int, float);

void flutter_kernel_k(double [], double [], int ,
                      int, int, float, float);


void  add_poisson_noise(float [], int, int, float);

void flutter_numerical(double [], int, int, int,
                       float, float, float);

void flutter_analog(double [], int, int,
                    int, float, float,
                    float);
