/*standard_routines.cpp*/
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
/**
 * @file standard_routines.cpp
 * @brief Routines for computing difference, RMSE, PNSR,
 * dynamic renormalization, Flutter Shutter Fourier transform
 *
 * @author Yohann Tendero <yohann.tendero@cmla.ens-cachan.fr>
 */


#include <math.h>
#include <stdio.h>
#define ABS(x)    (((x) > 0) ? (x) : (-(x))) //absolute value


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// DIFFERENCE BETWEEN 2 IMAGES /////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
* @fn image_difference(float *groundtruth,float *restored, int width, int height
* ,int channel_number, float *difference)
* @brief Compute the difference between two images
* Arguements:  groundtruth, restored,  width,  height,  difference
* The output is ''difference''
* @param *groundtruth : image 1
* @param *restored : image 2
* @param width of image1 & 2
* @param height of image1 & 2
* @param channel_number number of channels of image1 & 2
* @param *difference : output
*/
void image_difference(float *groundtruth,float *restored, int width, int height,
                      int channel_number, float *difference)
{

    for (int i=0; i<width*height*channel_number; i++)
    {
        difference[i]=groundtruth[i]-restored[i];

    }

}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// RMSE ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
* @fn RMSE(float *difference, int width, int height,int channel_number,
int border,int flag_rmse_ci)
* @brief Compute the RMSE given a difference image and print it (printf).
* Arguements:  difference, width, height, code_length, flag for printf
* The output is a printf
* @param *difference : image containing the difference between the groundtruth
* and the actual values.
* @param width of difference
* @param height of difference
* @param channel_number : number of channels of difference
* @param border : (type : integer) exclude the firsts and lasts "border" columns
* in order to avoid mesurements of berders effects
* @param flag_rmse_ci : only change the printf
*/
void RMSE(float *difference, int width, int height,int channel_number,int border
          ,int flag_rmse_ci)
{
    float rmse=0;
    int num_pixels=0;
    for (int k=0; k<channel_number; k++)
        for (int j=0; j<height; j++)
        {
            for (int i=border; i<width-border; i++)
            {
                rmse=rmse+difference[i+width*j+k*width*height]*
                     difference[i+width*j+k*width*height];
                num_pixels++;
            }
        }
    rmse=sqrt(rmse/(num_pixels));
    float psnr = 10.0f * log10f(255.0f * 255.0f / (rmse * rmse));
    //std::cerr << "RMSE=" <<  std::endl;
    //std::cerr << rmse <<  std::endl;
    //std::cerr << "PSRN=" <<  std::endl;
    //std::cerr << psnr <<  std::endl;
    if (flag_rmse_ci==0)
    {
        printf("RMSE: %f\n", rmse);
        printf("PSNR: %f\n", psnr);
    }
    else
    {
        printf("Contrast Invariant RMSE: %f\n", rmse);
        printf("Contrast Invariant PSNR: %f\n", psnr);
    }


}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// Dynamic renormalization to [0,255] //////////////////
////////////////////////////////////////////////////////////////////////////////
/**
* @fn dynamic_renormalization(float *image, int width, int height,
* int channel_number)
* @brief Change the dynamic to [0,255]
* Arguements:  image, width, height
* The output is ''image'
* @param *image : image
* @param width of image
* @param height of image
* @param channel_number number of channels of image
*/
void dynamic_renormalization(float *image, int width, int height,
                             int channel_number)
{

    ///Computing the mininum and max of image
    float min=image[0];
    float max=image[0];
    for (int i=1; i<width*height*channel_number; i++)
    {
        if (image[i]<min) min=image[i];
        if (image[i]>max) max=image[i];
    }

    ///Actually changing image-values
    for (int i=1; i<width*height*channel_number; i++)
    {
        image[i]=(255*(image[i]-min)/(max-min));
    }
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// MODULUS OF THE CODE COMPUTATION///////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * @fn abs_hat_alpha(float* code, int code_length, float xi)
 * @brief Given a code compute the modulus of the Fourier transform of
 * the Flutter-Shutter function defined by :
 * \alpha(t)=\sum_{k=0}^{L-1} \alpha_k \mathbb{1}_{[k\Deltat,(k+1)\Deltat[}(t)
 * @param float* code : array of floats where code[k]
 * contains the \alpha_k;
 * @param int code_length : length of the code;
 * @param float xi : frequency;
 * @param (positive float) deltat : the temporal sampling of the
 * flutter shutter function;
 * @return |\hat \alpha (\xi)|.
 */
float abs_hat_alpha(const float* code, int code_length, float xi, float deltat)
{
    ///Gives back the modulus of the Fourier transform
    //initialization
    float re_abs_hat_alpha=0; //real part
    float im_abs_hat_alpha=0; //imaginary part
    //Main loop
    for (int k=0; k<code_length; k++)
    {
        im_abs_hat_alpha =im_abs_hat_alpha+code[k]*sin(-xi*deltat*(k+0.5));
        re_abs_hat_alpha =re_abs_hat_alpha+code[k]*cos(-xi*deltat*(k+0.5));
    }

    if (ABS(xi)>0) //avoiding 0/0
    {
        im_abs_hat_alpha =im_abs_hat_alpha*sin(xi*deltat/2)
                          /(xi*deltat/2);//ELSE =1;
        re_abs_hat_alpha =re_abs_hat_alpha*sin(xi*deltat/2)
                          /(xi*deltat/2); //ELSE =1;
    }

    return(deltat*pow(im_abs_hat_alpha*im_abs_hat_alpha+
                      re_abs_hat_alpha*re_abs_hat_alpha,0.5));
    ///RETURNS the modulus at xi of the flutter shutter function given its code:
    /// \Deltat \frac{sin(\frac{\xi }{2})}{\frac \xi 2}
    /// \sum_{k=0}^{L-1} \alpha_k e^{i\xi  \frac{2k+1}{2}}
}
