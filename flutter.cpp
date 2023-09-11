/*flutter.cpp*/
/*
 * Copyright 2012, 2010 IPOL Image Processing On Line http://www.ipol.im/
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
 * @file flutter.cpp
 * @brief Routines for Flutter-Shutter simulation.
 * This file contains functions that are usefull for both
 * kinds of flutter shutter camera (numerical or analog)
 * The kernel definition, the deconvolution, and a renormalization function.
 * @author Yohann Tendero <yohann.tendero@cmla.ens-cachan.fr>
 */





#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "codes_flutter.cpp"
#include "flutter.h"
#include "fftw_routines.h"
#include "borders.h"
#include "mt19937ar.h"
#ifndef M_PI
/**
 * M_PI is a POSIX definition for Pi
 */
#define M_PI 3.14159265358979323846
#endif

#define ABS(x)    (((x) > 0) ? (x) : (-(x))) //absolute value
#define eps   0.000001//epsilon definition
#define POS(x) (((x) > 0) ? (x) : (0))



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FLUTTER_SHUTTER KERNEL COMPUTATION //////////////////
////////////////////////////////////////////////////////////////////////////////
/**
* @fn void flutter_kernel(float *kernel_real, float *kernel_imag,
* int code_number, int width, int height, float velocity, float deltat)
* @brief Given a code compute the Flutter-Shutter kernel (Fourier domain)
* Arguements: kernel_real,float kernel_imag,
*  number of the code (codes are stored in codes_flutter.cpp), width,
*  height, velocity
* The output are ''kernel_real'' & ''float kernel_imag''
* @param *kernel_real : image containing the real-part of the kernel
* (Fourier domain);
* @param *kernel_imag : image containing the imaginary-part of the kernel
*  (Fourier domain);
* @param (int) code_number : index of the code in the list;
* @param (int) width : width of the image used for simulation;
* @param (int) height : height of the image used for simulation;
* @param (float) velocity : relative velocity between the camera and
* the landscape;
* @param (float) deltat : exposure time.
*/
void flutter_kernel(double *kernel_real, double *kernel_imag, int code_number,
                    int width, int height, float velocity, float deltat)
{

///The kernel is
///k(\xi)=2 deltat \frac{sin(\frac{\xi v deltat} 2 }{\xi v deltat}
/// \sum_{k=0}^{N-1}\alpha_k exp(-i \xi v deltat \frac{2k+1}{2}
/// \xi \in [-\pi,\pi], (constant among lines).



    for (int j=0; j<height; j++)
    {
        float xi1=-width/2;
        for (int i=floor((width+1)/2); i<width; i++)
        {

            float xi=2*M_PI*velocity*deltat/width*xi1;
            xi1=xi1+1;

            for (int k=0; k<code_length; k++)
            {

                kernel_imag[i+width*j] =kernel_imag[i+width*j]+
                                        code[code_number][k]*
                                        sin(-xi*(k+0.5));
                kernel_real[i+width*j] =kernel_real[i+width*j]+
                                        code[code_number][k]*
                                        cos(-xi*(k+0.5));


            }
            if (ABS(xi)>eps)  //avoiding 0/0;else =deltat
            {
                kernel_imag[i+width*j] =kernel_imag[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
                kernel_real[i+width*j] =kernel_real[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
            }
            else
            {
                kernel_imag[i+width*j]=deltat*kernel_imag[i+width*j];
                kernel_real[i+width*j]=deltat*kernel_real[i+width*j];
            }

        }

        for (int i=0; i<floor((width+1)/2); i++)
        {

            float xi=2*M_PI*velocity*deltat/width*xi1;
            xi1=xi1+1;

            for (int k=0; k<code_length; k++)
            {

                kernel_imag[i+width*j] =kernel_imag[i+width*j]+
                                        code[code_number][k]*
                                        sin(-xi*(k+0.5));
                kernel_real[i+width*j] =kernel_real[i+width*j]+
                                        code[code_number][k]*
                                        cos(-xi*(k+0.5));


            }
            if (ABS(xi)>eps)  //avoiding 0/0;else =1
            {
                kernel_imag[i+width*j] =kernel_imag[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
                kernel_real[i+width*j] =kernel_real[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
            }
            else
            {
                kernel_imag[i+width*j]=deltat*kernel_imag[i+width*j];
                kernel_real[i+width*j]=deltat*kernel_real[i+width*j];
            }

        }



    }


}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FLUTTER_SHUTTER DECONVOLUTION ///////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * @fn flutter_restore(float *acquired,int code_number,int width, int height,
 * float velocity, float *restored, int border_size)
 * @brief Given the observed (noisy) image, and the device parameters it
 * performs the restoration - all parameters ar assumed to be known-
 * @param *acquired : acquired (observed) image
 * @param (int) code_number : index of the code in the list
 * @param (int) width of *Image
 * @param (int) height of *Image
 * @param (float) velocity  : relative velocity between the camera
 * and the landscape (counted in pixels per unit time $\Deltat$
 * @param *restored : output (restored image)
 * @param (int) border_size : the *acquired image is mirror symetrized of
 *  "border_size" columns before the deconvolution to reduce borders effects;
 * @param (float) deltat : temporal sampling.
 */

void flutter_restore(double *acquired,int code_number, int width, int height,
                     float velocity, double *restored,
                     int border_size, float deltat)
{
    ///Remark : not computing the Kernel one and for all
    ///because the size change (crop).


    int W=width+border_size; //new  width of image afer mirror
    double* Imsym = new double[W*height]; //new image allocation



    /// Classic Mirror symmetryzation
    borders(acquired,Imsym,width,height);



    ///Allocating DFT of input for deconvolution and compute the DFT of image
    double* reOut = new double[W*height];
    double* imOut = new double[W*height];

    forward_fftw_simple(Imsym, W, height, reOut, imOut);

    delete [] Imsym;
    ///Allocating DFT of input for deconvolution and compute
    ///the DFT of the kernel
    double* kernel_real = new double[W*height];
    double* kernel_imag = new double[W*height];

    for ( int i=0; i<W*height; i++)
    {
        kernel_real[i]=0.0f;
        kernel_imag[i]=0.0f;
    }

    flutter_kernel(kernel_real, kernel_imag, code_number, W, height,
                   velocity, deltat);


    ///DECONVOLUTION
    fftw_division(reOut, imOut, kernel_real, kernel_imag, W, height);
    delete [] kernel_real;
    delete [] kernel_imag;


    backward_fftw_simple(reOut, imOut, restored, W, height);
    delete [] reOut;
    delete [] imOut;
}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FLUTTER_SHUTTER INTEGRAL COMPUTATION ////////////////
//////////////////////////  (for normalization) ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * @fn integral_code(int code_number)
 * @brief Given a code compute it's integral (the integral of the flutter
 * shutter function associated).
 * Arguments: code_number (codes are stored in codes_flutter.cpp)
 * The output is a float : \int
* @param (int) code_number : index of the code in the list
* @return I the integral of the code.
 */
float integral_code(int code_number, float deltat)
{
    float I=0;
    for (int k=0; k<code_length; k++)
    {
        I=I+code[code_number][k];

    }
    return I*deltat;
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FLUTTER_SHUTTER_k  (NUMERICAL FLUTTER KERNEL) ///////
////////////////////////////////////////////////////////////////////////////////
/**
* @fn void flutter_kernel_k(float *kernel_real, float *kernel_imag, int k,
* float ak, int width, int height, float velocity, float deltat)
* @brief Given k computes the Flutter-Shutter kernel (Fourier domain)
* associated to the code (0,..., \alpha_k,0, ...0)
* used for numerical flutter shutter simulation;
* whose elements are zeroes except the k-th element equal to a_k
* Arguements: kernel_real,float kernel_imag,
*  number of the code (codes are stored in codes_flutter.cpp), width,
*  height, velocity
* The output are ''kernel_real'' & ''float kernel_imag''
* @param *kernel_real : image containing the real-part of the kernel
* (Fourier domain);
* @param *kernel_imag : image containing the imaginary-part of the kernel
*  (Fourier domain);
* @param (int) k : the index of the only non zero element
* @param (float) ak : value of the k-th element of the code
* @param (int) width : width of the image used for simulation;
* @param (int) height : height of the image used for simulation;
* @param (float) velocity : relative velocity between the camera and
* the landscape;
* @param (float) deltat : temporal sampling.
*/
void flutter_kernel_k(double *kernel_real, double *kernel_imag, int k,
                      int width, int height, float velocity, float deltat)
{

///The kth kernel is
///k(\xi)=2 deltat \frac{sin(\frac{\xi v deltat} 2 }{\xi v deltat}
/// exp(-i \xi v deltat \frac{2k+1}{2}
/// \xi \in [-\pi,\pi]


    for (int j=0; j<height; j++)
    {
        float xi1=-width/2;
        for (int i=floor((width+1)/2); i<width; i++)
        {

            float xi=2*M_PI*velocity*deltat/width*xi1;
            xi1=xi1+1;

            kernel_imag[i+width*j] =sin(-xi*(k+0.5));
            kernel_real[i+width*j] =cos(-xi*(k+0.5));



            if (ABS(xi)>eps) //avoiding 0/0;else =1
            {
                kernel_imag[i+width*j] =kernel_imag[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
                kernel_real[i+width*j] =kernel_real[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
            }
            else
            {
                kernel_imag[i+width*j]=deltat*kernel_imag[i+width*j];
                kernel_real[i+width*j]=deltat*kernel_real[i+width*j];
            }
        }
        for (int i=0; i<floor((width+1)/2); i++)
        {

            float xi=2*M_PI*velocity*deltat/width*xi1;
            xi1=xi1+1;

            kernel_imag[i+width*j] =sin(-xi*(k+0.5));
            kernel_real[i+width*j] =cos(-xi*(k+0.5));



            if (ABS(xi)>eps) //avoiding 0/0;else =1
            {
                kernel_imag[i+width*j] =kernel_imag[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
                kernel_real[i+width*j] =kernel_real[i+width*j]*
                                        ((2.0*deltat)*sin(xi/2.0))/xi;
            }
            else
            {
                kernel_imag[i+width*j]=deltat*kernel_imag[i+width*j];
                kernel_real[i+width*j]=deltat*kernel_real[i+width*j];
            }
        }

    }



}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FLUTTER_SHUTTER NOISY SIMULATION/////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * @fn void  add_poisson_noise(float *image, int width, int height,
                        float deltat, float *acquired, float snr)
 *
 * @brief Add Poisson noise to the ideal noise-less acquired image
 * @param *Image : input and output (void) image.
 * @param width of *Image
 * @param height of *Image
 * @param velocity  : relative speed between the camera and the landscape
 * @param *acquired : output ie given a code, a velecity, a SNR-level (snr) a
 *  Flutter-Shutter type and a landscape (*Image) simulate the acquired image by
 * the device
 * @param snr : set STD-DEV of the noise
 */
void  add_poisson_noise(double *image, int width, int height, float snr)
{

    int em;
    double lambda;
    bool rejected;
    double t;
    for (int i=0; i<width*height; i++)
    {

        lambda=POS(image[i]*snr*snr/100);
        //SNR FOR LEVEL 100 =>Poisson intensity to simulate
        // due to small errors mage[i]
        // can be \approx -10^{-6}
        // this avoids any trouble when taking the square root/
        ///contains the intensity lambda=\alpha \ast u (k)

        if ((image[i]<5000/(snr*snr)))
            ///Small intensity (after renormalization)

        {
            lambda=exp(-lambda);
            em=-1;
            t=1;
            rejected=true;
            while (rejected)
            {
                em=em+1;
                t=t*mt_genrand_res53();
                if (t<=lambda)
                {
                    image[i]=100*em/(snr*snr);//rescaling
                    rejected=false;
                }
            }
        }
        else ///Bigger intensity : approximation by a Gaussian
        {
            image[i]=round(lambda +
                           ((double)sqrt(-2.0 * log(mt_genrand_res53()))
                            * (double)cos(2.0 * M_PI * mt_genrand_res53())*
                            (double)(sqrt(lambda))))*100/(snr*snr);

        }

    }
}




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////// NUMERICAL FLUTTER_SHUTTER  SIMULATOR ////////////
////////////////////////////////////////////////////////////////////////////
/**
* void flutter_numerical(float *Image,int code_number,int width, int height,
* float velocity, float deltat, float *acquired, float snr)
* @brief Given a code it simulates a numerical flutter shutter
* @param *Image : image used for simulation
* @param code_number : : index of the code in the list
* @param width of *Image
* @param height of *Image
* @param velocity  : relative velocity between the camera and the landscape
* @param deltat : temporal sampling.
* @param *acquired : impage observed
* @param snr : noise level.
*/
void flutter_numerical(double *acquired, int code_number, int width, int height,
                       float velocity, float deltat, float snr)
{
    //A numerical flutter shutter is a sum of code_length analog flutter shutter
    // with codes : (0, ... , 0, ak , 0, ... , 0)
    // zeros except on k-th position : ak (=code[code_number][k])

    ///Allocating DFT of input
    double* reOut = new double[width*height];
    double* imOut = new double[width*height];

    for ( int i=0; i<width*height; i++)
    {
        reOut[i]=0.0f;
        imOut[i]=0.0f;
    }


    ///FFT FORWARD IMAGE
    forward_fftw_simple(acquired, width, height, reOut, imOut);

    /// COMPUTING THE FT OF THE CODE
    /// allocating, init

    for ( int i=0; i<width*height; i++)
    {
        acquired[i]=0.0f;
    }


    double** reOut_temp = new double*[code_length]; // to store the k-th
    double** imOut_temp = new double*[code_length]; // elementaryfourier transf.
    double** im_temp = new double*[code_length];//to store the k th observed
    double** kernel_real = new double*[code_length];// to store the k-th kernel
    double** kernel_imag = new double*[code_length];

    for (int k=0; k<code_length; k++)
    {
        reOut_temp[k] = new double[width*height];
        imOut_temp[k] = new double[width*height];
        im_temp[k] = new double[width*height];
        kernel_real[k] = new double[width*height];
        kernel_imag[k] = new double[width*height];
    }




    for ( int i=0; i<width*height; i++)
    {
        for (int k=0; k<code_length; k++)
        {
            reOut_temp[k][i]=reOut[i];
            imOut_temp[k][i]=imOut[i];
            im_temp[k][i]=0.0f;
            kernel_real[k][i]=0.0f;
            kernel_real[k][i]=0.0f;
        }
    }

    delete [] reOut;
    delete [] imOut;



    for (int k=0; k<code_length; k++)
    {

//filling the k th kernel
        flutter_kernel_k(kernel_real[k], kernel_imag[k], k,
                         width, height, velocity, deltat);
//convolution
        fftw_multiplication(reOut_temp[k], imOut_temp[k], kernel_real[k],
                            kernel_imag[k],
                            width , height);

        backward_fftw_simple(reOut_temp[k], imOut_temp[k], im_temp[k], width ,
                             height);
        //here im_temp[k] contains the k-th observed withtout noise
        //Noise:
        add_poisson_noise(im_temp[k], width, height, snr);
        //here im_temp contains the k-th observation
        //(and without flutter)
    }



    for (int k=0; k<code_length; k++)
    {

        for ( int i=0; i<width*height; i++)
        {
            //adding to get the observed and flutter effect (the code)
            acquired[i]=acquired[i]+code[code_number][k]*im_temp[k][i];

        }
        delete [] reOut_temp[k];
        delete [] imOut_temp[k];
        delete [] kernel_real[k];
        delete [] kernel_imag[k];
        delete [] im_temp[k];
    }

    delete [] reOut_temp;
    delete [] imOut_temp;
    delete [] kernel_real;
    delete [] kernel_imag;
    delete [] im_temp;


}




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////// FLUTTER_SHUTTER ANALOG //////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
* @fn flutter_anlog(float *acquired, int code_number, int width,
                             int height, float velocity, float deltat,
                             float snr)
* @brief Given a code  it simulates the analog flutter shutter observed.
* Arguements: Image, code_number, width,  height,  velocity, noise_less
* The output is ''noise_less''
* @param *acquired : image used for simulation
* @param (int) code_number : : index of the code in the list
* @param (int) width of *Image
* @param (int) height of *Image
* @param (float) velocity  : relative velocity between the camera and the
* landscape
* @param (float) deltat : temporal sampling
* @param (flaot) snr : noise level.
*/
void flutter_analog(double *acquired, int code_number, int width,
                    int height, float velocity, float deltat,
                    float snr)
{

    ///Allocating DFT of input
    double* reOut = new double[width*height];
    double* imOut = new double[width*height];
    for ( int i=0; i<width*height; i++)
    {
        reOut[i]=0.0f;
        imOut[i]=0.0f;
    }

    ///FFT FORWARD IMAGE
    forward_fftw_simple(acquired, width, height, reOut, imOut);

    /// COMPUTING THE FT OF THE CODE
    /// allocating, init
    double* kernel_real = new double[width*height];
    double* kernel_imag = new double[width*height];
    for ( int i=0; i<width*height; i++)
    {
        kernel_real[i]=0.0f;
        kernel_imag[i]=0.0f;

    }
    /// actual computation
    flutter_kernel(kernel_real, kernel_imag, code_number, width, height,
                   velocity, deltat);



    ///CONVOLUTION
    fftw_multiplication(reOut, imOut, kernel_real, kernel_imag, width , height);
    backward_fftw_simple(reOut, imOut, acquired, width , height);
    //Here noise_less contains the motion blurred landscape (no noise)
    delete [] reOut;
    delete [] imOut;
    delete [] kernel_real;
    delete [] kernel_imag;
    ///NOISE

    add_poisson_noise(acquired, width, height, snr);
}

