/*fftw_routines.cpp*/
/*
 * Copyright 2011, 2010 IPOL Image Processing On Line http://www.ipol.im/
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
 * @file fftw_routines.cpp
 * @brief Routines to use FFTW more simply.
 * @author Yohann Tendero <yohann.tendero@cmla.ens-cachan.fr>
 */

using namespace std;

#include <fftw3.h>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FFTW_FORWARD ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * @fn forward_fftw_simple(float* in,int width,int height,float* reOut,
 * float* imOut)
 * @brief Easy DFT FORWARD (call to the FFTW library)
 *
 * Arguments call: image, width, height,#code, Real_part_FFT,Im_part_FFT
 * Real_part_FFT,Im_part_FFT => contains the FFT of ''in''
 * @param *in : input image
 * @param width of *in
 * @param height of *in
 * @param *reOut : real-part of FFTW of *in
 * @param *imOut : imaginary-part of FFTW of *in
 */
void forward_fftw_simple(double* in,int width,int height,double* reOut,
                         double* imOut)
{
    ///init
    fftw_complex* spatial_repr;
    fftw_complex* frequency_repr;
    fftw_plan plan;
    spatial_repr= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*width*height);
    frequency_repr= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*width
                    *height);


    for (int i=0; i<width*height; i++)
    {
        spatial_repr[i][0] = in[i];
        spatial_repr[i][1] =  0.0f;
    }

    ///*Plan of the FFTW*/
    plan=fftw_plan_dft_2d(height, width, spatial_repr, frequency_repr,
                          FFTW_FORWARD, FFTW_ESTIMATE);

    ///*Execute the plan*/
    fftw_execute(plan);
    for(int i=0; i<width*height; i++)
    {
        reOut[i]=frequency_repr[i][0];
        imOut[i]=frequency_repr[i][1];
    }


    fftw_destroy_plan(plan);
    fftw_free(spatial_repr);
    fftw_free(frequency_repr);

}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FFTW_BACKWARD ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * @fn backward_fftw_simple(float* reIn,float* imIn,float* out,unsigned
 * int width, unsigned int height)
 * @brief  Easy DFT BACKWARD (call to the FFTW library)
 * Arguments call: Real_part_FFT,Im_part_FFT, image (output), width,height
 * The output is ''image''
 * @param *reIn : real-part of the input
 * @param *imIn : imaginary-part of the input
 * @param *out : output FFT^{-1}
 * @param widht of *reIn & *imIn
 * @param height of *reIn & *imIn
 */
void backward_fftw_simple(double* reIn,
                          double* imIn,
                          double* out,
                          unsigned int width,
                          unsigned int height)
{

    fftw_complex* spatial_repr;
    fftw_complex* frequency_repr;

    fftw_plan plan;

    spatial_repr=  (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*width*
                   height);
    frequency_repr=  (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*width*
                     height);



    for(unsigned i=0; i<width*height; i++)
    {
        frequency_repr[i][0]=reIn[i];
        frequency_repr[i][1]=imIn[i];
    }

    plan=fftw_plan_dft_2d(height, width, frequency_repr, spatial_repr,
                          FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    /// Dropping the complex part and renormalization (divide by width*height)
    for (unsigned i=0; i<width*height; i++)
    {
        out[i]=spatial_repr[i][0]/(width*height);
    }

    fftw_destroy_plan(plan);
    fftw_free(spatial_repr);
    fftw_free(frequency_repr);
}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FFTW_MULTIPLICATION /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
* @fn fftw_multiplication(float* reTabImage,float* imTabImage
* float* reTabFilter,float* imTabFilter,unsigned int width, unsigned int height)
* @brief  Implement the multiplication between two FT (convolution)
*
* Arguments call: Real_part_FFT1,Im_part_FFT1,Real_part_FFT2,Im_part_FFT2, width
* ,height
* The output are put in ''Real_part_FFT1'' & ''Im_part_FFT1''
* @param *reTabImage : real-part of the FFTW
* @param *imTabImage : imaginary-part of the FFTW
* @param *reTabFilter : real-part of the FFTW
* @param *imTabFilter : imaginary-part of the FFTW
* @param width of *reTabImage, *imTabImage, *reTabFilter, *imTabFilter
* @param height of *reTabImage, *imTabImage, *reTabFilter, *imTabFilter
*
*/
void fftw_multiplication(double* reTabImage,
                         double* imTabImage,
                         double* reTabFilter,
                         double* imTabFilter,
                         unsigned int width,
                         unsigned int height)
{

    float a,b,c,d;
    unsigned int i;

    for (i=0; i< width * height; i++)
    {
        a = reTabImage[i];
        b = imTabImage[i];
        c = reTabFilter[i];
        d = imTabFilter[i];

        ///*Computing Z1.Z2*/
        reTabImage[i] = a*c - b*d;
        imTabImage[i] = b*c + a*d;


    }

}







////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// FFTW_DIVISION ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * @fn fftw_division(float* reTabImage, float* imTabImage,float* reTabFilter,
 * float* imTabFilter,unsigned int width,unsigned int height)
* @brief  Implement the division between two FT (deconvolution)
* Arguments call: Real_part_FFT1,Im_part_FFT1,Real_part_FFT2,Im_part_FFT2, width
* ,height
* The output are put in ''Real_part_FFT1'' & ''Im_part_FFT1''
* @param *reTabImage : real-part of the FFTW
* @param *imTabImage : imaginary-part of the FFTW
* @param *reTabFilter : real-part of the FFTW
* @param *imTabFilter : imaginary-part of the FFTW
* @param width of *reTabImage, *imTabImage, *reTabFilter, *imTabFilter
* @param height of *reTabImage, *imTabImage, *reTabFilter, *imTabFilter
*
*/
void fftw_division(double* reTabImage,
                   double* imTabImage,
                   double* reTabFilter,
                   double* imTabFilter,
                   unsigned int width,
                   unsigned int height)
{

    float a,b,c,d,temp;
    unsigned int i;

    for (i=0; i< width * height; i++)
    {
        a = reTabImage[i];
        b = imTabImage[i];
        c = reTabFilter[i];
        d = imTabFilter[i];

        ///* Computing Z1/Z2 */
        temp=c*c+d*d;
        reTabImage[i] = (a*c+b*d)/temp;
        imTabImage[i] = (b*c-a*d)/temp;

    }

}


