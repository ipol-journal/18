/*demo_flutter.cpp*/
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
 * @mainpage Flutter-Shutter Simulator
 * @brief It simulates a Flutter Shutter camera :
 * given an image it simulates the observed image (motion blurred and noisy)
 * form the Poisson (photon) noise then deconvolve it.
 * Notice that motion blurs are made invertible by the Flutter Shutter apparatus
 * README.txt:
 * @verbinclude README.txt
 */

/**
 * @file demo_flutter.cpp
 * @brief Flutter-Shutter camera simulator main file
 * The input is (grayscale or color) png image, a "code" among the list,
 * a noise-level, a Flutter Shutter type (numerical or analog), a velocity $v$.
 * The output are 1) acquired image using a Flutter-Shutter acquisition strategy
 * 2) The reconstructed image.
 * 3) Groundtruth image
 * 4) The difference between the ideal reconstruction and the actual
 * reconstruction with dynamic on [0,255] by an affine contrast change.
 * 5) Txt file that contains the code coefficients
 * 6) Txt file that contains the Fourier transform of the flutter shutter
 * function. (used for invertibility check)
 * 7) The RMSEs (standard output; printf).
 *
 * IN/OUTputs are 8bits png (grayscale or color)
 *
 *******************************************************************************
 *******************************************************************************
 ********************** Sketch of the simulator ********************************
 *******************************************************************************
 *******************************************************************************
 *
 * Assume a code of length $L$, (w.l.o.g)
 * Given a noise level and a velocity $v$
 *
 ** I) Simulation of the oberved image
 *
 *    1) Analog Flutter Shutter case :
 *       a) Simulate noise less image (convolution against the Flutter Shutter
 *       function of the code at velocity $v$);
 *       b) Simulate the Poisson R.Vs
 *       (given the noise level);
 *       c) Crop to avoid periodization artifacts introduced by the DFT
 *       to get the observed: L/2 on each sides
 *
 *    2) Numerical Flutter Shutter case :
 *       a) Simulate the $N$ images of an Analog Flutter Shutter using
 *       the $N$ codes : (1,0,...,0), (0,1,0,...,0), (0,0,1,0,...,0)
 *       to (0,...,0,1) at the choosen velocity and noise level;
 *       b) Multiply the $k-th$ image obtained above by the coefficient $a_k$
 *       of the code then sum to get the the observed;
 *
 *
 ** II) Deconvolution of the observed image
 *    1) Apply a mirror symmetry in the direction the blur;
 *    2) Deconvolve (in Fourier domain);
 *    3) Crop to remove the symmetry; crop for borders effects L on each sides
 *    4) Comute the RMSE (omitting 3L firsts and lasts columns),
 *       write images, etc.
 *
 * Color image : each channel is processed independatly
 *
 * @author Yohann Tendero <yohann.tendero@cmla.ens-cachan.fr>
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "mt19937ar.h"
#include <math.h>
#ifndef M_PI
/**
 * M_PI is a POSIX definition for Pi
 */
#define M_PI 3.14159265358979323846
#endif
#include "io_png/io_png.h"
#include "flutter.h"
#include "standard_routines.h"
#include "midway.h"
#include "borders.h"
#include "codes_flutter.cpp"
using namespace std;


/**
 * @brief main function call
 *	arguemnts are "Image.png" "Flag_{code}" "SNR" "Flag_{NUM | ANALOG}"
 * "Velocity" "Deltat" "Imacquired.png" "ImReconstructed.png"
 * "Groundtruth.png" "Imdiff.png" "code_file_name.txt" "FT_code_file_name.txt"
 *
 * @param "Input.png" : file name of the input image (PNG file)
 * Image to use for the simulation (grayscale or color) PNG.
 * @param (int) "Flag_{code}" indicates the code to use.
 * Codes for the Flutter-Shutter Simulator
 * Code : 0 for Snapshot, 1 for Raskar-code, 2 for rand-code,
 * 3 for constant to 1 code, 4 for Sinc-code: sinc code is the best code given v
 * 5 for motion-invariant photography code;
 * @param (positive float) "SNR" (Signal-to-Noise Ratio)
 * for the level 100 : 100 (default)"
 * @param (int) " Flag_{NUM | ANALOG}" is "0" or "1"
 * Flutter-shutter type "0 for Analog Flutter OR 1 for Numerical Flutter"
 * @param (float) "Velocity" v is a number"
 * it is the relative velocity between the camera and the landscape
 * counted in pixels per time unit $\Deltat$
 * @param (positive float) \Deltat : the temporal sampling of the
 * flutter shutter function
 * @param "Imacquired.png"  : file name for observed image (PNG file)
 * @param "ImReconstructed.png ": file name for recovered from the observed
 * image (PNG file)
 * @param "Groundtruth.png" : file name for the ideal recovered image (PNG file)
 * @param "Imdiff.png" : image difference between the recovered and the
 * ideal recovered (PNG file)
 * @param "code_file_name.txt" : file name for file that contains the code
 * coefficients  $a_k$ (txt file), used for plot generation later on
 * @param "FT_code_file_name.txt" : file name for file that contains the Fourier
 * transform of the flutter Shutter function coming from the code, used for
 * plot generation later on
 * @return "Output"
 * "Imacquired: acquired image  in PNG. "
 * "ImReconstructed: reconstructed  in PNG. "
 * "Groundtruth (ideal reconstructed)   in PNG. "
 * "Imdiff: difference (groundtruth and actual reconstruction) image  in PNG. "
 * "code_file_name.txt" : txt files that contains the code coefficients  $a_k$
 * "FT_code_file_name.txt" : txt file that contains the Fourier
 * transform of the flutter Shutter function coming from the code
 * Standard output (printf) : RMSE and RMSE contrast invariant
 * (RMSE contrast invariant is the RMSE after contrast normalization)
 *
 * Color images : each channel is processed independantly
 * @verbinclude README.txt
 */
int main(int argc, char **argv)
{
//Check arguments : IN/OUT;
    if (argc != 13)
    {
        std::cerr << " ****************************************** " << std::endl
        << " ******************  Flutter-Shutter  *************** " << std::endl
        << " **************************************************** " << std::endl
        << "Usage: " << argv[0] << "Image.png Flag_{code} SNR     " << std::endl
        << "Flag_{NUM | ANALOG} Velocity Deltat Imacquired.png    " << std::endl
        << "ImReconstructed.png Groundtruth.png Imdiff.png        " << std::endl
        << "Code_file_name.txt  FT_code_file_name.txt             " << std::endl
        << "Input : Image.png                                     " << std::endl
        << "Image to use for the simulation                       " << std::endl
        << "(grayscale or color) PNG.                             " << std::endl
        << "  Codes for the Flutter-Shutter Simulator             " << std::endl
        << "Code : can be                                         " << std::endl
        << "0 for Snapshot (classic strategy)                     " << std::endl
        << "or 1 for the Raskar-code                              " << std::endl
        << "or 2 for rand-code, a code whose coefficients come    " << std::endl
        << "from a uniform R.V on [-1,1]                          " << std::endl
        << "3 for constant 1 code (accumulation,                  " << std::endl
        << "4 for Sinc-code (sinc code is the best code given v)  " << std::endl
        << "5 for motion-invariant photography code               " << std::endl
        << "SNR for the level 100 :100 by default                 " << std::endl
        << "Flag_{NUM | ANALOG} 0 for Numerical Flutter           " << std::endl
        << "or 1 for Analog Flutter                               " << std::endl
        << "Velocity v                                            " << std::endl
        << " Deltat : temporal sampling (float default 1)         " << std::endl
        << "Outputs :                                             " << std::endl
        << "Imacquired.png: acquired image  in PNG.               " << std::endl
        << "Imdiff: difference image (groundtruth and             " << std::endl
        << "actual reconstruction) image  in PNG.                 " << std::endl
        << "ImReconstructed: reconstructed  in PNG.               " << std::endl
        << "code_file_name: txt files that contains the           " << std::endl
        << "code coefficients  $a_k$                              " << std::endl
        << "FT_code_file :txt file that contains the Fourier      " << std::endl
        << "transform of the flutter Shutter function             " << std::endl
        << "****************************************************  " << std::endl
        << "******  Yohann Tendero, 2012  **********************  " << std::endl
        << "**************************************************** " << std::endl;
        return 1;
    }








    /** @brief The steps are the following
    * 1) Read image
    * 2) Simulate the observed
    * 3) Crop (realistic acquisition is non periodic) : of the observed
    * 4) Deconvolve (on the symetrized observed) then crop again
    * 5) Steps5- Compute RMSE, difference (residual noise), write output,etc.
    * */





    const int code_length=52; //Length of the code
    const int num_plot_points=100; //number of points used for plots
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Step1. READ IMAGE ///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    float * Image_color;  ///image
    size_t width, height, channel_number; /// width an height of the image
    if (NULL == (Image_color = read_png_f32(argv[1], &width, &height,
                                            &channel_number)))
    {
        std::cerr << "Unable to load  file " << argv[1] << std::endl;
        return 1;
    }


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //////////////////DECOMPOSE IMAGE TO DEAL WITH COLOR////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    ///Allocating & initializing
    float** restored_out = new float*[channel_number];
    float** acquired_out = new float*[channel_number];
    const int acquired_width=width-code_length;
    const int border_size=acquired_width;
    //Classic mirror symmetry adding 2*code_length in the direction of
    //the blur on each sides.
    const int restored_width=acquired_width+border_size;
    const int out_width=acquired_width-2*code_length;// avoid border effects.
    const int crop_translation=round(code_length/2);
    double** restored = new double*[channel_number];
    double* acquired = new double[width*height];
    double* acquired_crop= new double[acquired_width*height];

    // Read parameter
    int code_number=atoi(argv[2]);
    float snr=atof(argv[3]);
    int flag_num_ana_flutter=atoi(argv[4]);
    float velocity=atof(argv[5]);
    float deltat=atof(argv[6]);
    //Init of random generator here so
    //it will produce independant RV
    // for all k in the following loop
    mt_init_genrand((unsigned long int) time(NULL));
    unsigned k;



    for ( k=0; k<channel_number; k++) // LOOP OVER CHANNELS
    {



        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        /////////////// Step2. SIMULATION OF THE ACQUIRED IMAGE (NOISY) /////
        ///////////////////FOR NUMERICAL OR ANALOG FLUTTER//////////////////////
        ////////////////////////////////////////////////////////////////////////




        ///Allocating, init


        for (unsigned i=0; i<width*height; i++)
            acquired[i]=Image_color[k*width*height+i];
        ///Step2 Simulate observed
        if (flag_num_ana_flutter==0)//analog case
        {
            flutter_analog(acquired, code_number, width,
                           height, velocity, deltat, snr);
        }
        else
        {

            flutter_numerical(acquired, code_number, width, height,
                              velocity, deltat, snr);
        }

        ///Step3 Crop to avoid periodization (Step3)
        for (int column=0; column<acquired_width; column++)
        {
            for (unsigned line=0; line<height; line++)
            {
                acquired_crop[line*acquired_width+column]=
                    acquired[line*(width)+column+crop_translation];
            }
        }



        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ///////////////// Step4. SIMULATION OF THE RECONSTRUCTED IMAGE /////////
        ////////////////////////// (DECONVOLUTION) /////////////////////////////
        ////////////////////////////////////////////////////////////////////////






        restored[k] = new double[restored_width*height];

        flutter_restore(acquired_crop, code_number, acquired_width, height,
                        velocity, restored[k], border_size, deltat);



        ///Cropping for output images
        restored_out[k] = new float[out_width*height];
        acquired_out[k] = new float[out_width*height];

        for (int column=0; column<out_width; column++)
        {
            for (unsigned line=0; line<height; line++)
            {
                restored_out[k][line*out_width+column]=
                    restored[k][line*(acquired_width+border_size)+column+
                                2*crop_translation];
                acquired_out[k][line*out_width+column]=
                    acquired[line*(acquired_width+code_length)+column+
                             3*crop_translation];
            }
        }





    } //END OF THE LOO FOR RGB
    for (unsigned k=0; k<channel_number; k++)
        delete [] restored[k];
    delete [] acquired;
    delete [] acquired_crop;
    delete [] restored;

    float normalization_factor=integral_code(code_number, deltat);
    /// Transfer to standard color format and normalize the acquired
    /// for display (by $\int \alpha(t) dt$)
    float* acquired_color = new float[channel_number*out_width*height];
    float* restored_color = new float[channel_number*out_width*height];
    for (unsigned k=0; k<channel_number; k++)
        for (unsigned i=0; i<out_width*height; i++)
        {
            acquired_color[k*out_width*height+i]=acquired_out[k][i]
                                                 /normalization_factor;
            restored_color[k*out_width*height+i]=restored_out[k][i];
        }

    for (unsigned k=0; k<channel_number; k++)
    {
        delete [] acquired_out[k];
        delete [] restored_out[k];
    }
    delete [] acquired_out;
    delete [] restored_out;

// the outputs with colors are acquired_color and acquired_color



///CROP IMAGE FOR DIFFERENCE.
    float* Image_crop = new float[channel_number*out_width*height];
    /// The following realize a CROP
    for (unsigned k=0; k<channel_number; k++)
    {
        for (int column=0; column<out_width; column++)
        {
            for (unsigned line=0; line<height; line++)
            {
                Image_crop[line*out_width+column+k*height*out_width]=
                    Image_color[line*(width)+column+
                                3*crop_translation+k*height*width];
            }
        }
    }


    free(Image_color);
    ///Acquired
    write_png_f32(argv[7],acquired_color,out_width,height,channel_number);
    delete [] acquired_color;

    ///Reconstruced
    write_png_f32(argv[8],restored_color,out_width,height,channel_number);

    ///Groundtruth
    write_png_f32(argv[9],Image_crop,out_width, height, channel_number);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////////Step6. RMSE COMPUTATION //////////////////////
    ///////////////////////////////IMPOSING DYNAMIC[0,255]//////////////////////
    ////////////////////////////////////////////////////////////////////////////







    float* difference = new float[channel_number*out_width*height];
    image_difference(Image_crop,restored_color, out_width, height,
                     channel_number,difference);

    ///RMSE computation
    RMSE(difference, out_width, height,channel_number,3*code_length,0);

    ///IMPOSING DYNAMIC[0,255]
    dynamic_renormalization(difference, out_width, height,channel_number);




    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////////Step7. WRITING THE OUTPUT ////////////////////
    ////////////////////////////////////////////////////////////////////////////





    ///Difference
    write_png_f32(argv[10],difference,out_width, height, channel_number);
    delete [] difference;


    ///CSV-like files for kernel (flutter shutter fonction
    /// and its Fourier transform) to use with gnuplot or similar.

    /// 1) WRITE the code into file : code_file_name
    ofstream file_code(argv[11], ios::out | ios::trunc);
    // opening the file (erase if previous exists)
    if (file_code)
    {
        for (int k=0; k<code_length; k++)
        {
            file_code << code[atoi(argv[2])][k] <<";" << endl;
        }

        file_code.close();
    }
    else
        cerr << "Unable to open file !" << endl;

    /// 2) WRITE the Fourier transform code into file : code_file_name
    ofstream file_TFcode(argv[12], ios::out | ios::trunc);
    // opening the file (erase if previous exists)
    if (file_TFcode)
    {
        for (float xi=-M_PI;
                   xi<=M_PI;
                xi=xi+2*M_PI/num_plot_points)
        {
            file_TFcode << xi << " " <<
            abs_hat_alpha(code[atoi(argv[2])],  code_length, xi*velocity,
                          deltat) << endl;
        }

        file_TFcode.close();
    }
    else
        cerr << "Unable to open file !" << endl;



    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////////Step8. RMSE constrast invariant //////////////
    ////////////////////////////////////////////////////////////////////////////

    ///crop to speed up a little
    int difference_crop2_width=out_width-2*code_length;
    float** Image_crop2 = new float*[channel_number];
    float** restored_crop2 = new float*[channel_number];


    for (unsigned k=0; k<channel_number; k++)
    {
        Image_crop2[k] =  new float[difference_crop2_width*height];
        restored_crop2[k] = new float[difference_crop2_width*height];


        for (unsigned j=0; j<height; j++)
        {
            for (int i=0; i<difference_crop2_width; i++)
                ///AVOIDING 3 FIRSTs AND 3 LASTs COLUMNS TO AVOID
                ///BORDERS EFFECTS.
            {
                restored_crop2[k][i+(difference_crop2_width)*j]=
                    restored_color[i+out_width*j+code_length
                                   +k*height*out_width];
                Image_crop2[k][i+(difference_crop2_width)*j]=
                    Image_crop[i+out_width*j+code_length+k*height*out_width];

            }
        }

        midway(Image_crop2[k],restored_crop2[k],difference_crop2_width,height);
    }





    delete [] restored_color;
    delete [] Image_crop;


    float* difference_crop2=
        new float[channel_number*difference_crop2_width*height];
    float* Image_crop2_final =
        new float[channel_number*difference_crop2_width*height];
    float* restored_crop2_final =
        new float[channel_number*difference_crop2_width*height];



    for ( k=0; k<channel_number; k++)
        for (unsigned i=0; i<difference_crop2_width*height; i++)
        {
            Image_crop2_final[k*difference_crop2_width*height+i]=
                Image_crop2[k][i];
            restored_crop2_final[k*difference_crop2_width*height+i]=
                restored_crop2[k][i];
        }


    for (unsigned k=0; k<channel_number; k++)
    {
        delete [] Image_crop2[k];
        delete [] restored_crop2[k];
    }
    delete [] Image_crop2;
    delete [] restored_crop2;



    image_difference(Image_crop2_final,restored_crop2_final,
                     difference_crop2_width, height,channel_number,
                     difference_crop2);


    delete [] Image_crop2_final;
    delete [] restored_crop2_final;



    ///RMSE CI is the RMSE AFTER NORMALIZATION OF CONTRAST.
    ///For color images the constrast is normalized independantly
    RMSE(difference_crop2, difference_crop2_width, height, channel_number,0,1);

    delete [] difference_crop2;




    return 0;
}
