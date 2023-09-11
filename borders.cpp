/*borders.cpp*/
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

/** @file borders.cpp
* @brief Mirror symmetry among columns routines
* @author Yohann Tendero <yohann.tendero@cmla.ens-cachan.fr>
*/



#include "borders.h"




/**
* @fn borders(float *Image,float *modified, int w1, int h1)
* @brief classic mirror symmetry among columns only, symmetry axis width+0.5.
* The function doubles the width.
* @param *Image  : image
* @param *modified image
* @param w1 : width of image
* @param h1 : height of image
*
*/
void borders(double *Image,double *modified, int w1, int h1)
{




////////////////////////////////////////////////////////////////////////////////
////////////////////////// Copying the image on the left////////////////////////
////////////////////////////////////////////////////////////////////////////////


    for (int colonne=0;colonne<w1;colonne++)
    {
        for (int line=0;line<h1;line++)
        {
            modified[line*(2*w1)+colonne]=Image[line*w1+colonne];

        }
    }



////////////////////////////////////////////////////////////////////////////////
////////////////////////// Right side : symmetry/////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


    for (int colonne=w1;colonne<2*w1;colonne++)
    {
        for (int line=0;line<h1;line++)
        {
            modified[line*(2*w1)+colonne]=Image[line*w1+(2*w1-colonne-1)];
        }
    }

}
