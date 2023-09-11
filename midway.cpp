/*midway.cpp*/
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
* @file midway.cpp
* @brief Midway contrast equalization routine
* (used here for constrast invariant RMSE)
* @author Yohann Tendero <yohann.tendero@cmla.ens-cachan.fr>
*/

#include <vector>
#include <algorithm>

typedef std::pair<float,int> midway_pair;

bool pairs_order( midway_pair P1, midway_pair P2)
{
    return P1.first< P2.first;
}


/**
* @fn
* @brief Perform the midway contrast equalization between two images.
* @param I1 : image 1
* @param I2 : image 2
* @param width : width of image1 & 2
* @param height : height of image1 & 2
*/

void midway(float* I1,float* I2,int width,int height)
{


///Step1. Sorts using midway_pairs P1 & P2 (P1.first contains I1 pixel values,
/// and P1.second contains  pixel indexes (usefull after the sorts)
///(same for P2 with I2))
    std::vector<midway_pair> P1;
    std::vector<midway_pair> P2;
    P1.resize(width*height);
    P2.resize(width*height);

    for (int i=0; i<width*height; i++)
    {
        P1[i]=(std::make_pair(I1[i],i));
        P2[i]=(std::make_pair(I2[i],i));
    }

    std::sort(P1.begin(), P1.end(),pairs_order);
    std::sort(P2.begin(), P2.end(),pairs_order);


///Step2. Specification
    float temp;
    for (int i=0; i<width*height; i++)
    {
        temp=(P1[i].first+P2[i].first)/2;
        I1[P1[i].second]=temp;
        I2[P2[i].second]=temp;
    }

}

