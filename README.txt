# Readme file for The flutter shutter camera simulator.
Basic use :
(compile : make)

example : ./demo_flutter boat.png 1 100  0 1 1 imacquired.png imreconstructed.png
groundtruth.png imdiff.png code.txt FT_code.txt

./demo_flutter ImageIN.png Flag_{code} SNR  Flag_{NUM : 1 | ANALOG : 0} Velocity
Deltat Imacquired.png ImReconstructed.png Groundtruth.png Imdiff.png code.txt
FT_code.txt

Flag_{code} : 0 Snapshot, 1 for Agrawal, 2 for rand-code,
  3 for constant 1 code, 4 for Sinc-code (sinc code is the best code for
|v.deltat|=1), 5 for motion-invariant photography code;
Flag_{NUM : 1 | ANALOG : 0} : 1 for numerical flutter shutter, 0 for analog
flutter shutter;
SNR>0 : for level 100 (mean/std-dev) for signal dependant noise;
Velocity : of the landscape in pixel per deltat (float);
Deltat>0 : temporal sampling of the (coded) flutter shutter function;

Notice : |v.deltat|<2 for any flutter shutter code.


ImageIN.png must be PNG (color or grayscale) image (input iamge)
Output images : Imacquired.png, ImReconstructed.png, Groundtruth.png, Imdiff.png
 are  all grayscale images.
Output txt-files : code.txt FT_code.txt

example : ./demo_flutter boat.png 1 100  0 1 1 imacquired.png imrconstructed.png
groundtruth.png imdiff.png code.txt FT_code.txt

# ABOUT

* Author    : Yohann Tendero <tendero@cmla.ens-cachan.fr>
* Copyright : (C) 2012, IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides in implementation of the algorithm described in
http://www.ipol.im/pub/algo/mrt_flutter_shutter/

* the executable file is demo_flutter

This program reads and write PNG image.

- Compilation.
Automated compilation requires the make program.

- Library.
This code requires the libpng librarym and uses the io_png routines written by
Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
io_png.c is distributed under a GPL3+ or BSD license, at your
option. See the included copyright notice, conditions and disclaimer.
This code requires the Mersenne Twister pseudo-RNG code Matsumoto, Nishimura.
See the included copyright notice, conditions and disclaimer in mt19937ar.c.

for details.

# COMPILATION
1. Download the code package and extract it. Go to that directory.

2. Compile the source code (on Unix/Linux/Mac OS).
run : make

# USAGE

'demo_flutter' takes 12 parameters:
"Image.png Flag_{code} SNR Flag_{NUM | ANALOG} Velocity Deltat Imacquired.png
ImReconstructed.png Groundtruth.png Imdiff.png code.txt FT_code.txt"

(input) "Input.png" : file name of the input image (PNG file)
  Image to use for the simulation (grayscale or color) PNG.
(input) "Flag_{code}" is a number that indicate the code use.
  Codes for the Flutter-Shutter Simulator
  Code : 0 for Snapshot, 1 for Raskar-code, 2 for rand-code,
  3 for constant 1 code, 4 for Sinc-code (sinc code is the best code given v)
  5 for motion-invariant photography code.
(input) "SNR" is a positive float, (Signal-to-Noise Ratio)
  for the level 100 : 100 (default)"
(input) " Flag_{NUM | ANALOG}" is "0" or "1"
  Flutter-shutter type "0 for Analog Flutter OR 1 for Numerical Flutter"
(input) "Velocity" v is a float"
  it is the relative velocity between the camera and the landscape
  counted in pixel per time unit $\Deltat$
(input) "Deltat" (positive float) temporal sampling of the flutter shutter
function
(output) "Imacquired.png"  : file name for observed image (PNG file)
(output) "ImReconstructed.png ": file name for recovered from the observed
  image (PNG file)
(output) "Groundtruth.png" : file name for the ideal recovered image (PNG file)
(output) "Imdiff.png" : image difference between the recovered and the
  ideal recovered (PNG file)
(output) "code_file_name.txt" : file name for file that contains the code
  coefficients  $a_k$ (txt file), used for plot generation later on
(output) "FT_code_file_name.txt" : file name for file that contains the Fourier
  transform of the flutter Shutter function coming from the code, used for
(output) Standard output (printf) : RMSE and RMSE contrast invariant
  (RMSE contrast invariant is the RMSE after contrast normalization)



example : ./demo_flutter boat.png 1 100  0 1 1 imacquired.png imrconstructed.png
groundtruth.png imdiff.png code.txt FT_code.txt

# ABOUT THIS FILE

Copyright 2012 IPOL Image Processing On Line http://www.ipol.im/
Author: Yohann Tendero <tendero@cmla.ens-cachan.fr>

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.

