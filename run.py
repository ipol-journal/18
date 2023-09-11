#!/usr/bin/env python3

import subprocess
import argparse


# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("flutter_type", type=str)
ap.add_argument("code", type=int)
ap.add_argument("code1", type=int)
ap.add_argument("s3", type=int)
ap.add_argument("s4", type=float)
args = ap.parse_args()



if args.flutter_type == 'Numerical':
    p = ['demo_flutter', 'input_0.sel.png', str(args.code), str(args.s3), '1', str(args.s4), '1', \
        'Imacquired.png', 'ImReconstructed.png', 'Groundtruth.png', 'Imdiff.png', \
            'code.txt', 'TF_code.txt']
    subprocess.run(p)
else:
    p = ['demo_flutter', 'input_0.sel.png', str(args.code1), \
        str(args.s3), '0', str(args.s4), '1', \
        'Imacquired.png', 'ImReconstructed.png', 'Groundtruth.png', 'Imdiff.png', \
            'code.txt', 'TF_code.txt']
    subprocess.run(p)

p = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE,)
p.stdin.write( ' set lmargin 10 ;  set rmargin 10 \n')
p.stdin.write( ' set xlabel "k" \n')
p.stdin.write( ' set ylabel "Gain" \n')
p.stdin.write( ' set sample 1000  \n')
p.stdin.write( ' set terminal png \n')
p.stdin.write( ' set output "'+'code.png" \n')
p.stdin.write( ' set key on  below box title "Legend" \n')
p.stdin.write( ' set key left Left reverse \n ' )
p.stdin.write( ' plot "'+'code.txt" with steps title \
"Flutter shutter function" linewidth 1.5 \n' )
p.stdin.write(' exit \n')

p = subprocess.Popen(['gnuplot','-p'],
                shell=True,
                stdin=subprocess.PIPE,
                )
p.stdin.write( ' set lmargin 10 ;  set rmargin 10 \n')
p.stdin.write( ' set xlabel "xi" \n')
p.stdin.write( ' set ylabel "Fourier tranform (modulus)" \n')
p.stdin.write( ' set sample 1000  \n')
p.stdin.write( ' set terminal png \n')
p.stdin.write( ' set output "' +'TF_code.png" \n')
p.stdin.write( ' set key on  below box title "Legend" \n')
p.stdin.write( ' set key left Left reverse \n ' )
p.stdin.write( ' plot "' +'TF_code.txt" using 1:2 with \
lines title "Fourier transform (modulus) of the flutter function"\
linewidth 1.5 \n' )

