#!/usr/bin/env python3

import subprocess
import argparse


# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("flutter_type", type=int)
ap.add_argument("code", type=int)
ap.add_argument("code1", type=int)
ap.add_argument("s3", type=int)
ap.add_argument("s4", type=float)
args = ap.parse_args()


if args.flutter_type == 0:
    p = ['demo_flutter', 'input_0.png', str(args.code), str(args.s3), '1', str(args.s4), '1', \
        'Imacquired.png', 'ImReconstructed.png', 'Groundtruth.png', 'Imdiff.png', \
            'code.txt', 'TF_code.txt']
    subprocess.run(p)
else:
    p = ['demo_flutter', 'input_0.png', str(args.code1), \
        str(args.s3), '0', str(args.s4), '1', \
        'Imacquired.png', 'ImReconstructed.png', 'Groundtruth.png', 'Imdiff.png', \
            'code.txt', 'TF_code.txt']
    subprocess.run(p)

p = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE,)
p.stdin.write(b'set lmargin 10 ; set rmargin 10\n')
p.stdin.write(b'set xlabel "k" \n')
p.stdin.write(b'set ylabel "Gain" \n')
p.stdin.write(b'set sample 1000  \n')
p.stdin.write(b'set terminal png \n')
p.stdin.write(b'set output "'+'code.png" \n')
p.stdin.write(b'set key on  below box title "Legend" \n')
p.stdin.write(b'set key left Left reverse \n ' )
p.stdin.write(b'plot "'+'code.txt" with steps title \
"Flutter shutter function" linewidth 1.5 \n' )
p.stdin.write(b'exit \n')

p = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE,)
p.stdin.write(b'set lmargin 10 ;  set rmargin 10 \n')
p.stdin.write(b'set xlabel "xi" \n')
p.stdin.write(b'set ylabel "Fourier tranform (modulus)" \n')
p.stdin.write(b'set sample 1000  \n')
p.stdin.write(b'set terminal png \n')
p.stdin.write(b'set output "' +'TF_code.png" \n')
p.stdin.write(b'set key on  below box title "Legend" \n')
p.stdin.write(b'set key left Left reverse \n ' )
p.stdin.write(b'plot "' +'TF_code.txt" using 1:2 with \
lines title "Fourier transform (modulus) of the flutter function"\
linewidth 1.5 \n' )

