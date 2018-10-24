#!/usr/bin/env gnuplot

set terminal png
set output 'fft-08.png' 
set title 'Frequency Spectrum' 
set xlabel 'Frequency (Hz)' 
set ylabel "Amplitude"

#set logscale y
#set xr [0:2.5]
#set xtics (0, 0.5, 1, 1.5, 2, 2.5)
#set yr [0.0001:1]
#set ytics ('$10^{-4}$' 0.0001, '$10^{-3}$' 0.001, '$10^{-2}$' 0.01, '$10^{-1}$' 0.1, '$10^{0}$' 1)
 
#set key right top Left
#set tics nomirror
set grid

#set samples 1000

#set key spacing 1.5
unset key

#set key samplen 2 #spacing 3.5 font ",14"
#set key font ",20"
plot 'fft-8.txt' w lines

#datafile = 'standard_results_data'
#plot for [i=0:*] datafile index i using 1:2\
#with lines title columnheader(1)

