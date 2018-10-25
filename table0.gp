_rx!/usr/bin/env gnuplot

set terminal epslatex "enhanced" color size 15.9cm,10cm
set output "standard_results.tex"
set title 'Bit Error Rate Vs Signal to Noise Ratio' 
set xlabel 'Signal to Noise Ratio $E_b / N_0$ dB' 
#set ylabel offset 2,0
set ylabel "Bit Error Rate"

set logscale y
set xr [0:2.5]
#set xtics (0, 0.5, 1, 1.5, 2, 2.5)
set yr [0.0001:1]
set ytics ('$10^{-4}$' 0.0001, '$10^{-3}$' 0.001, '$10^{-2}$' 0.01, '$10^{-1}$' 0.1, '$10^{0}$' 1)
 
set key right top Left
set tics nomirror
set grid

set samples 1000

#set key spacing 1.5

set key samplen 2 #spacing 3.5 font ",14"
set key font ",20"
plot for [IDX=0:8] 'standard_results_data' i IDX u 1:2 w lines lw 5 title columnheader(1)

#datafile = 'standard_results_data'
#plot for [i=0:*] datafile index i using 1:2\
#with lines title columnheader(1)

