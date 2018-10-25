#!/usr/bin/gnuplot -persist

set title 'Frequency Spectrum'
set xlabel 'Frequency (MHz)'
set ylabel 'Amplitude (dB)'

#set xr [0:10]
set xtics 1

set grid

set samples 10000

plot 'fft.csv' w lines
