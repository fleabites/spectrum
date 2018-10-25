#!/bin/sh
make ad9361-iiostream-spectrum
./ad9361-iiostream-spectrum
./show-fft-plot.gp
