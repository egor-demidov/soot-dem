#!/usr/bin/env gnuplot

set term pdf colour enhanced
set output "experiment-and-simulation.pdf"
#set yrange [-5e-10:1e-10]

# set xtics 1e-7
# set format x "%.s*10^{%T}"

set xlabel "μm"
set ylabel "nN"

set datafile separator ","

set xrange [0:0.8]

plot "force-displacement-frac-95.csv" every 150 using ($1*1e6-0.16):($2*-1e9) with lines title "DEM",\
    "afm_data.csv" using 1:2 with lines title "AFM";