#!/usr/bin/env gnuplot

set term pdf colour enhanced
set output "gauss.pdf"
#set yrange [-5e-10:1e-10]

set xtics 1e-7

plot "force-displacement.csv" using 1:2 with lines;