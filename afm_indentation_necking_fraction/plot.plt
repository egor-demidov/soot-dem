#!/usr/bin/env gnuplot

set term pdf colour enhanced
set output "force-displacement.pdf"
#set yrange [-5e-10:1e-10]

# set xtics 1e-7
# set format x "%.s*10^{%T}"

set xlabel "μm"
set ylabel "nN"

plot "force-displacement-frac-92.csv" using ($1*1e6):($2*-1e9) with lines title "92% necking fraction",\
    "force-displacement-frac-95.csv" using ($1*1e6):($2*-1e9) with lines title "95% necking fraction",\
    "force-displacement-frac-98.csv" using ($1*1e6):($2*-1e9) with lines title "98% necking fraction";
