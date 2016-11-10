#set terminal epslatex color solid 8 standalone size 4,4
#set output "d2lattice_fft.tex"
#set format '$%g$'
unset ztics
unset key
set size square
set xrange [-1:1]
set xlabel '$\omega_x$'
set ylabel '$t$'

set multiplot layout 1,1 scale 1.1,1.1 offset 0,0.04

plot "1dlattice_fft.dat" u 1:2 w lp

unset multiplot
#set output

