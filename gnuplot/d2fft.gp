#set terminal epslatex color solid 8 standalone size 4,4
#set output "d2lattice_fft.tex"
#set format '$%g$'
unset ztics
unset key
set size square
set xrange [-1:1]
set yrange [-1:1]
set xlabel '$\omega_x$'
set ylabel '$\omega_y$'
set pm3d
set pm3d map
set palette model RGB defined ( 0.0 'blue', 0.6 'green', 1.2 'yellow', 1.8 'red' )

set multiplot layout 1,1 scale 1.1,1.1 offset 0,0.04

splot "2dlattice_fft.dat" u 1:2:(($3>0.00001)&&($3<0.5)?$3:1/0) w p pt 5 ps 1 palette

unset multiplot
#set output

