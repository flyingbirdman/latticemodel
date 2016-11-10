set terminal epslatex color solid 8 standalone size 8,3
set output "rect.tex"

set yrange [-0.5:0.5]
set xrange [-0.5:0.5]
set format '$%g$'
set xtics ('0' 0, '$\pi/2$' 0.5, '$-\pi/2$' -0.5, '$2\pi$' 2, '$3\pi/2$' 3, '$4\pi$' 4, '$5\pi' 5) offset 0,0
set ytics ('0' 0, '$\pi/2$' 0.5, '$-\pi/2$' -0.5, '$2\pi$' 2, '$3\pi/2$' 3, '$4\pi$' 4, '$5\pi' 5) offset 0,0
#set xlabel '\Large $k_0$' offset 0,0.2 font "Helvetica Bold"
#set ylabel '\Large optical lattice strength: $s$' offset -0.15,0 font "Helvetica Bold"
unset ztics
unset key

set pm3d
set pm3d map
unset surface
unset colorbox
#set view 0,0
set palette model RGB defined ( 0.0 'blue', 1.0 'white', 2.0 'red' )

set multiplot layout 1,3 scale 1.20,1.21 offset 0,0.04

set title '\Large (a) cubic' font "Helvetica" offset 0,-0.75
splot 'cubic_r_eigval.dat' u 1:2:3 w p

set title '\Large (b) hexagon' font "Helvetica" offset 0,-0.75
splot 'hexa_r_eigval.dat' u 1:2:3 w p

set title '\Large (c) kagome' font "Helvetica" offset 0,-0.75
splot 'kagome_r_eigval.dat' u 1:2:3 w p

unset multiplot
set output