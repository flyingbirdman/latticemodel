#set terminal epslatex color solid 8 standalone size 4,4
#set output "./tex/rect.tex"
#set format '$%g$'

reset
set yrange [-0.5:0.5]
set xrange [-0.5:0.5]
set xtics ('0' 0, '$\pi/2$' 0.5, '$-\pi/2$' -0.5, '$2\pi$' 2, '$3\pi/2$' 3, '$4\pi$' 4, '$5\pi' 5) offset 0,0
set ytics ('0' 0, '$\pi/2$' 0.5, '$-\pi/2$' -0.5, '$2\pi$' 2, '$3\pi/2$' 3, '$4\pi$' 4, '$5\pi' 5) offset 0,0
set xlabel '\Large $k_x$' offset 0,0.2 font "Helvetica Bold"
set ylabel '\Large $k_y$' offset 1.8,0 font "Helvetica Bold"
unset ztics
unset key

set pm3d
set pm3d map
unset surface
set colorbox
set palette model RGB defined ( 0.0 'blue', 1.0 'white', 2.0 'red' )

set multiplot layout 1,1 scale 1.20,1.21 offset 0,0.04
set title 'cubic lattice with $\sigma_x k_x+\sigma_y k_y$ spin-orbit coupling' font "Helvetica" offset 0,-0.75
splot 'd2test_zone_en.dat' u 1:2:3 w p
unset multiplot

#set output