#set terminal epslatex color solid 8 standalone size 7,2
#set output "d2lattice_dos.tex"
#set format '$%g$'

reset
set ytics
set format x ""
unset xlabel
unset ylabel
set style line 1 lt 3 lw 3 lc 'blue'
set style line 2 lt 1 lw 3 lc 'red'
set style line 3 lt 2 lw 3 lc 'green'
set xtics
set xrange [0:8]
set yrange [-5:10]
unset key

set multiplot layout 1,1 scale 1.,1. offset 0,0

#set xtics ('$\Gamma$' 0, 'X' 0.292893, 'M' 0.585786, '$\Gamma$' 1)
#set title '\Large (a) cubic' font "Helvetica" offset 0,-0.6

plot '1dlattice_dos.dat' u 1:($3+$5+$7+$9+$11+$13+$15+$17) w l ls 1, '' u 1:($2+$4+$6+$8+$10+$12+$14+$16) w l ls 2

unset multiplot

#set output