set terminal epslatex color solid 8 standalone size 3,3
set output "loop.tex"
set format '$%g$'

set ytics
set format x ""
unset xlabel
unset ylabel
set xrange [0:1]
set style line 1 lt 3 lw 3 lc 'blue'
set style line 2 lt 1 lw 3 lc 'red'
set style line 3 lt 2 lw 3 lc 'green'
set xtics
unset key

set multiplot layout 1,1 scale 1.,1. offset 0,0

set style line 12 lc rgb 'black' lt 1 lw 2
set grid xtics ls 12

set xtics ('$\Gamma$' 0, 'X' 0.292893, 'M' 0.585786, '$\Gamma$' 1)
set title 'Square lattice with $\sigma_xk_x+\sigma_yk_y$ spin orbit coupling' font "Helvetica" offset 0,-0.6

plot 'd2test_loop_en.dat' u 1:4 w l ls 1, '' u 1:5 w l ls 1 , '' u 1:6 w l ls 1, '' u 1:7 w l ls 1, '' u 1:8 w l ls 1, '' u 1:9 w l ls 1, '' u 1:10 w l ls 1, '' u 1:11 w l ls 1

unset multiplot

set output
