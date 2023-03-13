set term eps enhanced color
set output "test.eps"
#set size 0.8,0.8

psize=1.

plot [][] "solution-sph.txt" u 1:2  w lp pt 6 ps psize t 'SPH', "solution-fdm.txt" every ::0::6 u 1:2 w lp pt 6 ps psize t 'FDM', "solution-fdm.txt" every ::7::13 u 1:2 w lp pt 6 ps psize t 'FDM'


set term x11
replot
