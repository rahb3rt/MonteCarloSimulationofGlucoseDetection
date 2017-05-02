set term postscript color
set output "verdet_1.ps"
set style line  1 lt 1 lc rgb '#0c0887' # blue
set style line  2 lt 1 lc rgb '#4b03a1' # purple-blue
set style line  3 lt 1 lc rgb '#7d03a8' # purple
set style line  4 lt 1 lc rgb '#a82296' # purple
set style line  5 lt 1 lc rgb '#cb4679' # magenta
set style line  6 lt 1 lc rgb '#e56b5d' # red
set style line  7 lt 1 lc rgb '#f89441' # orange
set style line  8 lt 1 lc rgb '#fdc328' # orange
set style line  9 lt 1 lc rgb '#f0f921' # yellow
set key off
set xlabel "x-direction" 
set xlabel  offset character 3, 0, 0 font "" textcolor lt -1 norotate
set ylabel "y-direction" 
set ylabel  offset character -4, 0, 0 font "" textcolor lt -1 norotate
set zlabel "ICBS" 
set zlabel  offset character -3, 0, 0 font "" textcolor lt -1 norotate
splot 'verdet3.dat' using 9:10:3 with points palette pointsize 3 pointtype 7
set output "verdet_2.ps"
set view map
set size ratio .9

set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back
set object 1 rect fc rgb "black" fillstyle solid 1.0

splot "verdet3.dat" using 9:10:3 with points pointtype 5 pointsize 1 palette linewidth 30
