set terminal postscript

set title "Glucose Detection through Noninvasive Means"

set xLabel " "
set yLabel " "

set grid

plot "Linearoutput.dat" using 2:1
