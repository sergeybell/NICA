set term postscript dashed color
set output "Dyn.Apperture Qy.ps"

set pointsize 0.5
set xlabel 'Qy'
set ylabel 'Dyn. Aper'
set xrange [9.0:9.9]
set yrange [0:120]

set grid
plot 'Dyn.Apperture Qy' using 1:2 notitle with line, \
'Dyn.Apperture Qy' using 1:2 notitle with points pointtype 7
