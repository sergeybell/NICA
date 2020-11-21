set term postscript
set output "Distortion.ps"

set pointsize 1
set title "Tracking Original"
set xlabel 'x'
set ylabel 'MUX'
set xrange [-0.02:0.02]
set yrange [9.2:10.6]
set grid
plot 'tr-Distortion' using 2:($1==1 ? $3 : NaN) notitle with points pointtype 1, \
'tr-Distortion' using 2:($1==2 ? $3 : NaN) notitle with points pointtype 2, \
'tr-Distortion' using 2:($1==3 ? $3 : NaN) notitle with points pointtype 3
