set term postscript dashed color
set output "Distortion.ps"

set pointsize 1
set xlabel 'x'
set ylabel 'MUX'
set xrange [-0.02:0.02]
set yrange [9.2:9.8]
set grid
plot 'tr-Distortion' using 2:($1==1 ? $3 : NaN) notitle with points pointtype 1, \
'tr-Distortion' using 2:($1==2 ? $3 : NaN) notitle with points pointtype 2, \
'tr-Distortion' using 2:($1==4 ? $3 : NaN) notitle with points pointtype 4
