set term postscript dashed color
set output "Phase Space NICA Temp Full ypy.ps"

set pointsize 0.3
set xlabel 'y, m'
set ylabel 'py, rad'
set xrange [-0.02:0.02]
set yrange [-0.02:0.02]
set grid
plot 'tracking_temp_fullone' using 5:($1==1 ? $6 : NaN) notitle with points pointtype 1, \
'tracking_temp_fullone' using 5:($1==2 ? $6 : NaN) notitle with points pointtype 2, \
'tracking_temp_fullone' using 5:($1==3 ? $6 : NaN) notitle with points pointtype 3, \
'tracking_temp_fullone' using 5:($1==4 ? $6 : NaN) notitle with points pointtype 4, \
'tracking_temp_fullone' using 5:($1==5 ? $6 : NaN) notitle with points pointtype 5, \
'tracking_temp_fullone' using 5:($1==6 ? $6 : NaN) notitle with points pointtype 6, \
'tracking_temp_fullone' using 5:($1==7 ? $6 : NaN) notitle with points pointtype 7, \
'tracking_temp_fullone' using 5:($1==8 ? $6 : NaN) notitle with points pointtype 8, \
'tracking_temp_fullone' using 5:($1==9 ? $6 : NaN) notitle with points pointtype 9, \
'tracking_temp_fullone' using 5:($1==10 ? $6 : NaN) notitle with points pointtype 10, \
'tracking_temp_fullone' using 5:($1==11 ? $6 : NaN) notitle with points pointtype 11, \
'tracking_temp_fullone' using 5:($1==12 ? $6 : NaN) notitle with points pointtype 12, \
'tracking_temp_fullone' using 5:($1==13 ? $6 : NaN) notitle with points pointtype 13, \
'tracking_temp_fullone' using 5:($1==14 ? $6 : NaN) notitle with points pointtype 14, \
'tracking_temp_fullone' using 5:($1==15 ? $6 : NaN) notitle with points pointtype 15, \
'tracking_temp_fullone' using 5:($1==16 ? $6 : NaN) notitle with points pointtype 16, \
'tracking_temp_fullone' using 5:($1==17 ? $6 : NaN) notitle with points pointtype 17, \
'tracking_temp_fullone' using 5:($1==18 ? $6 : NaN) notitle with points pointtype 18, \
'tracking_temp_fullone' using 5:($1==19 ? $6 : NaN) notitle with points pointtype 19, \
'tracking_temp_fullone' using 5:($1==20 ? $6 : NaN) notitle with points pointtype 20, \
'tracking_temp_fullone' using 5:($1==21 ? $6 : NaN) notitle with points pointtype 21, \
'tracking_temp_fullone' using 5:($1==22 ? $6 : NaN) notitle with points pointtype 22, \
'tracking_temp_fullone' using 5:($1==23 ? $6 : NaN) notitle with points pointtype 23, \
'tracking_temp_fullone' using 5:($1==24 ? $6 : NaN) notitle with points pointtype 24, \
'tracking_temp_fullone' using 5:($1==25 ? $6 : NaN) notitle with points pointtype 25, \
'tracking_temp_fullone' using 5:($1==26 ? $6 : NaN) notitle with points pointtype 26, \
'tracking_temp_fullone' using 5:($1==27 ? $6 : NaN) notitle with points pointtype 27, \
'tracking_temp_fullone' using 5:($1==28 ? $6 : NaN) notitle with points pointtype 28, \
'tracking_temp_fullone' using 5:($1==29 ? $6 : NaN) notitle with points pointtype 29, \
'tracking_temp_fullone' using 5:($1==30 ? $6 : NaN) notitle with points pointtype 30, \
'tracking_temp_fullone' using 5:($1==31 ? $6 : NaN) notitle with points pointtype 31, \
'tracking_temp_fullone' using 5:($1==32 ? $6 : NaN) notitle with points pointtype 32, \
'tracking_temp_fullone' using 5:($1==33 ? $6 : NaN) notitle with points pointtype 33, \
'tracking_temp_fullone' using 5:($1==34 ? $6 : NaN) notitle with points pointtype 34, \
'tracking_temp_fullone' using 5:($1==35 ? $6 : NaN) notitle with points pointtype 35, \
'tracking_temp_fullone' using 5:($1==36 ? $6 : NaN) notitle with points pointtype 36, \
'tracking_temp_fullone' using 5:($1==37 ? $6 : NaN) notitle with points pointtype 37, \
'tracking_temp_fullone' using 5:($1==38 ? $6 : NaN) notitle with points pointtype 38, \
'tracking_temp_fullone' using 5:($1==39 ? $6 : NaN) notitle with points pointtype 39, \
'tracking_temp_fullone' using 5:($1==40 ? $6 : NaN) notitle with points pointtype 40, \
'tracking_temp_fullone' using 5:($1==41 ? $6 : NaN) notitle with points pointtype 41, \
'tracking_temp_fullone' using 5:($1==42 ? $6 : NaN) notitle with points pointtype 42, \
'tracking_temp_fullone' using 5:($1==43 ? $6 : NaN) notitle with points pointtype 43, \
'tracking_temp_fullone' using 5:($1==44 ? $6 : NaN) notitle with points pointtype 44, \
'tracking_temp_fullone' using 5:($1==45 ? $6 : NaN) notitle with points pointtype 45, \
'tracking_temp_fullone' using 5:($1==46 ? $6 : NaN) notitle with points pointtype 46
