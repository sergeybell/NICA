set term postscript
set output "Phase Space My NICA PTC ypy.ps"

set pointsize 0.4
set title "Tracking Original"
set xlabel 'y'
set ylabel 'py'
set xrange [-0.02:0.02]
set yrange [-0.02:0.02]
set grid
plot 'trackone' using 5:($1==1 ? $6 : NaN) notitle with points pointtype 1, \
'trackone' using 5:($1==2 ? $6 : NaN) notitle with points pointtype 2, \
'trackone' using 5:($1==3 ? $6 : NaN) notitle with points pointtype 3, \
'trackone' using 5:($1==4 ? $6 : NaN) notitle with points pointtype 4, \
'trackone' using 5:($1==5 ? $6 : NaN) notitle with points pointtype 5, \
'trackone' using 5:($1==6 ? $6 : NaN) notitle with points pointtype 6, \
'trackone' using 5:($1==7 ? $6 : NaN) notitle with points pointtype 7, \
'trackone' using 5:($1==8 ? $6 : NaN) notitle with points pointtype 8, \
'trackone' using 5:($1==9 ? $6 : NaN) notitle with points pointtype 9, \
'trackone' using 5:($1==10 ? $6 : NaN) notitle with points pointtype 10, \
'trackone' using 5:($1==11 ? $6 : NaN) notitle with points pointtype 11, \
'trackone' using 5:($1==12 ? $6 : NaN) notitle with points pointtype 12, \
'trackone' using 5:($1==13 ? $6 : NaN) notitle with points pointtype 13, \
'trackone' using 5:($1==14 ? $6 : NaN) notitle with points pointtype 14, \
'trackone' using 5:($1==15 ? $6 : NaN) notitle with points pointtype 15, \
'trackone' using 5:($1==16 ? $6 : NaN) notitle with points pointtype 16, \
'trackone' using 5:($1==17 ? $6 : NaN) notitle with points pointtype 17, \
'trackone' using 5:($1==18 ? $6 : NaN) notitle with points pointtype 18, \
'trackone' using 5:($1==19 ? $6 : NaN) notitle with points pointtype 19, \
'trackone' using 5:($1==20 ? $6 : NaN) notitle with points pointtype 20, \
'trackone' using 5:($1==21 ? $6 : NaN) notitle with points pointtype 21, \
'trackone' using 5:($1==22 ? $6 : NaN) notitle with points pointtype 22, \
'trackone' using 5:($1==23 ? $6 : NaN) notitle with points pointtype 23, \
'trackone' using 5:($1==24 ? $6 : NaN) notitle with points pointtype 24, \
'trackone' using 5:($1==25 ? $6 : NaN) notitle with points pointtype 25, \
'trackone' using 5:($1==26 ? $6 : NaN) notitle with points pointtype 26, \
'trackone' using 5:($1==27 ? $6 : NaN) notitle with points pointtype 27, \
'trackone' using 5:($1==28 ? $6 : NaN) notitle with points pointtype 28, \
'trackone' using 5:($1==29 ? $6 : NaN) notitle with points pointtype 29, \
'trackone' using 5:($1==30 ? $6 : NaN) notitle with points pointtype 30, \
'trackone' using 5:($1==31 ? $6 : NaN) notitle with points pointtype 31, \
'trackone' using 5:($1==32 ? $6 : NaN) notitle with points pointtype 32, \
'trackone' using 5:($1==33 ? $6 : NaN) notitle with points pointtype 33, \
'trackone' using 5:($1==34 ? $6 : NaN) notitle with points pointtype 34, \
'trackone' using 5:($1==35 ? $6 : NaN) notitle with points pointtype 35, \
'trackone' using 5:($1==36 ? $6 : NaN) notitle with points pointtype 36, \
'trackone' using 5:($1==37 ? $6 : NaN) notitle with points pointtype 37, \
'trackone' using 5:($1==38 ? $6 : NaN) notitle with points pointtype 38, \
'trackone' using 5:($1==39 ? $6 : NaN) notitle with points pointtype 39, \
'trackone' using 5:($1==40 ? $6 : NaN) notitle with points pointtype 40, \
'trackone' using 5:($1==41 ? $6 : NaN) notitle with points pointtype 41, \
'trackone' using 5:($1==42 ? $6 : NaN) notitle with points pointtype 42, \
'trackone' using 5:($1==43 ? $6 : NaN) notitle with points pointtype 43, \
'trackone' using 5:($1==44 ? $6 : NaN) notitle with points pointtype 44, \
'trackone' using 5:($1==45 ? $6 : NaN) notitle with points pointtype 45, \
'trackone' using 5:($1==46 ? $6 : NaN) notitle with points pointtype 46
