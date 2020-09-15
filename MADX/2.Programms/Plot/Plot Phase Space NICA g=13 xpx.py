set term postscript
set output "x_px.ps"

set pointsize 0.48
set title "Tracking Original"
set xlabel 'x'
set ylabel 'px'
set xrange [-0.03:0.03]
set yrange [-0.03:0.03]
set grid
plot 'basisone' using 3:($1==1 ? $4 : NaN) notitle with points pointtype 1, \
'basisone' using 3:($1==2 ? $4 : NaN) notitle with points pointtype 2, \
'basisone' using 3:($1==3 ? $4 : NaN) notitle with points pointtype 3, \
'basisone' using 3:($1==4 ? $4 : NaN) notitle with points pointtype 4, \
'basisone' using 3:($1==5 ? $4 : NaN) notitle with points pointtype 5, \
'basisone' using 3:($1==6 ? $4 : NaN) notitle with points pointtype 6, \
'basisone' using 3:($1==7 ? $4 : NaN) notitle with points pointtype 7, \
'basisone' using 3:($1==8 ? $4 : NaN) notitle with points pointtype 8, \
'basisone' using 3:($1==9 ? $4 : NaN) notitle with points pointtype 9, \
'basisone' using 3:($1==10 ? $4 : NaN) notitle with points pointtype 10, \
'basisone' using 3:($1==11 ? $4 : NaN) notitle with points pointtype 11, \
'basisone' using 3:($1==12 ? $4 : NaN) notitle with points pointtype 12, \
'basisone' using 3:($1==13 ? $4 : NaN) notitle with points pointtype 13, \
'basisone' using 3:($1==14 ? $4 : NaN) notitle with points pointtype 14, \
'basisone' using 3:($1==15 ? $4 : NaN) notitle with points pointtype 15, \
'basisone' using 3:($1==16 ? $4 : NaN) notitle with points pointtype 16, \
'basisone' using 3:($1==17 ? $4 : NaN) notitle with points pointtype 17, \
'basisone' using 3:($1==18 ? $4 : NaN) notitle with points pointtype 18, \
'basisone' using 3:($1==19 ? $4 : NaN) notitle with points pointtype 19, \
'basisone' using 3:($1==20 ? $4 : NaN) notitle with points pointtype 20, \
'basisone' using 3:($1==21 ? $4 : NaN) notitle with points pointtype 21, \
'basisone' using 3:($1==22 ? $4 : NaN) notitle with points pointtype 22, \
'basisone' using 3:($1==23 ? $4 : NaN) notitle with points pointtype 23, \
'basisone' using 3:($1==24 ? $4 : NaN) notitle with points pointtype 24, \
'basisone' using 3:($1==25 ? $4 : NaN) notitle with points pointtype 25, \
'basisone' using 3:($1==26 ? $4 : NaN) notitle with points pointtype 26, \
'basisone' using 3:($1==27 ? $4 : NaN) notitle with points pointtype 27, \
'basisone' using 3:($1==28 ? $4 : NaN) notitle with points pointtype 28, \
'basisone' using 3:($1==29 ? $4 : NaN) notitle with points pointtype 29, \
'basisone' using 3:($1==30 ? $4 : NaN) notitle with points pointtype 30, \
'basisone' using 3:($1==31 ? $4 : NaN) notitle with points pointtype 31, \
'basisone' using 3:($1==32 ? $4 : NaN) notitle with points pointtype 32, \
'basisone' using 3:($1==33 ? $4 : NaN) notitle with points pointtype 33, \
'basisone' using 3:($1==34 ? $4 : NaN) notitle with points pointtype 34, \
'basisone' using 3:($1==35 ? $4 : NaN) notitle with points pointtype 35, \
'basisone' using 3:($1==36 ? $4 : NaN) notitle with points pointtype 36, \
'basisone' using 3:($1==37 ? $4 : NaN) notitle with points pointtype 37, \
'basisone' using 3:($1==38 ? $4 : NaN) notitle with points pointtype 38, \
'basisone' using 3:($1==39 ? $4 : NaN) notitle with points pointtype 39, \
'basisone' using 3:($1==40 ? $4 : NaN) notitle with points pointtype 40, \
'basisone' using 3:($1==41 ? $4 : NaN) notitle with points pointtype 41, \
'basisone' using 3:($1==42 ? $4 : NaN) notitle with points pointtype 42, \
'basisone' using 3:($1==43 ? $4 : NaN) notitle with points pointtype 43, \
'basisone' using 3:($1==44 ? $4 : NaN) notitle with points pointtype 44, \
'basisone' using 3:($1==45 ? $4 : NaN) notitle with points pointtype 45, \
'basisone' using 3:($1==46 ? $4 : NaN) notitle with points pointtype 46
