set term postscript dashed color
set output "Phase Space X Xp.ps"

set pointsize 0.7
set xlabel 'x, m'
set ylabel 'px, rad'
set xrange [-0.02:0.02]
set yrange [-0.02:0.02]
set grid
plot 'tracking_my_nicaone' using 3:($1==1 ? $4 : NaN) notitle with points pointtype 1, \
'tracking_my_nicaone' using 3:($1==2 ? $4 : NaN) notitle with points pointtype 2, \
'tracking_my_nicaone' using 3:($1==3 ? $4 : NaN) notitle with points pointtype 3, \
'tracking_my_nicaone' using 3:($1==4 ? $4 : NaN) notitle with points pointtype 4, \
'tracking_my_nicaone' using 3:($1==5 ? $4 : NaN) notitle with points pointtype 5, \
'tracking_my_nicaone' using 3:($1==6 ? $4 : NaN) notitle with points pointtype 6, \
'tracking_my_nicaone' using 3:($1==7 ? $4 : NaN) notitle with points pointtype 7, \
'tracking_my_nicaone' using 3:($1==8 ? $4 : NaN) notitle with points pointtype 8, \
'tracking_my_nicaone' using 3:($1==9 ? $4 : NaN) notitle with points pointtype 9, \
'tracking_my_nicaone' using 3:($1==10 ? $4 : NaN) notitle with points pointtype 10, \
'tracking_my_nicaone' using 3:($1==11 ? $4 : NaN) notitle with points pointtype 11, \
'tracking_my_nicaone' using 3:($1==12 ? $4 : NaN) notitle with points pointtype 12, \
'tracking_my_nicaone' using 3:($1==13 ? $4 : NaN) notitle with points pointtype 13, \
'tracking_my_nicaone' using 3:($1==14 ? $4 : NaN) notitle with points pointtype 14, \
'tracking_my_nicaone' using 3:($1==15 ? $4 : NaN) notitle with points pointtype 15, \
'tracking_my_nicaone' using 3:($1==16 ? $4 : NaN) notitle with points pointtype 16, \
'tracking_my_nicaone' using 3:($1==17 ? $4 : NaN) notitle with points pointtype 17, \
'tracking_my_nicaone' using 3:($1==18 ? $4 : NaN) notitle with points pointtype 18, \
'tracking_my_nicaone' using 3:($1==19 ? $4 : NaN) notitle with points pointtype 19, \
'tracking_my_nicaone' using 3:($1==20 ? $4 : NaN) notitle with points pointtype 20, \
'tracking_my_nicaone' using 3:($1==21 ? $4 : NaN) notitle with points pointtype 21, \
'tracking_my_nicaone' using 3:($1==22 ? $4 : NaN) notitle with points pointtype 22, \
'tracking_my_nicaone' using 3:($1==23 ? $4 : NaN) notitle with points pointtype 23, \
'tracking_my_nicaone' using 3:($1==24 ? $4 : NaN) notitle with points pointtype 24, \
'tracking_my_nicaone' using 3:($1==25 ? $4 : NaN) notitle with points pointtype 25, \
'tracking_my_nicaone' using 3:($1==26 ? $4 : NaN) notitle with points pointtype 26, \
'tracking_my_nicaone' using 3:($1==27 ? $4 : NaN) notitle with points pointtype 27, \
'tracking_my_nicaone' using 3:($1==28 ? $4 : NaN) notitle with points pointtype 28, \
'tracking_my_nicaone' using 3:($1==29 ? $4 : NaN) notitle with points pointtype 29, \
'tracking_my_nicaone' using 3:($1==30 ? $4 : NaN) notitle with points pointtype 30, \
'tracking_my_nicaone' using 3:($1==31 ? $4 : NaN) notitle with points pointtype 31, \
'tracking_my_nicaone' using 3:($1==32 ? $4 : NaN) notitle with points pointtype 32, \
'tracking_my_nicaone' using 3:($1==33 ? $4 : NaN) notitle with points pointtype 33, \
'tracking_my_nicaone' using 3:($1==34 ? $4 : NaN) notitle with points pointtype 34, \
'tracking_my_nicaone' using 3:($1==35 ? $4 : NaN) notitle with points pointtype 35, \
'tracking_my_nicaone' using 3:($1==36 ? $4 : NaN) notitle with points pointtype 36, \
'tracking_my_nicaone' using 3:($1==37 ? $4 : NaN) notitle with points pointtype 37, \
'tracking_my_nicaone' using 3:($1==38 ? $4 : NaN) notitle with points pointtype 38, \
'tracking_my_nicaone' using 3:($1==39 ? $4 : NaN) notitle with points pointtype 39, \
'tracking_my_nicaone' using 3:($1==40 ? $4 : NaN) notitle with points pointtype 40, \
'tracking_my_nicaone' using 3:($1==41 ? $4 : NaN) notitle with points pointtype 41, \
'tracking_my_nicaone' using 3:($1==42 ? $4 : NaN) notitle with points pointtype 42, \
'tracking_my_nicaone' using 3:($1==43 ? $4 : NaN) notitle with points pointtype 43, \
'tracking_my_nicaone' using 3:($1==44 ? $4 : NaN) notitle with points pointtype 44, \
'tracking_my_nicaone' using 3:($1==45 ? $4 : NaN) notitle with points pointtype 45, \
'tracking_my_nicaone' using 3:($1==46 ? $4 : NaN) notitle with points pointtype 46, \
'tracking_my_nicaone' using 3:($1==47 ? $4 : NaN) notitle with points pointtype 47, \
'tracking_my_nicaone' using 3:($1==48 ? $4 : NaN) notitle with points pointtype 48, \
'tracking_my_nicaone' using 3:($1==49 ? $4 : NaN) notitle with points pointtype 49, \
'tracking_my_nicaone' using 3:($1==50 ? $4 : NaN) notitle with points pointtype 50, \
