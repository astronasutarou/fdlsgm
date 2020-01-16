#!/usr/bin/gnuplot

set xr [-1250:1250]
set yr [-1250:1250]

FILENAME='sample_xyt_output.txt'
set size ratio -1
set view equal xy
set view 0,0
set grid
splot \
  for [n=0:49] FILENAME i 0 ev ::(10*n)::(10*n+9) \
  u 1:2:3:($4-$1):($5-$2):($6-$3) w vec nohead lc (n%6+2) lw 1 not,\
  for [n=0:49] FILENAME i 0 ev ::(10*n)::(10*n+9) \
  u 4:5:6 w p pt 6 ps 0.7 lc (n%6+2) not,\
  for [n=50:99] FILENAME i 0 ev ::(10*n)::(10*n+9) \
  u 1:2:3:($4-$1):($5-$2):($6-$3) w vec nohead lc 0 lw 1 not,\
  for [n=50:99] FILENAME i 0 ev ::(10*n)::(10*n+9) \
  u 4:5:6 w p pt 4 ps 0.7 lc 0 not,\
  for [n=100:399] FILENAME i 0 ev ::(10*n)::(10*n+9) \
  u 1:2:3:($4-$1):($5-$2):($6-$3) w vec nohead lc 0 lw 0.5 not,\
  for [n=100:399] FILENAME i 0 ev ::(10*n)::(10*n+9) \
  u 4:5:6 w p pt 2 ps 0.7 lc 0 not,\
  FILENAME i 1 u 1:2:3:($4-$1):($5-$2):($6-$3) w vec lw 1 lc 7 not, \
  FILENAME i 1 u 1:2:(-1):($4-$1):($5-$2):(0) w vec lw 3 lc 2 not
