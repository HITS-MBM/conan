marg=0.15
#set bmargin at screen marg
#set tmargin at screen 1-marg

stats inputfile u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats inputfile u 2 nooutput
miny = STATS_min
maxy = STATS_max

if (domains==1) {
  load "domains.gnu"
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
} else {
  set grid xtics ytics nomxtics nomytics
  if (maxx - minx >= 20) {
    set mxtics 5
  }
  if (maxy - miny >= 20) {
     set mytics 5
  }
  set xlabel "Residue index"
  set ylabel "Residue index"
}

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set lmargin at screen marg
set rmargin at screen 1-marg-0.05
set cbrange [-maxz:maxz]
set zrange [-maxz:maxz]
set size square
set xtics out
set ytics out
set palette defined (-5 "#7f3b08",\
                     -4 "#b35806",\
                     -3 "#e08214",\
                     -2 "#fdb863",\
                     -1 "#fee0b6",\
                      0 "#ffffff",\
                      1 "#d8daeb",\
                      2 "#b2abd2",\
                      3 "#8073ac",\
                      4 "#542788",\
                      5 "#2d004b")
set pm3d map
if (maxx - minx <= 250 || maxy - miny <= 250) {
    set pm3d interpolate 5,5
}

set title title_time
set cblabel rotate by -90 "Distance difference (nm)"
set term pngcairo size 1200,1080 font "Arial, 28"
set output outputfile
splot inputfile u 1:2:( abs($3) > maxz ? maxz*$3/abs($3) : $3) notitle
