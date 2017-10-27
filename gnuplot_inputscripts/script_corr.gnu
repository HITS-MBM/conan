marg=0.15
#set bmargin at screen marg
#set tmargin at screen 1-marg

set xlabel "Residue index"
set ylabel "Residue index"

set grid xtics ytics nomxtics nomytics

stats 'aggregate/timeline.dat' u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats 'aggregate/timeline.dat' u 2 nooutput
miny = STATS_min
maxy = STATS_max
set grid xtics ytics nomxtics nomytics
if (maxx - minx >= 20) {
   set mxtics 5
}
if (maxy - miny >= 20) {
   set mytics 5
}

set size square
set lmargin at screen marg
set rmargin at screen 1-marg-0.05
set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax
set cbrange [-0.5:0.5]
set zrange [-0.5:0.5]
maxz=0.5
set cbtics ("< -0.5" -0.5, -0.25, 0, 0.25, "> 0.5" 0.5)
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
set pm3d interpolate 5,5
set title title_time
set cblabel rotate by -90 "Pearson correlation coefficient"
set term pngcairo size 1200,1080 font "Arial, 28"
set output outputfile
splot inputfile u 1:2:( abs($3) > maxz ? maxz*$3/abs($3) : $3) notitle


if (domains==1) {
  load "domains.gnu"
  unset xlabel
  unset ylabel
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
  set cbrange [-0.5:0.5]
  set zrange [-0.5:0.5]
  maxz=0.5
  set cbtics ("< -0.5" -0.5, -0.25, 0, 0.25, "> 0.5" 0.5)
  set xtics out
  set ytics out
  set pm3d map
  set pm3d interpolate 5,5
  set title title_time
  set cblabel rotate by -90 "Pearson correlation coefficient"
  set output outputfile[:strlen(outputfile)-3]."domains.png"
  splot inputfile u 1:2:( abs($3) > maxz ? maxz*$3/abs($3) : $3) notitle
}
