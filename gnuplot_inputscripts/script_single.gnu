marg=0.15
set lmargin at screen marg
set rmargin at screen 1-marg
set size square

stats inputfile u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats inputfile u 2 nooutput
miny = STATS_min
maxy = STATS_max

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

if (maxz>0) {
  set cbrange [0:maxz]
  set zrange [0:maxz]
}

set xtics out
set ytics out
set palette cubehelix
if (cblabel_str eq 'Lifetime (%)') {
    set palette negative
}
set pm3d map

if (maxx - minx <= 250 || maxy - miny <= 250) {
    set pm3d interpolate 5,5
}

set term pngcairo size 1200,1080 font "Arial, 28"
set cblabel rotate by -90 cblabel_str

set grid xtics ytics nomxtics nomytics
if (maxx - minx >= 20) {
   set mxtics 5
}
if (maxy - miny >= 20) {
   set mytics 5
}
set xlabel label_str
set ylabel label_str

set title title_str
set output outputfile

splot inputfile u 1:2:3 notitle

if (domains==1) {

  load "domains.gnu"
  unset xlabel
  unset ylabel
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
  set output outputfile[:strlen(outputfile)-3]."domains.png"
  replot

}
