marg=0.2
#set bmargin at screen marg
#set tmargin at screen 1-marg
set lmargin at screen marg
set rmargin at screen 1-marg

plot_file = sprintf("pca/pca_%i.png", i)

stats input_file u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats input_file u 2 nooutput
miny = STATS_min
maxy = STATS_max

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set xtics out
set ytics out
set pm3d map
set pm3d interpolate 5,5
set palette cubehelix
set term pngcairo size 1200,1080 font "Arial, 28"
set grid xtics ytics nomxtics nomytics
if (maxx - minx >= 50) {
   set mxtics 5
}
if (maxy - miny >= 50) {
   set mytics 5
}

title_str = sprintf("Principal component %i", i)
set title title_str
plot_file = sprintf("pca/pca_%i.png", i)
set output plot_file
set xlabel "Residue index"
set ylabel "Residue index"
set zrange [-maxz:maxz]
set cbrange [-maxz:maxz]
# Create output 

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

splot input_file u 1:2:3 notitle

if (domains==1) {
  load "domains.gnu"
  unset xlabel
  unset ylabel
  plot_file = sprintf("pca/pca_%i.domains.png", i)
  set output plot_file
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
# Create output 
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

  splot input_file u 1:2:3 notitle
}
