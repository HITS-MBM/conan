marg=0.2
#set bmargin at screen marg
#set tmargin at screen 1-marg
set lmargin at screen marg
set rmargin at screen 1-marg

stats 'aggregate/mdmat_average_rmsf.dat' u 3 nooutput
maxr = STATS_max
stats 'aggregate/timeline.dat' u 3 nooutput
mint = STATS_min*0.001
stats 'aggregate/timeline.dat' u 4 nooutput
maxt = STATS_max*0.001

stats 'aggregate/timeline.dat' u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats 'aggregate/timeline.dat' u 2 nooutput
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

set xlabel "Residue index"
set ylabel "Residue index"

# Create output 

set title "Time of the first encounter"
set cblabel rotate by -90 "Time (ns)"
set cbrange [mint:maxt]
set output 'aggregate/timeline_first_encounter.png'
splot 'aggregate/timeline.dat' u 1:2:(0.001*$3) notitle

set palette negative
set title "Time of the last encounter"
set cbrange [mint:maxt]
set cblabel rotate by -90 "Time (ns)"
set output 'aggregate/timeline_last_encounter.png'
splot 'aggregate/timeline.dat' u 1:2:(0.001*$4) notitle

set title "Number of encounters"
set cbrange [*:*]
set cblabel rotate by -90 "#"
set output 'aggregate/num_encounter.png'
splot 'aggregate/timeline.dat' u 1:2:7 notitle

set title "Average encounter time"
set cblabel rotate by -90 "Time (ns)"
set output 'aggregate/avg_encounter.png'
splot 'aggregate/timeline.dat' u 1:2:($7==0?0:0.001*$5*$6/$7) notitle

set title "Total interaction time"
set cblabel rotate by -90 "(%)"
set output 'aggregate/interaction_lifetime.png' 
splot 'aggregate/timeline.dat' u 1:2:(100*$5) notitle

set palette positive
set title "Average contact map"
set cbrange [0:maxr]
set cblabel rotate by -90 "Distance (nm)"
set output 'aggregate/avg_mdmat.png' 
splot 'aggregate/mdmat_average_rmsf.dat' u 1:2:3 notitle

set palette negative
set title "Standard deviation in contact map"
set cbrange [*:*]
set cblabel rotate by -90 "Distance (nm)"
set output 'aggregate/stddev_mdmat.png' 
splot 'aggregate/mdmat_average_rmsf.dat' u 1:2:4 notitle


if (domains==1) {
  load "domains.gnu"
  unset xlabel
  unset ylabel
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
# Create output 
  set palette positive
  set title "Time of the first encounter"
  set cblabel rotate by -90 "Time (ns)"
  set cbrange [mint:maxt]
  set output 'aggregate/timeline_first_encounter.domains.png'
  splot 'aggregate/timeline.dat' u 1:2:(0.001*$3) notitle

  set palette negative
  set cbrange [mint:maxt]
  set title "Time of the last encounter"
  set cblabel rotate by -90 "Time (ns)"
  set output 'aggregate/timeline_last_encounter.domains.png'
  splot 'aggregate/timeline.dat' u 1:2:(0.001*$4) notitle

  set title "Number of encounters"
  set cbrange [*:*]
  set cblabel rotate by -90 "#"
  set output 'aggregate/num_encounter.domains.png'
  splot 'aggregate/timeline.dat' u 1:2:7 notitle

  set title "Average encounter time"
  set cblabel rotate by -90 "Time (ns)"
  set output 'aggregate/avg_encounter.domains.png'
  splot 'aggregate/timeline.dat' u 1:2:($7==0?0:0.001*$5*$6/$7) notitle

  set title "Total interaction time"
  set cblabel rotate by -90 "(%)"
  set output 'aggregate/interaction_lifetime.domains.png'
  splot 'aggregate/timeline.dat' u 1:2:(100*$5) notitle

  set palette positive
  set title "Average contact map"
  set cbrange [0:maxr]
  set cblabel rotate by -90 "Distance (nm)"
  set output 'aggregate/avg_mdmat.domains.png'
  splot 'aggregate/mdmat_average_rmsf.dat' u 1:2:3 notitle

  set palette negative
  set cbrange [*:*]
  set title "Standard deviation in contact map"
  set cblabel rotate by -90 "Distance (nm)"
  set output 'aggregate/stddev_mdmat.domains.png'
  splot 'aggregate/mdmat_average_rmsf.dat' u 1:2:4 notitle
}
