max(x,y) = (x<y ? y:x)
stats 'average_rmsf.diff.dat' u 3 nooutput
maxr_meas = max (-STATS_min, STATS_max)
stats 'average_rmsf.diff.dat' u 4 nooutput
maxs_meas = max (-STATS_min, STATS_max)
marg=0.2

stats 'average_rmsf.diff.dat' u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats 'average_rmsf.diff.dat' u 2 nooutput
miny = STATS_min
maxy = STATS_max

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set lmargin at screen marg
set rmargin at screen 1-marg

set xtics out
set ytics out
set pm3d map
set pm3d interpolate 5,5
set palette cubehelix
set term pngcairo size 1200,1080 font "Arial, 28"
set grid xtics ytics nomxtics nomytics

if (maxx - minx >= 20) {
   set mxtics 5
}
if (maxy - miny >= 20) {
   set mytics 5
}
set xlabel "Residue index / ".title_A
set ylabel "Residue index / ".title_B

# Create output 

if (maxx - minx == maxy-miny) {

  set arrow from minx,miny to maxx,maxy lw 1 lt 3 lc rgb 'red' front nohead
  
  set title "Time of the first encounter"
  set cblabel rotate by -90 "Time (ns)"
  set output 'timeline_first_encounter.tri.png'
  splot 'timeline.tri.dat' u 1:2:(0.001*$3) notitle
  
  set title "Time of the last encounter"
  set cblabel rotate by -90 "Time (ns)"
  set output 'timeline_last_encounter.tri.png'
  splot 'timeline.tri.dat' u 1:2:(0.001*$4) notitle
  
  set palette negative
  set title "Total interaction time"
  set cblabel rotate by -90 "(%)"
  set output 'interaction_lifetime.tri.png' 
  splot 'timeline.tri.dat' u 1:2:(100*$5) notitle
  
  set palette positive
  set title "Average contact map"
  set cblabel rotate by -90 "Distance (nm)"
  set output 'avg_mdmat.tri.png' 
  splot 'average_rmsf.tri.dat' u 1:2:3 notitle
  
  set palette negative
  set title "Standard deviation in contact map"
  set cblabel rotate by -90 "Distance (nm)"
  set output 'stddev_mdmat.tri.png' 
  splot 'average_rmsf.tri.dat' u 1:2:4 notitle
  set palette positive
}

unset arrow

set xlabel "Residue index"
set ylabel "Residue index"

set cbrange [-maxt:maxt]
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
#set palette defined (-maxt "red", 0 "white", maxt "blue")
set title "Change in time of the first encounter"
set cblabel rotate by -90 "Time (ns)"
set output 'timeline_first_encounter.diff.png'
splot 'timeline.diff.dat' u 1:2:(0.001*$3) notitle

set title "Change in time of the last encounter"
set cblabel rotate by -90 "Time (ns)"
set output 'timeline_last_encounter.diff.png'
splot 'timeline.diff.dat' u 1:2:(0.001*$4) notitle

set cbrange [-100:100]
set title "Change in total interaction time"
set cblabel rotate by -90 "(%)"
set output 'interaction_lifetime.diff.png' 
splot 'timeline.diff.dat' u 1:2:(100*$5) notitle


set cbrange [-maxr_meas:maxr_meas]
set title "Change in average distances"
set cblabel rotate by -90 "Distance (nm)"
set output 'avg_mdmat.diff.png' 
splot 'average_rmsf.diff.dat' u 1:2:3 notitle

set cbrange [-maxs_meas:maxs_meas]
set title "Change in standard deviation in contact map"
set cblabel rotate by -90 "Distance (nm)"
set output 'stddev_mdmat.diff.png' 
splot 'average_rmsf.diff.dat' u 1:2:4 notitle                

if (domains==1) {
  load "domains.gnu"
  unset xlabel
  unset ylabel
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
# Create output 
  set cbrange [-maxt:maxt]                                       
  set title "Change in time of the first encounter"
  set cblabel rotate by -90 "Time (ns)"
  set output 'timeline_first_encounter.diff.domains.png'
  splot 'timeline.diff.dat' u 1:2:(0.001*$3) notitle

  set title "Change in time of the last encounter"
  set cblabel rotate by -90 "Time (ns)"
  set output 'timeline_last_encounter.diff.domains.png'
  splot 'timeline.diff.dat' u 1:2:(0.001*$4) notitle

  set cbrange [-100:100]
  set title "Change in total interaction time"
  set cblabel rotate by -90 "(%)"
  set output 'interaction_lifetime.diff.domains.png'
  splot 'timeline.diff.dat' u 1:2:(100*$5) notitle

  set cbrange [-maxr_meas:maxr_meas]
  set title "Change in average distances"
  set cblabel rotate by -90 "Distance (nm)"
  set output 'avg_mdmat.diff.domains.png'
  splot 'average_rmsf.diff.dat' u 1:2:3 notitle

  set cbrange [-maxs_meas:maxs_meas]
  set title "Change in standard deviation in contact map"
  set cblabel rotate by -90 "Distance (nm)"
  set output 'stddev_mdmat.diff.domains.png'
  splot 'average_rmsf.diff.dat' u 1:2:4 notitle
}
