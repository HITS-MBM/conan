marg=0.15
set bmargin at screen marg
set tmargin at screen 1-marg
set lmargin at screen marg
set rmargin at screen 1-marg
set size square

stats 'aggregate/interaction_types.dat' u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats 'aggregate/interaction_types.dat' u 2 nooutput
miny = STATS_min
maxy = STATS_max

set xrange [minx-1:maxx+1]
set yrange [miny-1:maxy+1]
set cbrange [0:3]
set zrange [0:3]

if (maxx - minx >= 50) {
   set mxtics 5
}
if (maxy - miny >= 50) {
   set mytics 5
}
set xlabel "Residue index"
set ylabel "Residue index"

if (maxx - minx >= maxy - miny) {
   point_size = 100.0/(maxx-minx)
} else {
   point_size = 100.0/(maxy-miny)
}


set xtics out
set ytics out
#set pm3d map
set title "Interaction types"
set grid xtics ytics nomxtics nomytics
set term pngcairo size 1440,1080 font "Arial, 28"
unset colorbox
set key outside right reverse Left spacing 1.5
set output "aggregate/interaction_types.png"
plot "aggregate/interaction_types.dat" u 1:($3==4? $2 : 1/0) w p ps 3 pt 5 lc rgb 'orange'  title "Hydrophobic",\
     "aggregate/interaction_types.dat" u 1:($3==4? $2 : 1/0) w p ps 3 pt 5 lc rgb 'green'   title "H-bond",\
     "aggregate/interaction_types.dat" u 1:($3==4? $2 : 1/0) w p ps 3 pt 5 lc rgb 'magenta' title "Salt bridge",\
     "aggregate/interaction_types.dat" u 1:($3==1? $2 : 1/0) w p ps point_size pt 5 lc rgb 'orange'  notit,\
     "aggregate/interaction_types.dat" u 1:($3==2? $2 : 1/0) w p ps point_size pt 5 lc rgb 'green'   notit,\
     "aggregate/interaction_types.dat" u 1:($3==3? $2 : 1/0) w p ps point_size pt 5 lc rgb 'magenta' notit

if (domains==1) {
  load "domains.gnu"
  set grid noxtics noytics mxtics mytics
  set tics scale 0,1
  set xtics out
  set ytics out
  #set pm3d map
  set title "Interaction types"
  set term pngcairo size 1440,1080 font "Arial, 28"
  unset colorbox
  unset xlabel
  unset ylabel
  set key outside right reverse Left spacing 1.5
  set output "aggregate/interaction_types.domains.png"
  plot "aggregate/interaction_types.dat" u 1:($3==4? $2 : 1/0) w p ps 3 pt 5 lc rgb 'orange'  title "Hydrophobic",\
       "aggregate/interaction_types.dat" u 1:($3==4? $2 : 1/0) w p ps 3 pt 5 lc rgb 'green'   title "H-bond",\
       "aggregate/interaction_types.dat" u 1:($3==4? $2 : 1/0) w p ps 3 pt 5 lc rgb 'magenta' title "Salt bridge",\
       "aggregate/interaction_types.dat" u 1:($3==1? $2 : 1/0) w p ps point_size pt 5 lc rgb 'orange'  notit,\
       "aggregate/interaction_types.dat" u 1:($3==2? $2 : 1/0) w p ps point_size pt 5 lc rgb 'green'   notit,\
       "aggregate/interaction_types.dat" u 1:($3==3? $2 : 1/0) w p ps point_size pt 5 lc rgb 'magenta' notit
}
