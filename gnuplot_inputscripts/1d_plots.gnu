set xlabel "Time (ns)"
set ylabel "RMSD in contact map wrt previous frame (nm)"
set term pngcairo size 1440,1080 font "Arial, 28"
set output "aggregate/rmsd_contact_map_previous.png"
set format y "%.3f"

plot "aggregate/time_mdmat_rmsd.dat" u 1:2 w l lw 2.5 lc rgb 'red' notitle

set ylabel "RMSD in contact map wrt initial frame (nm)"
set output "aggregate/rmsd_contact_map_initial.png"
set format y "%.3f"
plot "aggregate/time_mdmat_rmsd.dat" u 1:3 w l lw 2.5 lc rgb 'red' notitle


asymm=system(" [ ! -e aggregate/total_interactions.y.dat ]; echo $?")

if (asymm) {
   print "Plotting interactions from an asymmetric run!"

   reset
   set term pngcairo size 1440,1080 font "Arial, 28"

   stats "aggregate/total_interactions.x.dat" u 1:2 nooutput
   minx = STATS_min_x
   maxx = STATS_max_x

   set xlabel "Residue index"
   set ylabel "Average total interactions formed"
   set xrange [minx-1:maxx+1]
   set output "aggregate/total_interactions.x.png"
   plot "aggregate/total_interactions.x.dat" u 1:2 w lp lw 2.5 lc rgb 'red' notitle

   reset
   stats "aggregate/total_interactions.y.dat" u 1:2 nooutput
   miny = STATS_min_x
   maxy = STATS_max_x

   set xlabel "Residue index"
   set ylabel "Average total interactions formed"
   set xrange [miny-1:maxy+1]
   set output "aggregate/total_interactions.y.png"
   plot "aggregate/total_interactions.y.dat" u 1:2 w lp lw 2.5 lc rgb 'red' notitle

   if (domains==1) {
     load "domains_1d.x.gnu"
     set grid noytics nomytics noxtics mxtics
     set tics scale 0,1
     unset xlabel
     set xrange [minx-1:maxx+1]
     set output "aggregate/total_interactions.x.domains.png"
     plot "aggregate/total_interactions.x.dat" u 1:2 w lp lw 2.5 lc rgb 'red' notitle

     load "domains_1d.y.gnu"
     set xrange [miny-1:maxy+1]
     set output "aggregate/total_interactions.y.domains.png"
     plot "aggregate/total_interactions.y.dat" u 1:2 w lp lw 2.5 lc rgb 'red' notitle
}

} else {

stats "aggregate/local_interactions.dat" u 1 nooutput
minres = STATS_min
maxres = STATS_max

set xlabel "Residue index"
set ylabel "Average local interaction time (%)"
set format y "%.2f"
set xrange [minres-1:maxres+1]
set output "aggregate/local_interactions.png"
plot "aggregate/local_interactions.dat" u 1:(100*$2) w lp lw 2.5 lc rgb 'red' notitle


if (domains==1) {
  load "domains_1d.gnu"
  set grid noytics nomytics noxtics mxtics
  set tics scale 0,1
  unset xlabel
  set xrange [minres-1:maxres+1]
  set output "aggregate/local_interactions.domains.png"
  plot "aggregate/local_interactions.dat" u 1:(100*$2) w lp lw 2.5 lc rgb 'red' notitle
}

}
