marg=0.15
#set bmargin at screen marg
#set tmargin at screen 1-marg
set lmargin at screen marg
set rmargin at screen 1-marg
set size square
set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

stats 'aggregate/timeline.dat' u 1 nooutput
minx = STATS_min
maxx = STATS_max
stats 'aggregate/timeline.dat' u 2 nooutput
miny = STATS_min
maxy = STATS_max



set cbrange [0:maxz]
set zrange [0:maxz]

set xtics out
set ytics out
set palette cubehelix
set pm3d map
set term pngcairo size 1200,1080 font "Arial, 28"
set cblabel rotate by -90 "Distance (nm)"

set grid xtics ytics nomxtics nomytics
if (maxx - minx >= 20) {
   set mxtics 5
}
if (maxy - miny >= 20) {
   set mytics 5
}

if (maxx - minx <= 250 || maxy - miny <= 250) {
    set pm3d interpolate 5,5
}

set xlabel "Residue index"
set ylabel "Residue index"

time = begin*0.001
delta= 0.001*dt

do for [i=0:nframes-1] { 
  title_time = sprintf ("t = %7.3f (ns)",time)
  set title title_time
  outputfile = sprintf('frames/%05d.png',i)
  set output outputfile
  inputfile = sprintf('matrices/%05d.dat',i)
  splot inputfile u 1:2:3 notitle
  time = time + delta
}

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
set cblabel rotate by -90 "Distance difference (nm)"
set zrange [-maxdr:maxdr]
set cbrange [-maxdr:maxdr]

time = begin*0.001
delta= 0.001*dt

if (dr_mode%2 == 1) {
  do for [i=1:nframes-1] {
    time = time + delta
    title_time = sprintf ("Change up to t = %7.3f (ns)",time)
    set title title_time
    outputfile = sprintf('frames/%05d_dr.init.png',i)
    set output outputfile
    inputfile = sprintf('matrices/%05d_dr.init.dat',i)
    splot inputfile u 1:2:( abs($3) > maxdr ? maxdr*$3/abs($3) : $3) notitle
  }
}

time = begin*0.001
delta= 0.001*dt

if (dr_mode > 1) {
  do for [i=1:nframes-1] {
    time = time + delta
    title_time = sprintf ("Change at t = %7.3f (ns)",time)
    set title title_time
    outputfile = sprintf('frames/%05d_dr.prev.png',i)
    set output outputfile
    inputfile = sprintf('matrices/%05d_dr.prev.dat',i)
    splot inputfile u 1:2:( abs($3) > maxdr ? maxdr*$3/abs($3) : $3) notitle
  }
}


if (domains==1) {
  set autoscale xfixmin
  set autoscale xfixmax
  set autoscale yfixmin
  set autoscale yfixmax
  unset xlabel
  unset ylabel
  set cbrange [0:maxz]
  set zrange [0:maxz]

  set xtics out
  set ytics out
  set palette cubehelix

  load "domains.gnu"
  set xtics out
  set ytics out
  set palette cubehelix
  set pm3d map
  set term pngcairo size 1200,1080 font "Arial, 28"
  set cblabel rotate by -90 "Distance (nm)"

  set grid noxtics noytics mxtics mytics
  set tics scale 0,1

  time = begin*0.001
  delta= 0.001*dt

  do for [i=0:nframes-1] {
    title_time = sprintf ("t = %7.3f (ns)",time)
    set title title_time
    outputfile = sprintf('frames/%05d.domains.png',i)
    set output outputfile
    inputfile = sprintf('matrices/%05d.dat',i)
    splot inputfile u 1:2:3 notitle
    time = time + delta
  }

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
  set cblabel rotate by -90 "Distance difference (nm)"
  set zrange [-maxdr:maxdr]
  set cbrange [-maxdr:maxdr]

  time = begin*0.001
  delta= 0.001*dt

  if (dr_mode%2 == 1) {
    do for [i=1:nframes-1] {
      time = time + delta
      title_time = sprintf ("Change up to t = %7.3f (ns)",time)
      set title title_time
      outputfile = sprintf('frames/%05d_dr.init.domains.png',i)
      set output outputfile
      inputfile = sprintf('matrices/%05d_dr.init.dat',i)
      splot inputfile u 1:2:( abs($3) > maxdr ? maxdr*$3/abs($3) : $3) notitle
    }
  }

  if (dr_mode > 1) {
    do for [i=1:nframes-1] {
      time = time + delta
      title_time = sprintf ("Change at t = %7.3f (ns)",time)
      set title title_time
      outputfile = sprintf('frames/%05d_dr.prev.domains.png',i)
      set output outputfile
      inputfile = sprintf('matrices/%05d_dr.prev.dat',i)
      splot inputfile u 1:2:( abs($3) > maxdr ? maxdr*$3/abs($3) : $3) notitle
    }
  }

}
