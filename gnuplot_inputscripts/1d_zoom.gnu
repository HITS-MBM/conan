set xlabel "Time (ns)"
set ylabel "Distance (nm)"
set term pngcairo size 1440,1080 font "Arial, 28"
set output outputfile
set format y "%.2f"
set title title_plot
set key  outside
plot inputfile u ($1/1000):2 w l lw 1 lc rgb 'black' notitle,\
     ""        u ($1/1000):($3==1?$2:1/0) w l lw 2.5 lc rgb "red" title "interaction",\
     inter_high w l lw 1 lc rgb "blue" dt 1 title "cutoffs", inter_low w l lw 1 lc rgb "blue" dt 1 notitle
