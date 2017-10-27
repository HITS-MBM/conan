set term pngcairo size 1200,1080 font "Arial, 28"
set output "pca/pca_variance.png"
set xlabel "PC #"
set ylabel "Variance (nm^2)"
set key top left
stats "pca/pca_output.txt" u 1 nooutput
minx = STATS_min
maxx = STATS_max
set xtics minx,1,maxx
plot "pca/pca_output.txt" u 1:2 w lp ps 3 pt 5 lw 2 lt 1 lc rgb "red" title "Eigenvalue",\
     "pca/pca_output.txt" u 1:3 w lp ps 3 pt 5 lw 2 lt 1 lc rgb "blue" title "Running total",\
     "pca/pca_output.txt" u 1:($3/$4) w l lw 2 lt 1 lc rgb "green" title "Total"

