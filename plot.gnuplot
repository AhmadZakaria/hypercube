reset
set style fill solid
set term x11 0
plot "lengths_hist.dat" using 1:2 with boxes ti "Distribution of Eulidean Lengths"
set term x11 1
plot "dists_hist.dat" using 1:2 with boxes ti "Distribution of Distances from Space Diagonal"
set term x11 2
plot "angles_hist.dat" using 1:2 with boxes ti "Distribution of Angles (radians)"

pause -1  "Hit return to continue"
