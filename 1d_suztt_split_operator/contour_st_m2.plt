# Gnuplot script to plot contour of f(x, x') at fixed beta

set terminal pngcairo size 1000,800 enhanced font 'Arial,12'
set output 'contour_st_m2_w0.png'

#set title "Contour plot of f(x, x') for \\beta = 2"
set xlabel "x"
set ylabel "x'"
set zlabel "f(x,x')"
set view map
set contour base
set cntrparam levels auto 20
unset surface
set pm3d at b
set palette rgb 21,22,23

# Force gridding from scattered data
set dgrid3d 101,101,2  # (nx, ny, smoothing)
splot 'test_st_m2_w0.txt' using 1:2:3 with pm3d notitle
