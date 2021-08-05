set terminal pngcairo size 3500,2000 font 'Sans,8'
 
set output 'P_VE.png'


# Plot 3D de la pression
#set isosamples 50
#set hidden3d
#set cntrparam levels 10 

set dgrid3d
set pm3d; 
unset surf; 
set view map
set autoscale
set contour base; 
set cntrparam levels incremental 0,0.250,5

splot "P_GNU.txt" using 1:2:3 