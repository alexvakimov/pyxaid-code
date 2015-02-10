#
set terminal png
set output "spectr.png"

plot "ab_spectrx.dat"  using 1:2 w l  lt 1 ,\
     "ab_spectrx_conv.dat" using 1:($2/5.0) w l lt 2 
