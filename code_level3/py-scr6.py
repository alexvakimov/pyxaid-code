from PYXAID import *
# Example of usage
opt = 1
scl1 = 13.60569253 # Ry to eV
scl2 = 1000.0 # some scaling to get plottable data
HOMO = 102
minE = 0.0
maxE = 10.0
dE = 0.05
tmin = 0
tmax = 250

# Note a change with regard to earlier version!!! We use only imaginary part of Hprime to compute spectra. This is because Hprime 
# now contains a Hamiltonian for field-matter interaction (transition dipole moment scaled by purely imaginary constant). So "x_re" is now
# changed to "x_im", and so on
excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","x_im",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectrx.dat",HOMO,minE,maxE,dE)
excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","y_im",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectry.dat",HOMO,minE,maxE,dE)
excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","z_im",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectrz.dat",HOMO,minE,maxE,dE)


# Accordingly, we print only imaginary part of the Hprime for plotting, because the real part is zero!
# Energy, couplings and H' in space of orbital indexes
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,"spectr1/ave_Ham_re.dat")
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,"spectr1/ave_Ham_im.dat")
excitation_spectrum.ham_map("res/0_Hprime_",tmin,tmax,"x_im",opt,scl2,"spectr1/ave_Hprime_x_im.dat")
excitation_spectrum.ham_map("res/0_Hprime_",tmin,tmax,"y_im",opt,scl2,"spectr1/ave_Hprime_y_im.dat")
excitation_spectrum.ham_map("res/0_Hprime_",tmin,tmax,"z_im",opt,scl2,"spectr1/ave_Hprime_z_im.dat")

# Same, but in space of orbital energies
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,scl1,"spectr1/1ave_Ham_re.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,scl1,"spectr1/1ave_Ham_im.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Hprime_",tmin,tmax,"x_im",opt,scl1,scl2,"spectr1/1ave_Hprime_x_im.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Hprime_",tmin,tmax,"y_im",opt,scl1,scl2,"spectr1/1ave_Hprime_y_im.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Hprime_",tmin,tmax,"z_im",opt,scl1,scl2,"spectr1/1ave_Hprime_z_im.dat")

