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

excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","x_re",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectrx.dat",HOMO,minE,maxE,dE)
excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","y_re",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectry.dat",HOMO,minE,maxE,dE)
excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","z_re",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectrz.dat",HOMO,minE,maxE,dE)



# Energy, couplings and H' in space of orbital indexes
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,"spectr1/ave_Ham_re.dat")
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,"spectr1/ave_Ham_im.dat")
excitation_spectrum.ham_map("res/0_Hprime_",tmin,tmax,"x_re",opt,scl2,"spectr1/ave_Hprime_x_re.dat")
excitation_spectrum.ham_map("res/0_Hprime_",tmin,tmax,"y_re",opt,scl2,"spectr1/ave_Hprime_y_re.dat")
excitation_spectrum.ham_map("res/0_Hprime_",tmin,tmax,"z_re",opt,scl2,"spectr1/ave_Hprime_z_re.dat")

# Same, but in space of orbital energies
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,scl1,"spectr1/1ave_Ham_re.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,scl1,"spectr1/1ave_Ham_im.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Hprime_",tmin,tmax,"x_re",opt,scl1,scl2,"spectr1/1ave_Hprime_x_re.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Hprime_",tmin,tmax,"y_re",opt,scl1,scl2,"spectr1/1ave_Hprime_y_re.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Hprime_",tmin,tmax,"z_re",opt,scl1,scl2,"spectr1/1ave_Hprime_z_re.dat")

