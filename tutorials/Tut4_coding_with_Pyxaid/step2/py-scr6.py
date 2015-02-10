from PYXAID import *
# Example of usage
opt = 1
scl1 = 13.60569253 # Ry to eV
scl2 = 1000.0 # some scaling to get plottable data
HOMO = 6
minE = 0.0
maxE = 10.0
dE = 0.05
tmin = 0
tmax = 50

#excitation_spectrum.
[exE, exI] = excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","x_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectrx.dat",HOMO,minE,maxE,dE)
#excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","y_re",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectry.dat",HOMO,minE,maxE,dE)
#excitation_spectrum.calculate("res/0_Ham_","_re","res/0_Hprime_","z_re",tmin,tmax,opt,scl1,scl2,"spectr1/ab_spectrz.dat",HOMO,minE,maxE,dE)


# Energy, couplings and H' in space of orbital indexes
#excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,"spectr/ave_Ham_re.dat")
#excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,"spectr/ave_Ham_im.dat")

# Same, but in space of orbital energies
#excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,scl1,"spectr/1ave_Ham_re.dat")
#excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,scl1,"spectr/1ave_Ham_im.dat")



#===================== User scripting goes here ===================

# Convert the 'naked' data in to the format for convolution
XP = []
i = 0
while i<len(exI):
    XP.append([exE[i], exI[i]])
    i = i + 1

# Now do the convolution
dx0 = dE  # original grid spacing
dx = dE  # this is new grid spacing
var = 2.0*dE  # new variance

#print XP

TM,R = pdos.convolve(XP,dx0,dx,var)

#print R
pdos.printout("spectr/ab_spectrx_conv.dat",TM,2,R)

