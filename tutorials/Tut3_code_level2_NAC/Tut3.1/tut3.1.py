import os
import sys
from PYXAID.pyxaid_core import *

def run_test():

    # Parameters
    dt = 0.1 # Angstrom
    # Set #1
    minband = 1  - 1
    nocc = 3 - 1
    maxband = 10 - 1

    # Set #2
#    minband = 1  - 1
#    nocc = 3 - 1
#    maxband = 5 - 1

    # Set #3
#    minband = 3 - 1  
#    nocc = 3 - 1
#    maxband = 4 - 1
    

    EXE_CONVERT = "/software/group/oprezhdo_group/espresso-5.0.2/bin/iotk"

    # Containers for WFCs
    curr_wfc0 = wfc()
    curr_tmp0 = wfc()
    next_wfc0 = wfc()
    next_tmp0 = wfc()


    # Convert exportent wavefunction from binary to text format - using IOTK utility of QE
    os.system("%s convert x0.export/wfc.1 x0.export/wfc.1.xml" % EXE_CONVERT )
    os.system("%s convert x1.export/wfc.1 x1.export/wfc.1.xml" % EXE_CONVERT )

    # Read one wavefunction - at time t
    curr_tmp0.QE_read_acsii_index("x0.export/index.xml" )
    curr_tmp0.QE_read_acsii_wfc("x0.export/wfc.1.xml" )


    # Read the next wavefunction - at time t + dt
    next_tmp0.QE_read_acsii_index("x1.export/index.xml" )
    next_tmp0.QE_read_acsii_wfc("x1.export/wfc.1.xml" )

    # Trim each wavefunction to the active space chosen
    curr_wfc0 = wfc(curr_tmp0,minband,nocc,curr_tmp0,nocc+1,maxband)
    next_wfc0 = wfc(next_tmp0,minband,nocc,next_tmp0,nocc+1,maxband)


    #============================================
    # Try different pre-conditioning schemes:  
    do_complete = 1  
    # just choose one of the following options:
    wfc_preprocess = ""
#    wfc_preprocess = "normalize"
#    wfc_preprocess = "complete"
#    wfc_preprocess = "restore"

    if wfc_preprocess=="normalize": # just normalize
        curr_wfc0.normalize()
        next_wfc0.normalize()

    # This option is made defaut in the latest version
    elif wfc_preprocess=="complete": # complete wfc with complex conjugate part,
                                     # the result is then normalized
        curr_wfc0.complete()
        next_wfc0.complete()

    # The following option exists - but may be very slow or not work   
    # a better diagonalization procedure is needed
    # It is disabled in the new version
    elif wfc_preprocess=="restore":  # restore real wfc from the projection, the result is normailzed
                                     # optionally the wfc can be completed with complex conjugate
                                     # before restoring
        curr_wfc0.restore(0,do_complete) # first 0 - k-point, second 1 - do complete wfc
        next_wfc0.restore(0,do_complete)

   
    #============================================
    # Compute Hamiltonian
    ham(curr_wfc0,next_wfc0,0,0, 0,maxband-minband, dt,"Ham_original_"+wfc_preprocess )
#    ham(curr_wfc0,next_wfc0,0,0, 0,maxband-minband, dt,"Ham_original_"+wfc_preprocess+"_5x5" )
#    ham(curr_wfc0,next_wfc0,0,0, 0,maxband-minband, dt,"Ham_original_"+wfc_preprocess+"_2x2" )



run_test()



