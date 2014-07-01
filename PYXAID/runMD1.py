#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
from pyxaid_core import *


def get_value(params,key,default,typ):
# Function to extract parameter from the dictionary
    # Try to get value from the dictionary
    str_val = "None"
    if key in params:
        if params[key]!=None:
            str_val = params[key]

    # If nothing found - use default value
    if str_val!="None":
        pass  # All is ok
    else: 
        str_val = default
        print "Warning: Parameter with key = %s does not exist in dictionary" % key
        print "Using the default value of %s" % default

    # Convert string to desired data type
    if typ=="s":
        return str_val
    elif typ=="f":
        return float(str_val)
    elif typ=="i":
        return int(float(str_val))


def runMD(params):
#------------ Read the parameters -----------------
# Parameters meaning
# pp_type - pseudopotential type: US - ultra-soft, NC - norm-conserving, PAW - projector-augmented waves
# wd - working directory, where all output (working) files will be written
# rd - results directory, where all final results (energy, NAC, H', etc.) will be written by default it will be set to wd
# This MD uses corrected NAC method

    print "Starting runMD"

    # Now try to get parameters from the input
    BATCH_SYSTEM = get_value(params,"BATCH_SYSTEM","srun","s")  # either "srun" (for SLURM) or "mpirun" (for PBS)
    NP = get_value(params,"NP","1","i")
    EXE = get_value(params,"EXE","","s")
    EXE_EXPORT = get_value(params,"EXE_EXPORT","","s")
    EXE_CONVERT = get_value(params,"EXE_CONVERT","","s")  # this is the path to iotk executable
    start_indx = get_value(params,"start_indx","0","i")
    stop_indx = get_value(params,"stop_indx","1","i")
    dt = get_value(params,"dt","1.0","f") # time step in fs - rescale NAC if actual dt is different
    pp_type = get_value(params,"pp_type","NC","s")
    wd = get_value(params,"wd","wd","s")
    rd = get_value(params,"rd",wd,"s")
    minband = get_value(params,"minband",1,"i")
    maxband = get_value(params,"maxband",2,"i")
    nocc = get_value(params,"nocc",1,"i")
    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    prefix0 = get_value(params,"prefix0","x0.scf","s")
    prefix1 = get_value(params,"prefix1","x1.scf","s")

    wfc_preprocess = get_value(params,"wfc_preprocess","restore","s") # variants are: normalize, complete, restore
    do_complete = get_value(params,"do_complete",1,"i") # argument for option "restore"
    
    compute_Hprime = get_value(params,"compute_Hprime",0,"i") # transition dipole moments


    # Sanity/Convention check
    if(minband<=0): 
        print "Error: minband should be >0, current value of minband = ",minband
        sys.exit(0)
    if(minband>maxband):
        print "Error: minband must be smaller or equal to maxband. Current values: minband = ",minband," maxband = ",maxband
        sys.exit(0)
    if(nocc>maxband):
        print "Error: nocc must be smaller or equal to maxband. Current values: nocc = ",nocc," maxband = ",maxband
        sys.exit(0)
    if(nocc<minband):
        print "Error: nocc must be larger or equal to minband. Current values: nocc = ",nocc," minband = ",minband
        sys.exit(0)

    # Convert minband, maxband and nocc from external (QE-consistent, starting from 1) convetion
    # to the internal one (starting from 0)    
    minband = minband - 1
    maxband = maxband - 1
    nocc = nocc - 1



    # Initialize variables
    curr_index = start_indx - 1
    t = start_indx

    print "In runMD: current working directory for python: ",os.getcwd()
    print "In runMD: current working directory for sh:",os.system("echo pwd")

    os.system("mkdir %s" % wd)  # Create the working directory where all output files will be written
                                # results directory should already exist

    while t<=stop_indx:
        print "t= ", t
        print "initializing variables"
        # Initialize variables
        # Wavefunctions for the neutral system
        curr_wfc0 = wfc()
        next_wfc0 = wfc()
        # Auxiliary wavefunction for the neutral system
        curr_tmp0 = wfc()
        next_tmp0 = wfc()

        # Wavefunctions for the charged system - will be used if nac_method==1
        curr_wfc1 = wfc()
        next_wfc1 = wfc()
        # Auxiliary wavefunction for the charged system
        curr_tmp1 = wfc()
        next_tmp1 = wfc()
 

        if t==start_indx:
            print "Starting first point in this batch"
            # Run calculations
            os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix0,t,prefix0,t) )
            os.system( "%s -n %s %s < x0.exp.in > x0.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )
            # Create temporary directory
            os.system("mkdir %s/curr0" % wd )
            # Copy some results to that directory
            os.system( "mv %s.%d.out %s/curr0" % (prefix0,t, wd) )
            os.system( "mv *.wfc* %s/curr0" % wd )
            os.system( "mv x0.export %s/curr0" % wd ) # "x0" - corresponds to x0 as a prefix in input files


            if nac_method>=1:  # In addition compute the charged system
                os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix1,t,prefix1,t) )
                os.system( "%s -n %s %s < x1.exp.in > x1.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

                os.system("mkdir %s/curr1" % wd )
                os.system( "mv %s.%d.out %s/curr1" % (prefix1,t, wd) )
                os.system( "mv *.wfc* %s/curr1" % wd )
                os.system( "mv x1.export %s/curr1" % wd ) # "x1" - corresponds to x1 as a prefix in input files
        

            #sys.exit(0)

        elif t>start_indx:
            print "Continuing with other points in this batch"
            # Run calculations
            os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix0,t,prefix0,t) )
            os.system( "%s -n %s %s < x0.exp.in > x0.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )
            # Create temporary directory
            os.system("mkdir %s/next0" % wd )
            # Copy some results in that directory
            os.system( "mv %s.%d.out %s/next0" % (prefix0,t, wd) )
            os.system( "mv *.wfc* %s/next0" % wd )
            os.system( "mv x0.export %s/next0" % wd )

            if nac_method>=1:
                os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix1,t,prefix1,t) )
                os.system( "%s -n %s %s < x1.exp.in > x1.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

                os.system("mkdir %s/next1" % wd )
                os.system( "mv %s.%d.out %s/next1" % (prefix1,t, wd) )
                os.system( "mv *.wfc* %s/next1" % wd )
                os.system( "mv x1.export %s/next1" % wd )

            #sys.exit(0)  # for test purpuses

        else:
            pass


        # Now general part - from current and next wavefunctions calculate NACs:
        if curr_index>=start_indx:
            print "Generate NAC from WFCs at two adjacent points"
            # Read the N-electron wavefunction descriptions
            if nac_method>=0:
                # Convert binary files to xml - this is needed in newest version of QE
                # becuase pw_export will produce binary files in any case (bug?)
                os.system("%s convert %s/curr0/x0.export/wfc.1 %s/curr0/x0.export/wfc.1.xml" % (EXE_CONVERT,wd,wd))
                os.system("%s convert %s/next0/x0.export/wfc.1 %s/next0/x0.export/wfc.1.xml" % (EXE_CONVERT,wd,wd))

                curr_tmp0.QE_read_acsii_index("%s/curr0/x0.export/index.xml" % wd)
                curr_tmp0.QE_read_acsii_wfc("%s/curr0/x0.export/wfc.1.xml" % wd )

                next_tmp0.QE_read_acsii_index("%s/next0/x0.export/index.xml" % wd)
                next_tmp0.QE_read_acsii_wfc("%s/next0/x0.export/wfc.1.xml" % wd )

                os.system("-rf %s/curr0/x0.export/wfc.1.xml" % wd)
                os.system("-rf %s/next0/x0.export/wfc.1.xml" % wd)


                curr_wfc0 = wfc(curr_tmp0,minband,nocc,curr_tmp0,nocc+1,maxband)
                next_wfc0 = wfc(next_tmp0,minband,nocc,next_tmp0,nocc+1,maxband)
                print "nac_method>=0: 2 WFCs are read and trimmed"



            if nac_method>=1:
                # In addition read info for N+1 electron wavefnctions

                # Convert binary files to xml - this is needed in newest version of QE
                # becuase pw_export will produce binary files in any case (bug?)
                os.system("%s convert %s/curr1/x1.export/wfc.1 %s/curr1/x1.export/wfc.1.xml" % (EXE_CONVERT,wd,wd))
                os.system("%s convert %s/next1/x1.export/wfc.1 %s/next1/x1.export/wfc.1.xml" % (EXE_CONVERT,wd,wd))

                curr_tmp1.QE_read_acsii_index("%s/curr1/x1.export/index.xml" % wd)
                curr_tmp1.QE_read_acsii_wfc("%s/curr1/x1.export/wfc.1" % wd )

                next_tmp1.QE_read_acsii_index("%s/next1/x1.export/index.xml" % wd)
                next_tmp1.QE_read_acsii_wfc("%s/next1/x1.export/wfc.1" % wd )

                os.system("-rf %s/curr1/x1.export/wfc.1.xml" % wd)
                os.system("-rf %s/next1/x1.export/wfc.1.xml" % wd)


                # Combine wavefunctions
                # Careful - we take the second wavefunction not from nocc+1 but from nocc+2
                # that is we simply skip the orbital nocc+1 from the second (charged) wavefunction
                # this is because in charged system this is a HOMO, while we need LUMO
                curr_wfc1 = wfc(curr_tmp0,minband,nocc,curr_tmp1,nocc+2,maxband)
                next_wfc1 = wfc(next_tmp0,minband,nocc,next_tmp1,nocc+2,maxband)
                print "nac_method>=1: 2 WFCs are read and trimmed"



            #-----------------------------------------------------------------
            if nac_method>=0:              
                # Options for wfc preprocessing
                if wfc_preprocess=="normalize": # just normalize
                    curr_wfc0.normalize()
                    next_wfc0.normalize()

                elif wfc_preprocess=="complete": # complete wfc with complex conjugate part, 
                                                 # the result is then normalized
                    curr_wfc0.complete()
                    next_wfc0.complete()

                elif wfc_preprocess=="restore":  # restore real wfc from the projection, the result is normailzed
                                                 # optionally the wfc can be completed with complex conjugate
                                                 # before restoring                    
                    curr_wfc0.restore(0,do_complete) # first 0 - k-point, second 1 - do complete wfc
                    next_wfc0.restore(0,do_complete)

                print "nac_method>=0: WFCs are preconditioned"

                # Finally compute Hamiltonian and the overlap matrix
                # In this case - any reasonable value for nocc leads to the same results
                # Keep in mind, the curr_wfc0 - is already only a subset of the whole wfc that has been 
                # computed, so all bands - from 0 to maxband - minband will be active!
                ham(curr_wfc0,next_wfc0,0,0, 0,maxband-minband, dt,"%s/0_Ham_%d" % (rd, curr_index) )
                print "nac_method>=0: Hamiltonian is computed and printed"
 
                if compute_Hprime==1:
                    os.system("%s convert %s/curr0/x0.export/grid.1 %s/curr0/x0.export/grid.1.xml" % (EXE_CONVERT,wd,wd))
                    curr_wfc0.QE_read_acsii_grid("%s/curr0/x0.export/grid.1.xml" % wd)
                    curr_wfc0.compute_Hprime(0,maxband-minband,"%s/0_Hprime_%d" % (rd, curr_index) )
                    os.system("-rf %s/curr0/x0.export/grid.1.xml" % wd)


            if nac_method>=1:
                # Options for wfc preprocessing
                if wfc_preprocess=="normalize": # just normalize
                    curr_wfc1.normalize()
                    next_wfc1.normalize()

                elif wfc_preprocess=="complete": # complete wfc with complex conjugate part,
                                                 # the result is then normalized
                    curr_wfc1.complete()
                    next_wfc1.complete()

                elif wfc_preprocess=="restore":  # restore real wfc from the projection, the result is normailzed
                                                 # optionally the wfc can be completed with complex conjugate
                                                 # before restoring
                    curr_wfc1.restore(0,do_complete) # first 0 - k-point, second 1 - do complete wfc
                    next_wfc1.restore(0,do_complete)

                print "nac_method>=1: WFCs are preconditioned"


                # Finally compute Hamiltonian and the overlap matrix
                # WORKS ONLY IN CASE: curr_wfc0.nspin==2 and curr_wfc1.nspin==2:
                # other cases yet to be implemented
                # here we use maxband-1 because 1 orbitals is skipped

                ham(curr_wfc1,next_wfc1,0,0, 0,maxband-1-minband, dt,"%s/1_Ham_%d" % (rd, curr_index) )
                print "nac_method>=1: Hamiltonian is computed and printed"


                if compute_Hprime==1:
                    os.system("%s convert %s/curr1/x0.export/grid.1 %s/curr1/x0.export/grid.1.xml" % (EXE_CONVERT,wd,wd))
                    curr_wfc1.QE_read_acsii_grid("%s/curr1/x1.export/grid.1" % wd)
                    curr_wfc1.compute_Hprime(0,maxband-1-minband,"%s/1_Hprime_%d" % (rd, curr_index) )
                    os.system("-rf %s/curr1/x1.export/grid.1.xml" % wd)



            #-----------------------------------------------------------------

            # Remove current run, make next run to be the current one
            os.system("rm -rf %s/curr0" % wd )
            os.system("mv %s/next0 %s/curr0" % (wd, wd) )
            if nac_method==1:
                os.system("rm -rf %s/curr1" % wd )
                os.system("mv %s/next1 %s/curr1" % (wd, wd) )

            print "old files deleted, new have become old"


# ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
# after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        print "End of step t=", t
        t = t + 1

#================= End of runMD function =============================
