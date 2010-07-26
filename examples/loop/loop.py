#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from qecalc.qetask.pwtask import PWTask

configString = """
[Launcher]
#paraPrefix:   mpiexec -n 2
#paraPostfix: -npool 2

outdir: temp/
"""

if __name__ == "__main__":
    pw = PWTask(configString = configString)
    #parse inputs and sync with Settings:
    pw.syncSetting()
    lat_params = [7.5, 7.6, 7.7]
    for a in lat_params:
        # whole lattice and structure will be auto updated on change in 'a' according
        # to the lattice symmetry (ibrav):
        pw.input.structure.lattice.a = a
        pw.input.save()
        pw.launch()
        print 'Stress = ', pw.output.property('stress')

    ecut_wfc_list = [15.5, 16, 17.5]
    for ecut_wfc in ecut_wfc_list:
        # if the variable did not exist, it will be created, othervise overwritten
        pw.input.namelist('system').add('ecutwfc', ecut_wfc)
        pw.input.save()
        pw.launch()
        print 'Total Energy = ', pw.output.property('total energy')

    import os
    os.system('rm -rf scf.out temp/')

