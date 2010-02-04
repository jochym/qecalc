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
paraPrefix:   mpiexec -n 8
paraPostfix: -npool 8

outdir: temp/
"""

if __name__ == "__main__":
    pw = PWTask(configString = configString)
    #parse inputs and sync with Settings:
    pw.syncSetting()
    lat_params = [5.5, 5.6, 5.7]
    for a in lat_params:
        # whole lattice and structure will be auto updated on change in a according
        # to lattice symmetry:
        pw.input.structure.lattice.a = a
        # changes in structure should be propagated into the parsing object:
        pw.input.structure.updatePWInput()
        pw.input.save()
        # or just use pw.input.structure.save()
        pw.launch()
        print 'Stress = ', pw.output.property('stress')

    ecut_wfc_list = [15, 16, 17.5]
    for ecut_wfc in ecut_wfc_list:
        # if the variable did not exist, it will be created, othervise overwritten
        pw.input.namelist('system').add('ecutwfc', ecut_wfc)
        pw.input.save()
        pw.launch()
        print 'Total Energy = ', pw.output.property('total energy')

