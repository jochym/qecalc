#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



from qecalc.qetask.pwtask import PWTask
from qeutils import kmesh
from qeutils import bandstobxsf

if __name__ == "__main__":

    
    pw = PWTask('config.ini')
    pw.input.parse()
    pw.output.parse()
    print pw.output.property('fermi energy')

    kpoints = [17,17,17]

    bands = kmesh.packBands(kpoints,pw.output.property('bands')[1])

    bandstobxsf(7.8071, \
                pw.input.structure.lattice.reciprocalBase(), \
                bands, 'test.bxsf')


__author__="Nikolay Markovskiy"
__date__ ="$Dec 4, 2009 7:52:03 PM$"
