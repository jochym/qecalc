#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from diffpy.Structure.lattice import Lattice, cosd
from math import sqrt, degrees
import numpy
from qecalc.qetask.qeparser.qeinput import QEInput
from qecalc.qetask.qeparser.qelattice import QELattice
from qecalc.qetask.qeparser.qestructure import QEStructure

if __name__ == '__main__':

    pwInput = QEInput('scf.in', type = 'pw')
    pwInput.parse()
    qeLattice = QELattice(qeConf = pwInput)
    print qeLattice.latticeParams()
    qeLattice.a = 13.0
    qeLattice.b = 24.0
    qeLattice.c = 3.
    qeLattice.ibrav = 4
    print qeLattice.b, qeLattice.c,
    print qeLattice.latticeParams()
    qeLattice.saveLatticeToPWSCF('./scf_2.in')

    pwInput = QEInput('scf.in', type = 'pw')
    pwInput.parse()
    myStruct = QEStructure(qeConf = pwInput)
    myStruct.lattice.ibrav = 4
    print myStruct.lattice.a
    print myStruct.lattice.c
    myStruct.lattice.a = 43
    #pwInput.save()
    myStruct.saveStructureToPWSCF('scf_3.in')
    print myStruct.structure

#    qeLattice2 = QELattice()
#    qeLattice2.setLatticeFromPrimitiveVectors(qeLattice.ibrav, qeLattice.lattice().base )
    #print qeLattice.lattice().base
    #testLattice = Lattice(5,5,5,90,90,90)
    #print testLattice.base
