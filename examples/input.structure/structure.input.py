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

from qecalc.calcpw import PWCalc

if __name__ == '__main__':

    calc = PWCalc('config.ini')

    calc.pw.input.parse()
    print calc.pw.input.structure.lattice.toString()
    calc.pw.input.structure.lattice.ibrav = 8 # make orthorombic
    calc.pw.input.structure.lattice.a = 3
    calc.pw.input.structure.lattice.b = 4
    calc.pw.input.structure.lattice.c = 5

    print
    print calc.pw.input.structure.lattice.toString()

    print '\nLattice parameters: ', calc.pw.input.structure.lattice.latticeParams()

    #Write/update lattice into PW config file, if the file does not exist,
    # a new one will be created
    calc.pw.input.structure.lattice.save('./scf_2.in')

    #Write/update structure into PW config file, if the file does not exist,
    # a new one will be created
    print '\nPrinting structure:'
    print calc.pw.input.structure.diffpy()

    calc.pw.input.structure.save('./scf_44.in')



