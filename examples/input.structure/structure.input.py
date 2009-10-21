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

#from qecalc.calcpw import PWCalc
from qecalc.qetask.taskpw import PWTask
from qecalc.setting import Setting

from qecalc.qetask.qeparse.pwinput import PWInput

if __name__ == '__main__':

    #calc = PWCalc('config.ini')
    setting = Setting('config.ini')
    pwInput = PWInput(setting)
    pwInput.parse()
#    pw = PWTask(setting)
#    pw.input.parse()
#    calc.pw.input.parse()
#    calc.pw.input.structure.lattice.a = 13
    calc.pw.input.structure.lattice.b = 14
    calc.pw.input.structure.lattice.c = 3

    print calc.pw.input.structure.lattice.latticeParams()

    calc.pw.input.structure.lattice.saveLatticeToPWSCF('./scf_2.in')

    print calc.pw.input.structure.structure

    calc.pw.input.structure.saveStructureToPWSCF('./scf_3.in')



