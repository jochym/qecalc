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

from qecalc.multiphononcalc import MultiPhononCalc
from qecalc.pwcalc import PWCalc
from qecalc.singlephononcalc import SinglePhononCalc
from qecalc.d3calc import D3Calc

import os

def testCalc(calc):
    calc.launch()
    for task in calc.taskList:
        print task.output.properties()

if __name__ == "__main__":

    pwcalc = PWCalc()
    sphon = SinglePhononCalc(sectionList = ['pw.x', 'ph.x', 'dynmat.x'] )
    mphon = MultiPhononCalc(sectionList = ['pw.x', 'ph.x multi', 'q2r.x', 'matdyn.x'])
    d3calc = D3Calc()

    calcList = [pwcalc, sphon]

    for calc in calcList:
        testCalc(calc)


    os.system('rm -rf ' + pwcalc.pw.setting.get('outdir'))

    for task in mphon.taskList:
        task.syncSetting()
        task.output.parse()

    print mphon.lookupProperty('multi phonon')
    print mphon.matdyn.setting.get('flvec')
    mphon.dispersion.launch('M', 'Gamma', 'A', 'L', 50, 50, 50)
    mphon.dispersion.plot()

    print "Test completed successfully";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 18, 2009 3:17:57 PM$"
