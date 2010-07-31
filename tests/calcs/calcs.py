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

configString = """
# all the relevant input files must be preconfiguered for specific tasks
# before using this class

[Launcher]
# parallelization parameters
# if this section is empty - serial mode is used
paraPrefix:   mpiexec -n 2
paraPostfix: -npool 2

serialPrefix: mpiexec -n 1
serialPostfix:


outdir: temp/


[pw.x]
# pwscf input/output files
pwfInput:  scf.in
pwOutput: scf.out


[ph.x]
#ph.x input/ouput, relevant to all phonon calculations:
phInput:  ph.in
phOutput: ph.out

[ph.x multi]
#ph.x input/ouput, relevant to all phonon calculations:
phInput:  ph_disp.in
phOutput: ph_disp.out

[ph.x d3]
#ph.x input/ouput, relevant to all phonon calculations:
phInput:  ph.in
phOutput: ph.out
fildrho: mgb2.drho_G


[dynmat.x]
#dynmat.x input/output files relevant to single phonon calculation
dynmatInput:  dynmat.in
dynmatOutput: dyn.out


[q2r.x]
# input/output files relevant to multiple phonon calculation
q2rInput:      q2r.in
q2rOutput:     q2r.out


[matdyn.x]
# input/output files relevant to multiple phonon calculation
matdynInput:   matdyn.in
matdynOutput:  matdyn.out

[d3.x]
d3Input:  d3.in
d3Output: d3.out
"""



def testCalc(calc):
    calc.launch()
    for task in calc.taskList:
        print task.output.properties()

if __name__ == "__main__":

    pwcalc = PWCalc(configString = configString)
    sphon = SinglePhononCalc(configString = configString, sectionList = ['pw.x', 'ph.x', 'dynmat.x'] )
    mphon = MultiPhononCalc(configString = configString, sectionList = ['pw.x', 'ph.x multi', 'q2r.x', 'matdyn.x'])
    d3calc = D3Calc(configString = configString, sectionList = ['pw.x', 'ph.x d3', 'd3.x'])

    calcList = [pwcalc, sphon, mphon, d3calc]

    try:
        for calc in calcList:
            testCalc(calc)
        mphon.dispersion.launch('M', 'Gamma', 'A', 'L', 50, 50, 50)
        mphon.dispersion.plot()

        print "Test completed successfully";

    except:
        print "Test did NOT complete successfully";
    os.system('rm -rf *.out *.freq *.modes *.dyn*  *.mold *.axsf *.fc temp/ *_G CRASH')

    

__author__="Nikolay Markovskiy"
__date__ ="$Dec 18, 2009 3:17:57 PM$"
