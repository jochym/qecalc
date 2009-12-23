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
paraPrefix:   mpiexec -n 8
paraPostfix: -npool 8

serialPrefix: mpiexec -n 1
serialPostfix:

#useTorque: True # default: False
#paraPrefix: mpirun --mca btl openib,sm,self
#paraPostfix: -npool 900
#
#serialPrefix: mpiexec
#serialPostfix:

#Name of a script to execute a command on multiple nodes
#relevant if outdir is not located on Parallel/Network File system.
#Default value is empty
#paraRemoteShell: bpsh -a

# this string will be passed to qsub, -d workingDir -V are already there:
paraTorqueParams: -l nodes=2:ppn=12 -N myjob -j oe
serialTorqueParams: -l nodes=1:ppn=1 -N myjob -j oe

#outdir: /scratch/temp


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
    d3calc = D3Calc(configString = configString)

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
