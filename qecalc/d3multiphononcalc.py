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
from qecalc.d3calc import D3Calc
from qeutils import kmesh


qgrid = [2,2,2]

mphon = MultiPhononCalc('config.ini')

mphon.ph.input.parse()
mphon.ph.input.qpoints.setAutomatic(qgrid)
mphon.ph.input.save()

mphon.launch([mphon.pwph, mphon.q2r])

#irred_qpoints = mphon.ph.output.property('qpoints')

grid, qpoints_indep, qpoints_full = mphon.ph.output.property('qpoints')

print qpoints_indep
print qpoints_full


#generate qpoint grid:


#[nq1,nq2,nq3], qpoints_indep, qpoints_full

#qpoints = kmesh.kMeshCart([2,2,2],mphon.pw.input.structure.lattice.reciprocalBase())

mphon.matdyn.input.parse()
mphon.matdyn.input.qpoints.set(qpoints_full)
mphon.matdyn.input.save()
mphon.matdyn.launch()

modes, freqs, qpoints =  mphon.matdyn.output.property('multi phonon')

print modes, freqs, qpoints

d3calc = D3Calc('config.ini')

d3data = []
tasks = [d3calc.ph, d3calc.d3]
for qpoint in qpoints_indep:
    for task in tasks:
        task.input.parse()
        task.input.qpoints.set(qpoint)
        task.input.save()
        task.launch()

    print qpoint
    print d3calc.d3.output.property('d3 tensor')

    d3data.append(d3calc.d3.output.property('d3 tensor'))

# print d3data


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 8, 2009 11:20:41 AM$"
