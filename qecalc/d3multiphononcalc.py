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
# from qeutils import kmesh


qgrid = [2,2,2]
fildrho = "'si.drho'"

mphon = MultiPhononCalc('config.ini')

mphon.ph.input.parse()
mphon.ph.input.namelist('inputph').remove('fildrho')
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

# calculation at Gamma:
d3calc.ph.input.parse()
d3calc.d3.input.parse()
#fildrho = d3calc.ph.input.namelist('inputph').param('fildrho')
fildrhoG = "'" + fildrho.strip("'") + '_G' + "'"
d3calc.ph.input.namelist('inputph').set('fildrho', fildrhoG)
d3calc.d3.input.namelist('inputph').set('fild0rho', fildrhoG)
d3calc.d3.input.namelist('inputph').set('fildrho', fildrhoG)
d3calc.ph.input.save()
d3calc.d3.input.save()
d3calc.pw.launch()
for task in tasks:
    task.input.parse()
    task.input.qpoints.set([0.0, 0.0, 0.0])
    task.input.save()
    task.launch()
print '[0,0,0]'
print d3calc.d3.output.property('d3 tensor')    

# non Gamma points    
d3calc.ph.input.namelist('inputph').set('fildrho', fildrho)
d3calc.d3.input.namelist('inputph').set('fildrho', fildrho)
d3calc.ph.input.save()
d3calc.d3.input.save()
for qpoint in qpoints_indep[1:]:
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
