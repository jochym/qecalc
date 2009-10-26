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
from thermo.thermodyn import PhononThermodynamics
import numpy
import os
from qeutils.volufit import ValueFit
from qeutils.voluphon import VoluPhon

# AlB2:
a_range = numpy.array([5.6724469859199997, 5.6750274863032173, 5.6785814939910439, 5.6815649668271124, 5.6845383973710142, 5.6876674822014746, 5.6907180972966769, 5.6936814257211736, 5.6968461022751251, 5.6997600959854884, 5.7027596495990052, 5.7057350314770057, 5.7087138263269255, 5.7116526400046892, 5.7146061023692676, 5.7175703988694977, 5.7205909099531231, 5.723237255174543, 5.7264457055705025, 5.7292626217675133, 5.7321587312193687, 5.7350400741208398])
c_a_range = numpy.array([1.0895502669999999, 1.0902407788521165, 1.0903670852940217, 1.090818910093668, 1.091273295187184, 1.0916348386375583, 1.0920383848325765, 1.0924889987431634, 1.0928205534256967, 1.0932932254548433, 1.09371350862399, 1.0941445604290785, 1.0945705200579436, 1.0950163576406637, 1.0954506600760952, 1.0958756352290124, 1.0962652071369421, 1.0968668161284669, 1.0971422049606323, 1.0976396117859046, 1.0980884501524621, 1.0985427361942841])
e_total = numpy.array([-17.173580479999998, -17.173578419999998, -17.173568249999999, -17.17354963, -17.17352326, -17.17348866, -17.173445699999998, -17.17339467, -17.17333614, -17.17326946, -17.17319517, -17.173113239999999, -17.173023669999999, -17.172927080000001, -17.1728223, -17.172710599999998, -17.172591560000001, -17.17246548, -17.172331499999999, -17.172185509999998, -17.172038229999998, -17.17188402])

fc_name = 'alb2888.fc'
matdynfldos = 'alb2888.phdos'

indexRange = [0,2,4,6,8,10,12]
prcntVol = 2.0*numpy.array(indexRange)/1000.0

# Volume range after fitting of energies
nVolumePoints = 1000.0
finePrcntVolMax = 42.0/1000.0
finePrcntVol = numpy.linspace(prcntVol[0], finePrcntVolMax, nVolumePoints)

temperatures = [0, 1, 10, 50,100,150,200,250,300,450,625,785]

class FreeEnergy():
# total free energy class at some temperature
    def __init__(self, prcntVol, e_total, phonon_energy):
        self.prcntVol = prcntVol
        self.fitEnergy = ValueFit(prcntVol,e_total[indexRange], 'polynom', 4, True)
        self.fitPhononEnergy = ValueFit(prcntVol,phonon, 'polynom', 1, True)
        
    def totalFreeEnergy(self,v, *params):
        return self.fitEnergy.fittedValue(v) + self.fitPhononEnergy.fittedValue(v)

    def minimizeTotalFreeEnergy(self,minVol, maxVol):
        import scipy.optimize
        brentOut=scipy.optimize.brent(self.totalFreeEnergy,(),
           (minVol, maxVol),\
                          tol = 1.e-5, full_output = 1)
#        print '\n\nBrentOut:'
#        print brentOut
        optPrcntVolume = brentOut[0]
        optFreeEnergy = brentOut[1]
        return optPrcntVolume




if __name__ == "__main__":
    mphonCalc = MultiPhononCalc("config.ini")
    voluPhon = VoluPhon("config.ini", prcntVol)
    voluPhon.setA(a_range[indexRange], 'polynom', 1, True)
    c_range = a_range*c_a_range
    voluPhon.setC(c_range[indexRange], 'polynom', 1, True)
    # Generate phonon doses from different fc files
    for i in indexRange:
        os.system('cp ' + str(i) + '_' + fc_name + ' ' + fc_name)
        mphonCalc.matdyn.launch()
        os.system('cp ' + matdynfldos + ' ' + str(i) + '_' + matdynfldos )

    #plot free energies:
    for temperature in temperatures:
        lines = ''
        phonon = []
        for i,v in zip(indexRange, prcntVol):
            axis, dos = qecalc.getPhononDOS(str(i) + '_' + matdynfldos)
            phonTherm = PhononThermodynamics(axis, dos)
            Cv = phonTherm.Cv(temperature)
            phonon.append(phonTherm.freeEnergy(temperature))
            totalFreeEnergy = e_total[i] + phonon[-1]
            lines = lines + '%f    %f    %f    %f    %f\n'%(v, \
            totalFreeEnergy, e_total[i], phonon[-1], Cv)
        file = open(str(temperature) + '_free_energy.out', 'w')
        file.write(lines)
        file.close()

        freeEnergy = FreeEnergy(prcntVol, e_total, phonon)
        optPrcntVolume = freeEnergy.minimizeTotalFreeEnergy(finePrcntVol[0], finePrcntVol[-1])
        print "\nOptimized lattice parameters:"
        a = voluPhon.a.fittedValue(optPrcntVolume)
        c = voluPhon.c.fittedValue(optPrcntVolume)
        print "Temperature = %fK"%temperature
        print 'Optimized %Volume = ', optPrcntVolume*100.0
        print 'Optimized a = ', a
        print 'Optimized c = ', c
        print '\n'
        lines = ''
        for v in finePrcntVol:
            e = freeEnergy.fitEnergy.fittedValue(v)
            p = freeEnergy.fitPhononEnergy.fittedValue(v)
            totalFreeEnergy = e + p
            lines = lines + '%f    %f    %f    %f\n'%(v, \
            totalFreeEnergy, e, p)
        file = open(str(temperature) + '_fitted_free_energy.out', 'w')
        file.write(lines)
        file.close()




__author__="Nikolay Markovskiy"
__date__ ="$Oct 12, 2009 3:01:51 PM$"
