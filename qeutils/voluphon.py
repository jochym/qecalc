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
import numpy
import volufit
#from qephon import QEPhon, QEPhonQHA

from qecalc.multiphononcalc import MultiPhononCalc

# This class should be initialized with matdyn.modes from different
# volume expansios and then it is fed to QHA program by Eyvaz Isaev
# to obtain total DOSes
# a,b,c,energy - fitting objects. fitted value is accesible through fittedValue
class VoluPhon():
    def __init__(self, fname, prcntVolume):
        #will  need this unless QE has improved
        #self.mphonQHA = QEPhonQHA(fname)
        self.mphon = MultiPhononCalc(fname)
        self.__prcntVolume = prcntVolume


    def setA(self, values, fitter):
        self.a = volufit.ValueFit(self.__prcntVolume, values, fitter)

        
    def setC(self, values, fitter):
        self.c = volufit.ValueFit(self.__prcntVolume, values, fitter)


    def setEnergy(self, values, fitter):
        self.energy = volufit.ValueFit(self.__prcntVolume, values, fitter)


    def setPhonons(self, indexRange, fitter):
        """Will read freqs from x_matdyn.modes files and fit the freqs"""
        matdynModesName = self.mphon.matdyn.setting.matdynModes
        self.mphon.matdyn.setting.matdynModes = str(indexRange[0]) + '_' + matdynModesName
        self.mphon.matdyn.output.parse()
        Pol, Omega, qPoints = self.mphon.matdyn.output.property('multi phonon')
        volOmega = numpy.zeros(shape=(len(indexRange), numpy.shape(Omega)[0], \
                                                      numpy.shape(Omega)[1]  ) )
        volOmega[0] = Omega
        for i in range(1,len(indexRange)):
            self.mphon.matdyn.setting.matdynModes = str(indexRange[i]) + '_' + matdynModesName
            self.mphon.matdyn.output.parse()
            Pol, Omega, qPoints = self.mphon.matdyn.output.property('multi phonon')
            volOmega[i] = Omega
        self.mphon.matdyn.setting.matdynModes = matdynModesName
        self.freqs = volufit.FreqFit(self.__prcntVolume, volOmega,fitter)


    def gammaDispersion(self, *pathNPoints):
        if self.freqs.fitter.type != 'polynom':
            raise Exception('This method is only relevant for polynomial fit')
        self.mphon.dispersion.setPath(*pathNPoints)
        self.mphon.dispersion.setValues(-self.freqs.coeff()[:,:,-1])
        self.mphon.dispersion.plot()


    def prcntVolume (self):
        return self.__prcntVolume


if __name__ == "__main__":
    indexRange = range(0,7,2)
    prcntVol = array(indexRange)/1000.0
    voluPhon = VoluPhon('config.ini', prcntVol)
    voluPhon.setPhonons(indexRange)
    voluPhon.gammaDispersion('Gamma','K', 'M', 'Gamma', 'A', 'H', 'L', 'A', \
                              100, 100, 100, 100, 100, 100, 100)

__author__="markovsk"
__date__ ="$Sep 16, 2009 2:04:08 PM$"