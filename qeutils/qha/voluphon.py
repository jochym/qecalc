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

from qecalc.qetask import PWTask
from qecalc.qetask import MatdynTask
from qeutils.phdispersion import PHDispersion

class VoluPhon():
    def __init__(self, prcntVolume, filename = None ):
        """
        This class requires a set of  "n_matdyn.modes" files 
        corresponding to different volume expansions 
        prcntVolume - array with 
        """
        self.pw = PWTask(filename = filename)
        self.matdyn = MatdynTask(filename = filename)
        self.__prcntVolume = prcntVolume


    def setA(self, values, fitter):
        self.a = volufit.ValueFit(self.__prcntVolume, values, fitter)

        
    def setC(self, values, fitter):
        self.c = volufit.ValueFit(self.__prcntVolume, values, fitter)


    def setEnergy(self, values, fitter):
        self.energy = volufit.ValueFit(self.__prcntVolume, values, fitter)


    def setPhonons(self, indexRange, fitter):
        """Will read freqs from x_matdyn.modes files and fit the freqs"""
        matdynModesName = self.matdyn.setting.get('flvec')
        self.matdyn.setting.set('flvec', str(indexRange[0]) + '_' + matdynModesName)
        self.matdyn.output.parse()
        Pol, Omega, qPoints = self.matdyn.output.property('multi phonon')
        volOmega = numpy.zeros(shape=(len(indexRange), numpy.shape(Omega)[0], \
                                                      numpy.shape(Omega)[1]  ) )
        volOmega[0] = Omega
        for i in range(1,len(indexRange)):
            self.matdyn.setting.set('flvec', str(indexRange[i]) + '_' + matdynModesName)
            self.matdyn.output.parse()
            Pol, Omega, qPoints = self.matdyn.output.property('multi phonon')
            volOmega[i] = Omega
        self.matdyn.setting.set('flvec', matdynModesName)
        self.freqs = volufit.FreqFit(self.__prcntVolume, volOmega,fitter)


    def gammaDispersion(self, *pathNPoints):        
        #self.pw.input.parse()
        #self.matdyn.input.parse()
        self.dispersion = PHDispersion( self.pw.input.structure.lattice, self.matdyn)
        if self.freqs.fitter.type != 'polynom':
            raise Exception('This method is only relevant for polynomial fit')
        self.dispersion = PHDispersion( self.pw.input.structure.lattice, self.matdyn)
        self.dispersion.setPath(*pathNPoints)
        self.dispersion.setValues(-self.freqs.coeff()[:,:,-1])
        self.dispersion.save('gamma_disp')
        self.dispersion.plot()


    def prcntVolume (self):
        return self.__prcntVolume


if __name__ == "__main__":
    indexRange = range(0,5,2)
    prcntVol = array(indexRange)/1000.0
    voluPhon = VoluPhon(prcntVolume = prcntVol, filename = 'config.ini')
    voluPhon.setPhonons(indexRange)
    voluPhon.gammaDispersion('Gamma','K', 'M', 'Gamma', 'A', 'H', 'L', 'A', \
                              100, 100, 100, 100, 100, 100, 100)

