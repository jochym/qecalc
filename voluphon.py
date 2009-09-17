__author__="markovsk"
__date__ ="$Sep 16, 2009 2:04:08 PM$"

import numpy
import volufit
from qephon import QEPhon, QEPhonQHA

# This class should be initialized with matdyn.modes from different
# volume expansios and then it is fed to QHA program by Eyvaz Isaev
# to obtain total DOSes
class VoluPhon():
    def __init__(self, fname, prcntVolume):
        #will  need this unless QE has improved
        self.phonQHA = QEPhonQHA(fname)
        self.phon = QEPhon(fname)
        self.__prcntVolume = prcntVolume

        self.__latParFitDirective = ('polynom', 1)
        self.__energyFitDirective = ('polynom', 4)
        self.__freqFitDirective = ('polynom', 2)

    def setA(self, values):
        self.a = volufit.ValueFit(self.__prcntVolume, values, \
                                                    *self.__latParFitDirective)

        
    def setC(self, values):
        self.c = volufit.ValueFit(self.__prcntVolume, values, \
                                                    *self.__latParFitDirective)


    def setEnergy(self, values):
        self.energy = volufit.ValueFit(self.__prcntVolume, values, \
                                                    *self.__energyFitDirective)


    def setPhonons(self, indexRange):
        """Will read freqs from x_matdyn.modes files"""
        fname = str(indexRange[0]) + '_' + self.phon.matdynModes
        Pol, Omega, qPoints = self.phon.getMultiPhonon(fname)
        volOmega = numpy.zeros(shape=(len(indexRange), numpy.shape(Omega)[0], \
                                                      numpy.shape(Omega)[1]  ) )
        volOmega[0] = Omega
        for i in range(1,len(indexRange)):
            fname = str(indexRange[i]) + '_' + self.phon.matdynModes
            Pol, Omega, qPoints = self.phon.getMultiPhonon(fname)
            volOmega[i] = Omega

        self.freqs = volufit.FreqFit(self.__prcntVolume, volOmega, \
                                                     *self.__freqFitDirective)


    def gammaDispersion(self, *pathNPoints):
        if self.freqs.fitter.type != 'polynom':
            raise Exception('This method is only relevant for polynomial fit')
        self.phon.dispersion.setPath(*pathNPoints)
        self.phon.dispersion.setValues(-self.freqs.coeff()[:,:,-1])
        self.phon.dispersion.plot()


if __name__ == "__main__":
    print "Hello World";
