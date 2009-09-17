# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="markovsk"
__date__ ="$Sep 11, 2009 12:42:48 PM$"

from scipy import *
from scipy.optimize import *
from pylab import *
#from dos_utils import *
#from matdyn import *
from qecalc import QECalc

class QuasiHarmonic():
    def __init__(self, fname, indexRange, prcntVolume, modePrefix):
        """indexRange is for files' indexing. prcntVolume - coresponding volume
       expansions"""
        self.qecalc = QECalc(fname)
        self.__modePrefix = modePrefix
        self.__indexRange = indexRange
        self.__prcntVolume = prcntVolume
        if len(indexRange) != len(prcntVolume):
            raise Exception("QuasiHarmonic dimensions: indexRange != prcntVolume")
        self.loadPhonons()
        self.__polyOrder = 1
        self.fitMatdyn()

    def fitFunc(self, p, x):
        return p[0] * x
#+ p[1] * x * x
        #+ p[2] * x * x * x + p[3] * x * x * x * x

    def fitOmega(self, prcntOmega):
    #    fitfunc = lambda p, x: p[0] * x + p[1] * x * x
        errFunc = lambda p, x, y: (y - self.fitFunc(p, x))
      # initial parameters:
        pinit = [prcntOmega[1]/self.__prcntVolume[1]]
#        pinit = [1,0]
        v, success = leastsq(errFunc, pinit, args=(self.__prcntVolume, prcntOmega))
        if success < 1 or success > 4:
            print success
            print prcntOmega
            print v
    #   assert success != 1, "fitOmega: Fitting was not successful"
        return v


    def getQuasiFreqs(prcntVol):
        """Returns array of freqs, corresponding to the provided volume expansion"""
        nPoints = shape(self.__coeffOmega)[0]
        nModes = shape(self.__coeffOmega)[1]
        fittedFreqs = zeros(shape = (nPoints, nModes) )
        for i in range( nPoints ):
            for j in range( nModes ):
                p = self.__coeffOmega[i,j]
                fittedFreqs[i,j] = p[0]*prcntVol
        return (fittedFreqs+1.0)*self.__omega0

    def getFittedFreqs(prcntVol):
        nPoints = shape(self.__coeffOmega)[0]
        nModes = shape(self.__coeffOmega)[1]
        fittedFreqs = zeros(shape = (nPoints, nModes) )
        for i in range( nPoints ):
            for j in range( nModes ):
                p = self.__coeffOmega[i,j]
                fittedFreqs[i,j] = fitFunc(p, prcntVol)
        return (fittedFreqs+1.0)*self.__omega0

    def loadPhonons(self):
        fname = self.__modePrefix + str(self.__indexRange[0]) + '.modes'
        Pol, Omega, qPoints = self.qecalc.getMultiPhonon(fname)
        self.__volOmega = zeros(shape=(len(self.__indexRange), shape(Omega)[0], shape(Omega)[1]  ) )
        self.__volPol = zeros(shape=(len(self.__indexRange), shape(Pol)[0], shape(Pol)[1], shape(Pol)[2], shape(Pol)[3] ) )
        self.__volPol[0] = Pol
        self.__volOmega[0] = Omega
        for i in range(1,len(indexRange)):
            fname = self.__modePrefix + str(self.__indexRange[i]) + '.modes'
            Pol, Omega, qPoints = self.qecalc.getMultiPhonon(fname)
            self.__volPol[i] = Pol
            self.__volOmega[i] = Omega
#        return volPol, volOmega, qPoints

    def fitMatdyn(self):
    #    volPol, volOmega, qPoints = loadPhonons(indexRange, prefix)
    #    prcntVol = array(indexRange)/1000.0
        # percent change in Omegas relative to equilibrium
        self.__volPrcntOmega = zeros(shape(self.__volOmega))
        # introduce Omega0 to get rid of 0s in denominator
        Omega0 = copy(self.__volOmega[0])
        self.__omega0 = Omega0
        for i in range(shape(Omega0)[0]):
            for j in range(shape(Omega0)[1]):
                if Omega0[i,j] == 0.:
                    Omega0[i, j] = 1
        for i in range(len(self.__indexRange)):
            self.__volPrcntOmega[i] = (self.__volOmega[i] - self.__volOmega[0])/Omega0
        self.__coeffOmega = zeros( shape = ( shape(Omega0)[0], shape(Omega0)[1], self.__polyOrder) )
        nonlinearity = 0
        for i in range( shape(Omega0)[0] ):
            for j in range( shape(Omega0)[1] ):
                self.__coeffOmega[i,j] = self.fitOmega(self.__volPrcntOmega[:,i,j])
                if self.__polyOrder > 1:
                    nonlinearity = nonlinearity + abs(self.__coeffOmega[i,j][1])
#    print volPrcntOmega[:,1000,3]
#    print fitOmega(prcntVol, volPrcntOmega[:,1000,3])
#    plot(volPrcntOmega[:,1000,3])
#    show()

#    return prcntVol, coeffOmega, volOmega, nonlinearity

    def gammaDispersion(self, *pathNPoints):
        self.qecalc.dispersion.setPath(*pathNPoints)
        self.qecalc.dispersion.setValues(-self.__coeffOmega[:,:,0])
        self.qecalc.dispersion.plot()


if __name__ == "__main__":
    indexRange = range(0,15,2)
    mphon = QEPhon('config.ini')
    mphon.dispersion.setPhononPath('Gamma','K', 'M', 'Gamma', 'A', 'H', 'L', 'A', \
                              100, 100, 100, 100, 100, 100, 100)
    mphon.dispersion.plot()
    prcntVol = array(indexRange)/1000.0
    qh = QuasiHarmonic('config.ini', indexRange, prcntVol, 'matdyn_')
    qh.gammaDispersion('Gamma','K', 'M', 'Gamma', 'A', 'H', 'L', 'A', \
                              100, 100, 100, 100, 100, 100, 100)


if __name__ == "__main__":
    print "Hello World";
