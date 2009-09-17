# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="markovsk"
__date__ ="$Sep 16, 2009 12:42:48 PM$"

import numpy
import scipy.optimize

class VoluFit():
    def __init__(self, *fitDirective):
        if fitDirective[0] != 'polynom':
            raise Exception('This fit is not implemented')
        self.type = fitDirective[0]
        polyOrder = fitDirective[1]
        if polyOrder < 1 or polyOrder > 4:
            raise Exception('this polyOrder is not supported')
        self.order = self.__polyOrder = polyOrder


    def func(self, p, x):
        if self.__polyOrder == 1:
            return p[0]*x
        if self.__polyOrder == 2:
            return p[1]*x + p[0]*x*x
        if self.__polyOrder == 3:
            return p[2]*x + p[1]*x*x + p[0]*x*x*x
        if self.__polyOrder == 4:
            return p[3]*x + p[2]*x*x + p[1]*x*x*x + p[0]*x*x*x*x


    def fit(self, xdata, ydata):
        errFunc = lambda p, x, y: (y - self.func(p, x))
      # initial parameters:
        pinit = numpy.zeros(self.order)
        pinit[0] = (ydata[1]-ydata[0])/(xdata[1] - xdata[0])
        v, success = scipy.optimize.leastsq(errFunc, pinit, args=(xdata, ydata))
        if success < 1 or success > 4:
            print success
            print sef.__ydata
            print v
    #   assert success != 1, "fitFreq: Fitting was not successful"
        return v



class FreqFit(VoluFit):
    def __init__(self, prcntVolume, freqs, *fitDirective):
        """indexRange is for files' indexing. prcntVolume - coresponding volume
        expansions"""
        self.fitter = VoluFit(*fitDirective)
#        Fit.__init__(self, *fitDirective)
        self.__prcntVolume = prcntVolume
        self.__freqs = freqs
        self.__fit()
        

    def __fit(self):
        import copy
        self.__prcntFreqs = numpy.zeros(self.__freqs.shape)
        # introduce Omega0 to get rid of 0s in denominator
        Omega0 = copy.copy(self.__freqs[0])
        self.__freqs0 = Omega0
        for i in range(Omega0.shape[0]):
            for j in range(Omega0.shape[1]):
                if Omega0[i,j] == 0.:
                    Omega0[i, j] = 1
        for i in range(self.__freqs.shape[0]):
            self.__prcntFreqs[i] = (self.__freqs[i] - self.__freqs[0])/Omega0
        self.__coeff = numpy.zeros( shape = ( Omega0.shape[0], Omega0.shape[1],\
                                              self.fitter.order) )
#        nonlinearity = 0
        for i in range( Omega0.shape[0] ):
            for j in range( Omega0.shape[1] ):
                self.__coeff[i,j] = self.fitter.fit(self.__prcntVolume, \
                                                 self.__prcntFreqs[:,i,j])
#                if self.fitter.__polyOrder > 1:
#                    nonlinearity = nonlinearity + abs(self.__coeff[i,j][1])


    def quasiFreqs(prcntVol):
        """Returns array of freqs, corresponding to the provided volume expansion"""
        if self.fitter.type != 'polynom':
            raise Exception('Does not have any meaning for the given type of fit')
        nPoints = shape(self.__coeff)[0]
        nModes = shape(self.__coeff)[1]
        fittedFreqs = zeros(shape = (nPoints, nModes) )
        for i in range( nPoints ):
            for j in range( nModes ):
                p = self.__coeff[i,j]
                fittedFreqs[i,j] = p[-1]*prcntVol
        return (fittedFreqs+1.0)*self.__freqs0


    def fittedFreqs(prcntVol):
        """Returns array of freqs, corresponding to the provided volume expansion"""
        nPoints = shape(self.__coeff)[0]
        nModes = shape(self.__coeff)[1]
        fittedFreqs = zeros(shape = (nPoints, nModes) )
        for i in range( nPoints ):
            for j in range( nModes ):
                p = self.__coeff[i,j]
                fittedFreqs[i,j] = self.fitter.func(p, prcntVol)
        return (fittedFreqs+1.0)*self.__freqs0


    def freqs(self):
        """Returns volume expansion"""
        return self.__freqs


    def coeff(self):
        return self.__coeff



class ValueFit():
    def __init__(self, prcntVolume, values, *fitDirective):
        """indexRange is for files' indexing. prcntVolume - coresponding volume
        expansions"""
        self.fitter = VoluFit(*fitDirective)
#        Fit.__init__(self, *fitDirective)
        self.__prcntVolume = numpy.array(prcntVolume)
        self.__values = numpy.array(values)
        self.__fit()


    def __fit(self):
        if self.__values[0] == 0:
            raise Exception('Denominator is 0!')
        self.__prcntValues = (self.__values - self.__values[0])/self.__values[0]
        self.__coeff = self.fitter.fit(self.__prcntVolume, \
                                            self.self.__prcntValues)


    def fittedValue(prcntVol):
        return (self.fitter.func(self.__coeff, prcntVol) + 1.0)*self.__values[0]


    def values():
        return(self.__values)


    def coeff(self):
        return self.__coeff


if __name__ == "__main__":
    print "Hello World";
