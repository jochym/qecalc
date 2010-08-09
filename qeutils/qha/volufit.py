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
import scipy.optimize

class Fit():
    def __init__(self):
        self._coeff = None

    def func(self, p, x):
        """
        p - coefficients, x - values
        """
        pass

    def fit(self): pass


class PolynomialFit(Fit):
    def __init__(self, *fitDirective):
        Fit.__init__(self)
        self.type = 'polynom'
        polyOrder = fitDirective[0]
        hasConstantTerm = fitDirective[1]
        if hasConstantTerm == True:
            self.constantTerm = 1
        else:
            self.constantTerm = 0

        if polyOrder < 1 or polyOrder > 4:
            raise Exception('this polyOrder is not supported')
        self._polyOrder = polyOrder
        self.order = polyOrder + 1

    def fit(self, xdata, ydata):
        errFunc = lambda p, x, y: (y - self.func(p, x))
      # initial parameters:
        pinit = numpy.zeros(self.order)
        pinit[-1] = (ydata[1]-ydata[0])/(xdata[1] - xdata[0])
        v, success = scipy.optimize.leastsq(errFunc, pinit, args=(xdata, ydata))
        if success < 1 or success > 4:
            print success
            print ydata
            print v
    #   assert success != 1, "fitFreq: Fitting was not successful"
        return v

    def func(self, p, x):
        if self._polyOrder == 1:
            return p[1]*self.constantTerm + p[0]*x
        if self._polyOrder == 2:
            return  p[2]*self.constantTerm + p[1]*x + p[0]*x*x
        if self._polyOrder == 3:
            return p[3]*self.constantTerm + p[2]*x + p[1]*x*x + p[0]*x*x*x
        if self._polyOrder == 4:
            return p[4]*self.constantTerm + p[3]*x + p[2]*x*x + p[1]*x*x*x + \
                                                              p[0]*x*x*x*x



class BirchMurnaghanFit(Fit):
    def __init__(self):
        Fit.__init__(self)
        
    def func(self, p, x):
        """
        x is assumed to be a percent volume,
        p[0] - E0
        p[1] - V0*B0
        p[2] - B0'
        """
        val = p[0] + 9.0/16.0*p[1]*( (1./(x+1.0)**(2./3.) - 1.)**3*p[2] + \
            + ( (1./(x+1.)**(2./3.) - 1.))**2*(6.0-4.0/(x+1.)**(2./3.)) )

        return val

class BirchMurnaghanFit(Fit):
    def __init__(self):
        Fit.__init__(self)

    def func(self, p, x):
        """
        x is assumed to be a percent volume,
        p[0] - E0
        p[1] - V0*B0
        p[2] - B0'
        """
        val = p[0] + 9.0/16.0*p[1]*( (1./(x+1.0)**(2./3.) - 1.)**3*p[2] + \
            + ( (1./(x+1.)**(2./3.) - 1.))**2*(6.0-4.0/(x+1.)**(2./3.)) )

        return val


    def fit(self, xdata, ydata):
        errFunc = lambda p, x, y: (y - self.func(p, x))
      # initial parameters:
        pinit = numpy.zeros(3)
        pinit[1] = (ydata[1]-ydata[0])/(xdata[1] - xdata[0])
        v, success = scipy.optimize.leastsq(errFunc, pinit, args=(xdata, ydata))
        if success < 1 or success > 4:
            print success
            print ydata
            print v
    #   assert success != 1, "fitFreq: Fitting was not successful"
        return v

#class VoluFit(PolynomialFit):
#    def __init__(self, *fitDirective):
#        Fit.__init__(self)
#        if fitDirective[0] != 'polynom':
#            raise Exception('This fit is not implemented')
#        self.type = fitDirective[0]
#        polyOrder = fitDirective[1]
#        hasConstantTerm = fitDirective[2]
#        if hasConstantTerm == True:
#            self.constantTerm = 1
#        else:
#            self.constantTerm = 0
#
#        if polyOrder < 1 or polyOrder > 4:
#            raise Exception('this polyOrder is not supported')
#        self._polyOrder = polyOrder
#        self.order = polyOrder + 1


class FreqFit():
    def __init__(self, prcntVolume, freqs, fitter):
        """indexRange is for files' indexing. prcntVolume - coresponding volume
        expansions"""
        self.fitter = fitter
#        Fit.__init__(self, *fitDirective)
        self._prcntVolume = prcntVolume
        self._freqs = freqs
        self._fit()
        

    def _fit(self):
        import copy
        self._prcntFreqs = numpy.zeros(self._freqs.shape)
        # introduce Omega0 to get rid of 0s in denominator
        Omega0 = copy.copy(self._freqs[0])
        self._freqs0 = Omega0
        for i in range(Omega0.shape[0]):
            for j in range(Omega0.shape[1]):
                if Omega0[i,j] == 0.:
                    Omega0[i, j] = 1
        for i in range(self._freqs.shape[0]):
            self._prcntFreqs[i] = (self._freqs[i] - self._freqs[0])/Omega0
        self._coeff = numpy.zeros( shape = ( Omega0.shape[0], Omega0.shape[1],\
                                              self.fitter.order) )
#        nonlinearity = 0
        for i in range( Omega0.shape[0] ):
            for j in range( Omega0.shape[1] ):
                self._coeff[i,j] = self.fitter.fit(self._prcntVolume, \
                                                 self._prcntFreqs[:,i,j])
#                if self.fitter.__polyOrder > 1:
#                    nonlinearity = nonlinearity + abs(self._coeff[i,j][1])


    def quasiFreqs(self, prcntVol):
        """Returns array of freqs, corresponding to the provided volume expansion"""
        if self.fitter.type != 'polynom':
            raise Exception('Does not have any meaning for the given type of fit')
        nPoints = self._coeff.shape[0]
        nModes = self._coeff.shape[1]
        fittedFreqs = numpy.zeros(shape = (nPoints, nModes) )
        for i in range( nPoints ):
            for j in range( nModes ):
                p = self._coeff[i,j]
                fittedFreqs[i,j] = p[-1]*prcntVol
        return (fittedFreqs+1.0)*self._freqs0


    def fittedFreqs(self, prcntVol):
        """Returns array of freqs, corresponding to the provided volume expansion"""
        nPoints = self._coeff.shape[0]
        nModes = self._coeff.shape[1]
        fittedFreqs = numpy.zeros(shape = (nPoints, nModes) )
        for i in range( nPoints ):
            for j in range( nModes ):
                p = self._coeff[i,j]
                fittedFreqs[i,j] = self.fitter.func(p, prcntVol)
        return (fittedFreqs+1.0)*self._freqs0


    def freqs(self):
        """Returns volume expansion"""
        return self._freqs


    def coeff(self):
        return self._coeff



class ValueFit():
    def __init__(self, prcntVolume, values, fitter):
        """indexRange is for files' indexing. prcntVolume - coresponding volume
        expansions"""
        self.fitter = fitter
#        Fit.__init__(self, *fitDirective)
        self._prcntVolume = numpy.array(prcntVolume)
        self._values = numpy.array(values)
        self._fit()


    def _fit(self):
        if self._values[0] == 0:
            raise Exception('Denominator is 0!')
        prcntValues = (self._values - self._values[0])/self._values[0]
        self._coeff = self.fitter.fit(self._prcntVolume, \
                                            prcntValues)


    def fittedValue(self,prcntVol):
        return (self.fitter.func(self._coeff, prcntVol) + 1.0)*self._values[0]

    def chi2(self):
        fitted = []
        for v in self._prcntVolume:
            fitted.append(self.fittedValue(v))
        fitted = numpy.array(fitted)
        diff = fitted - self._values
        return numpy.sqrt(numpy.dot(diff, diff))

    def values(self ):
        return(self._values)


    def coeff(self):
        return self._coeff


class Fit1D():
    def __init__(self, xValues, yValues, fitter):
        """indexRange is for files' indexing. prcntVolume - coresponding volume
        expansions"""
        self.fitter = fitter
#        Fit.__init__(self, *fitDirective)
        self._xVal = xValues
        self._yVal = yValues
        #self._prcntVolume = numpy.array(prcntVolume)
        #self._values = numpy.array(values)
        self._fit()


    def _fit(self):        
        self._coeff = self.fitter.fit(self._xVal, self._yVal)


    def fittedValue(self, xVal):
        return self.fitter.func(self._coeff, xVal)

    def chi2(self):
        fitted = []
        for v in self._xVal:
            fitted.append(self.fittedValue(v))
        fitted = numpy.array(fitted)
        diff = fitted - self._yVal
        return numpy.sqrt(numpy.dot(diff, diff))

    def values(self ):
        return(self._yVal)


    def coeff(self):
        return self._coeff

if __name__ == "__main__": pass

__author__="markovsk"
__date__ ="$Sep 16, 2009 12:42:48 PM$"