# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="markovsk"
__date__ ="$Sep 17, 2009 10:16:24 AM$"

import numpy

class Thermodyn():
    def __init__(self, axis, dos, units):
        """dos is 2d array; x:y format"""
        #units
        self.h = 6.62606896e-34
        self.cmToHz = 29979245800.0
        self.kb = 1.3806504e-23
        self.avogadro = 6.022214179e23
        self.Ry = 4.587420897e+17
        self.energy_units = {
                              'cm-1' : 29979245800.0,
                              'meV'  : 241798840766.20228
        }
        self.setDOS(axis, dos, units)


    def setDOS(self, axs, dos, units):
        """dos is assumed to be normalized by phonon degrees of freedom in
        the unit cell (e.g. its norm is 9 for MgB2)"""
        if dos != None and units != None and axs != None:
            if units not in self.energy_units:
                raise Exception('Energy units are not known')
        else:
            raise Exception('dos, axis,  and units should be provided!')
        self.deltaE = axs[1] - axs[0]
        self.units = units
        self.axis = axs[1:]*self.energy_units[self.units]
        self.g = dos[1:]
#            else:
#                raise Exception('dos = None!')




class PhononThermodynamics(Thermodyn):
    def __init__(self, axis, dos):
        """dos is assumed to be normalized by phonon degrees of freedom in
        the unit cell (e.g. its norm is 9 for MgB2)"""
        Thermodyn.__init__(self, axis, dos, 'cm-1')

    def zeroPointEnergy(self, axis = None, dos = None):
        self._checkSetDOS(axis, dos)
        zpEnergy = self.h*self.axis*0.5*self.g
        return zpEnergy.sum()*self.Ry*self.deltaE

    def _checkSetDOS(self, axis, dos):
        if dos == None:
            if axis != None:
                raise Exception('Should set axis together with dos')
            # use self._g and self._axis
        else:
            self.setDOS(axis, dos, self.units)

    def freeEnergy(self, T, axis = None, dos = None):
        self._checkSetDOS(axis, dos)
        if T == 0.0:
            return self.zeroPointEnergy()
        else:
            arg = (self.h/self.kb)*self.axis/T
#        print 1.0 - 1.0/numpy.exp(1.0/arg)
            F = self.kb*T*( 0.5*arg + numpy.log(1.0 - 1.0/numpy.exp(arg)) )*self.g
            return F.sum()*self.Ry*self.deltaE
        
    def Cv(self, T, axis = None, dos = None):
        self._checkSetDOS(axis, dos)
        arg = (self.h/self.kb)*self.axis/T
        expArg = numpy.exp(arg)
        Cv = self.g*arg**2*(expArg/(expArg-1.0)**2)
        return Cv.sum()*self.kb*6.022214179e23*self.deltaE
        
if __name__ == "__main__":
    print "Hello World";
