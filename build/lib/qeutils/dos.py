#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To change this template, choose Tools | Templates
# and open the template in the editor.


class PhononDos:
    def __init__(self):
        self._freqs = None
        self._modes = None
        self._qpts = None


#    def loadPhonons(self, fname = None):
#        self._modes, self._freqs, self._qpts =  \
#                                    self.matdyn.output.property('multi phonon')

    def getPhonons(self):
        return self._modes, self._freqs, self._qpts

    def qpoints(self):
        return self._qpts

    def freqs(self):
        return self._freqs

    def modes(self):
        return self._modes


    def setRange(self, minOmega, maxOmega, deltaOmega):
        if minOmega == None:
            minOmega = numpy.min(self._freqs)
        if maxOmega == None:
            maxOmega = numpy.max(self._freqs)
        if minOmega > maxOmega: minOmega, maxOmega = maxOmega, minOmega
        if deltaOmega == None:
            deltaOmega = (maxOmega - minOmega)/200.0
        return minOmega, maxOmega, deltaOmega

    def set(self, modes, freqs, qpts):
        self._modes, self._freqs, self._qpts =   modes, freqs, qpts


    def DOS(self, minOmega = None, maxOmega = None, deltaOmega = None):
        minOmega, maxOmega, deltaOmega =  \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histOmega = numpy.zeros(nPoints)
        norm = 0.0
        for cell_freqs in self._freqs:
            for omega in cell_freqs:
                idx = int( (abs(omega) - minOmega)/deltaOmega )
                if idx < len(histOmega):
                    histOmega[idx] = histOmega[idx] + 1.0
                    norm = norm + 1.0

        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        return histOmega/norm, axis

    def partDOS(self, atomSymbol, minOmega = None, maxOmega = None, deltaOmega = None):
        from numpy import real
        minOmega, maxOmega, deltaOmega   =      \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histPartOmega = numpy.zeros(nPoints)
        norm = 0.0
        for iAtom, atom in enumerate(self.structure.diffpy()):
            if atomSymbol == atom.element:
                for cell_freqs, vectors in zip(self._freqs, self._modes):
                    for omega, vector in zip(cell_freqs, vectors[:,iAtom,:]):
                        idx = int( (abs(omega) - minOmega)/deltaOmega )
                        if idx < len(histPartOmega):
                            weight = (real(vector*vector.conjugate())).sum()
                            histPartOmega[idx] = histPartOmega[idx] + weight
                            norm = norm + weight
        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        return histPartOmega/norm, axis

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Nov 16, 2009 1:27:22 PM$"
