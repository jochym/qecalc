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

from phdispersion import PHDispersion

from qecalc import QECalc

from qetask.pwtask import PWTask
from qetask.phtask import PHTask
from qetask.q2rtask import Q2RTask
from qetask.matdyntask import MatdynTask
from qetask.pwphmerger import PWPHMerger

class MultiPhononCalc(QECalc):
    """ Calc for multi phonon calculations:
    
    Task list:
      pw     -- PWTask

      ph     -- PHTask

      q2r    -- Q2RTask

      matdyn -- MatdynTask

      pwph   -- PWPHMerger - task merger used for submission of pw.x and ph.x
      commands in a single command string

      taskList = [pwph, q2r, matdyn]
      
    Example:
      
      >>> mphonCalc = MultiPhononCalc('config.ini')
      >>> mphon.launch()
      >>> print mphon.pw.output.property('total energy')
      >>> print mphon.matdyn.output.listParsers()
      >>> print mphon.matdyn.output.property('phonon dos')
      >>> polVecs, freqs, qpoints =  mphon.lookupProperty('multi phonon')
            
    """
    def __init__(self, filename):
        QECalc.__init__(self)
        self._freqs = None
        self._modes = None
        self._qpts = None
        self.pw = PWTask(filename)
        self.ph = PHTask(filename)
        self.q2r = Q2RTask(filename)
        self.matdyn = MatdynTask(filename)
        self.pwph = PWPHMerger(self.pw,self.ph, cleanOutDir = self.pw.input.outDir())
        self.dispersion = PHDispersion(self.pw.input.structure, self.matdyn)
        self.taskList = [self.pwph, self.q2r, self.matdyn]


    def syncInputs(self):
        for task in self.taskList:
            task.input.parse()

        # remove amass from phonon input
        for param in self.ph.input.namelist('input').params:
            if 'amass(' in param:
                self.ph.input.namelist('input').remove(param)
        # initialize amass based on PW input
        for i, atom in enumerate(self.pw.input.structure.atomicSpecies):
            amass = 'amass(' + str(i+1) + ')'
            self.ph.input.namelist('input').add(amass,atom.mass)
        #syncronise outdir based on PW input
        self.ph.input.namelist('input').add('outdir', \
                             self.pw.input.namelist('control').param('outdir'))


    def loadPhonons(self, fname = None):
        self._modes, self._freqs, self._qpts =  \
                                    self.matdyn.output.property('multi phonon')

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
__date__ ="$Oct 19, 2009 1:37:29 PM$"
