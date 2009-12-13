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
from phdos import PhononDOS

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
        self.dispersion = PHDispersion(self.pw.input.structure.lattice, self.matdyn)
        self.dos = PhononDOS(self.matdyn)
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


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 1:37:29 PM$"
