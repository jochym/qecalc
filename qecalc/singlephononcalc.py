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

from qecalc import QECalc
import numpy

from qetask.pwtask import PWTask
from qetask.phtask import PHTask
from qetask.dynmattask import DynmatTask
from qetask.pwphmerger import PWPHMerger

class SinglePhononCalc(QECalc):
    """ Calc for multi phonon calculations:
        Task list:
          pw     -- PWTask
          ph     -- PHTask
          dynmat -- DynmatTask
          pwph   -- PWPHMerger - task merger used for submission of pw.x and ph.x
                    commands in a single command string
          taskList = [pwph, dynmat]
          Example:
            >>> phonCalc = SinglePhononCalc('config.ini')
            >>> phonCalc.launch()
            >>> print phonCalc.dynmat.output('single phonon')
    """
    def __init__(self, filename):
        QECalc.__init__(self, filename)
        self.pw = PWTask(self.setting)
        self.ph = PHTask(self.setting)
        self.dynmat = DynmatTask(self.setting)
        self.pwph = PWPHMerger(self.pw,self.ph, cleanOutDir = True)
        self.taskList = [self.pwph, self.dynmat]

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 2:52:30 PM$"
