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

from qecalc import QECalc
import numpy

from qetask.taskpw import PWTask
from qetask.taskph import PHTask
from qetask.taskq2r import Q2RTask
from qetask.taskdynmat import DynmatTask
from qetask.mergerpwph import PWPHMerger

class SinglePhononCalc(QECalc):
    def __init__(self, fname):
        QECalc.__init__(self, fname)
        self.pw = PWTask(self.setting)
        self.ph = PHTask(self.setting)
        self.dynmat = DynmatTask(self.setting)
        self.pwph = PWPHMerger(self.pw,self.ph, cleanOutDir = True)
        self.taskList = [self.pwph, self.dynmat]

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 2:52:30 PM$"
