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

from qetask import QETask
from qeparser.qeinput import QEInput
from qeparser.qeoutput import QEOutput

class MatdynTask(QETask):
    def __init__(self, setting, cleanOutDir = False):
        QETask.__init__(self, setting, cleanOutDir)
        self.input = QEInput(filename = self.setting.matdynInput, type = 'matdyn')
        self.output = QEOutput(self.setting, type = 'matdyn')
        self.cmdStr = "matdyn.x -inp " + self.setting.matdynInput + " > " + \
        self.setting.matdynOutput
        self.name = 'matdyn.x'

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 5:19:13 PM$"
