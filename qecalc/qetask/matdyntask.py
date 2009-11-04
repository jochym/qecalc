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
from qeparser.matdyninput import MatdynInput
from qeparser.qeoutput import QEOutput

class MatdynTask(QETask):
    def __init__(self, filename, cleanOutDir = None):
        QETask.__init__(self, filename, cleanOutDir)

        self.name = 'matdyn.x'
        
        configDic = {
        'matdynInput': 'matdyn.in',
        'matdynOutput': 'matdyn.out',
        'matdynModes': 'matdyn.modes',
        'matdynFreqs': 'matdyn.freq',
        'matdynfldos': 'matdyn.phdos'
        }
        self.setting.section(self.name, configDic)

        self.input = MatdynInput(filename = self.setting.matdynInput)
        self.output = QEOutput(self.setting, type = 'matdyn')
        self._cmdStr = "matdyn.x -inp " + self.setting.matdynInput + " > " + \
                        self.setting.matdynOutput

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 5:19:13 PM$"
