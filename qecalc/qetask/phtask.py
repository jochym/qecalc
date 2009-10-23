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

class PHTask(QETask):
    def __init__(self, setting, cleanOutDir = False):
        QETask.__init__(self, setting, cleanOutDir)
        self.input = QEInput(self.setting.phInput, type = 'ph')
        self.output = QEOutput(self.setting, type='pw')
        self.cmdStr = self.setting.paraPrefix + " ph.x " +  \
                      self.setting.paraPostfix + " -inp " + \
                      self.setting.phInput + " > " + \
                      self.setting.phOutput + "< /dev/null"                 
        self.name = 'ph.x'