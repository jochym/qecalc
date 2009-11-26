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

class DynmatTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.name = 'dynmat.x'

        configDic = {
        'dynmatInput': 'dynmat.in',
        'dynmatOutput': 'dynmat.out'
        }
        self.setting.section(self.name, configDic)
        
        self.input = QEInput(filename = self.setting.dynmatInput, type = 'dynmat')
        self.output = QEOutput(self.setting, type = 'dynmat')
        self._cmdStr = "dynmat.x < " + self.setting.dynmatInput
        #+ " > " + \
        #               self.setting.dynmatOutput
    def _syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.namelist('input').remove('filout')
        self.input.namelist('input').add('filout', self.setting.dynmatOutput)


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 3:24:51 PM$"
