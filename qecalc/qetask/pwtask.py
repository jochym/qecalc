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
from qeparser.pwinput import PWInput
from qeparser.qeoutput import QEOutput

class PWTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        #self.name = 'pw.x'

        # pwscf input and output
        configDic = {
        'pwscfInput': 'scf.in',
        'pwscfOutput': 'scf.out',
        }
        self.setting.section(self.name(), configDic)
        self.input = PWInput(self.setting.pwscfInput)
        self.output = QEOutput(self.setting, type='pw')
#        self._cmdStr = self.setting.paraPrefix + " pw.x " +  \
#                       self.setting.paraPostfix + " -inp " + \
#                       self.setting.pwscfInput + " > " + \
#                       self.setting.pwscfOutput + "< /dev/null"


    def cmdLine(self):
        return self.setting.paraPrefix + " pw.x " +  \
                       self.setting.paraPostfix + " -inp " + \
                       self.setting.pwscfInput + " > " + \
                       self.setting.pwscfOutput + "< /dev/null"


    def name(self):
        return 'pw.x'