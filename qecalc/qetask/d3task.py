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

class D3Task(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        #self.name = 'ph.x'

        configDic = {
        'd3Input': 'd3.in',
        'd3fildyn': 'd3dyn',
        'd3Output': 'd3.out'
        }
        self.setting.section(self.name(), configDic)

        self.input = QEInput(filename = self.setting.d3Input, type = 'd3')
        self.output = QEOutput(self.setting, type='d3')
#        self._cmdStr = self.setting.paraPrefix + " ph.x " +  \
#                       self.setting.paraPostfix + " -inp " + \
#                       self.setting.phInput + " > " + \
#                       self.setting.phOutput + "< /dev/null"


    def cmdLine(self):
        return  "d3.x < " +  \
                       self.setting.d3Input + " > " + \
                       self.setting.d3Output + "< /dev/null"


    def name(self):
        return 'd3.x'