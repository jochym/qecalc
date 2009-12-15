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
# To change this template, choose Tools | Templates
# and open the template in the editor.

from qetask import QETask
from qeparser.qeinput import QEInput
from qeparser.qeoutput import QEOutput

class DOSTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        #self.name = 'dos.x'

        # pwscf input and output
        configDic = {
        'dosInput' : 'dos.in',
        'dosOutput': 'dos.out',
        'fldos'      : 'fldos.dos'
        }
        self.setting.section(self.name(), configDic)
        self.input = QEInput(self.setting.dosInput, type = 'dos')
        self.output = QEOutput(self.setting, type='dos')
#        self._cmdStr = self.setting.paraPrefix + " dos.x " +  \
#                       " -inp " + \
#                       self.setting.dosInput + " > " + \
#                       self.setting.dosOutput + "< /dev/null"


    def cmdLine(self):
        return self.setting.paraPrefix + " dos.x " +  \
                       " -inp " + \
                       self.setting.dosInput + " > " + \
                       self.setting.dosOutput + "< /dev/null"


    def name(self):
        return 'dos.x'


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Nov 10, 2009 11:37:18 AM$"
