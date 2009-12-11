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
from qeparser.phinput import PHInput
from qeparser.qeoutput import QEOutput

class PHTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        #self.name = 'ph.x'
        
        configDic = {
        'phInput': 'ph.in',
        'phOutput': 'ph.out',
        'phFildyn'  : 'matdyn'
        }
        self.setting.section(self.name(), configDic)

        self.input = PHInput(filename = self.setting.phInput)
        self.output = QEOutput(self.setting, type='ph')


    def cmdLine(self):
        return  self.setting.paraPrefix + " ph.x " +  \
                       self.setting.paraPostfix + " -inp " + \
                       self.setting.phInput + " > " + \
                       self.setting.phOutput + "< /dev/null"

    def _syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.namelist('inputph').add('fildyn', \
                                              "'" + self.setting.phFildyn + "'")
        self.input.namelist('inputph').add('outdir', \
                                               "'" +  self.setting.outDir + "'")

    def name(self):
        return 'ph.x'