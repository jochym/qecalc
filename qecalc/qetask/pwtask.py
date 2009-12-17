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
    def __init__(self, filename = None,configString = None, cleanOutDir = False):
        QETask.__init__(self, filename, configString, cleanOutDir)

        # pwscf main input and output
        configDic = {
        'pwscfInput': 'scf.in',
        'pwscfOutput': 'scf.out',
        }

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'outdir': ''
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

    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """

        self.input.parse()
        
        self.setting.syncPathInNamelist('outdir', 'control', 'outdir', \
                                                self.input, self._path_defaults)
        #print 'Preved'
        #print self.setting.outDir
