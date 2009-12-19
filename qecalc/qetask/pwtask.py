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
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        # pw main input and output
        configDic = {
        'pwInput': 'scf.in',
        'pwOutput': 'scf.out',
        }

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'outdir': './'
        }

        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)
        self.input = PWInput( self.setting.get('pwInput') )
        self.output = QEOutput(self.setting, type='pw')
#        self._cmdStr = self.setting.paraPrefix + " pw.x " +  \
#                       self.setting.paraPostfix + " -inp " + \
#                       self.setting.pwscfInput + " > " + \
#                       self.setting.pwscfOutput + "< /dev/null"


    def cmdLine(self):
        return self.setting.get('paraPrefix') + " pw.x " +  \
                       self.setting.get('paraPostfix') + " -inp " + \
                       self.setting.get('pwInput') + " > " + \
                       self.setting.get('pwOutput') + "< /dev/null"


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
