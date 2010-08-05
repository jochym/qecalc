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
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setParallel()
        
        configDic = {
        'phInput': 'ph.in',
        'phOutput': 'ph.out',
        'fildyn'  : None,
#        'fildrho' : None,
#        'fildvscf' : None
        }

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'fildyn': 'matdyn',
        'fildrho': '',
        'fildvscf': '',
        'outdir': './',
        'prefix': 'pwscf'
        }
        
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)

        self.input = PHInput(filename = self.setting.get('phInput'))
        # add pointer to setting for input filenames synchronization 
        self.input._setting = self.setting        
        self.output = QEOutput(self.setting, type='ph')


    def cmdLine(self):
        return self._getCmdLine('ph.x', 'phInput', 'phOutput')
    

    def name(self):
        return 'ph.x'

    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.parse()

        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'inputph', varName, \
                                                self.input, self._path_defaults)

        #self._syncPathInNamelist('fildyn', 'inputph', 'phfildyn')
        #self._syncPathInNamelist('fildrho', 'inputph', 'phfildrho')
        #self._syncPathInNamelist('fildvscf', 'inputph', 'phfildvscf')
        #self._syncPathInNamelist('outdir', 'inputph', 'outDir')