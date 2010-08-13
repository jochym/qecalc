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
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'flfrc': None,
        'flvec': 'matdyn.modes',
        'flfrq': 'matdyn.freq',
        'fldos': 'matdyn.dos'
        }

        self.readSetting(filename, configString, sectionName)
        
#        self._cmdStr = "matdyn.x -inp " + self.setting.matdynInput + " > " + \
#                        self.setting.matdynOutput

                        
    def cmdLine(self):
        return self._getCmdLine('matdyn.x', 'matdynInput', 'matdynOutput')


    def name(self):
        return 'matdyn.x'


    def readSetting(self, filename = None, configString = None, sectionName = None):
        """
        Initializes Setting, QEInput and QEOutout classes 
        and synchronizes with QEInput object
        """
        
        configDic = {
        'matdynInput': 'matdyn.in',
#        'flfrc': None,
        'matdynOutput': 'matdyn.out',
#        'flvec': None,
#        'flfrq': None,
#        'fldos': None
        }
        
        QETask.readSetting(self, filename, configString)
        
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)

        self.input = MatdynInput(filename = self.setting.get('matdynInput'))
        # add pointer to setting for input filenames synchronization 
        self.input._setting = self.setting        
        self.output = QEOutput(self.setting, type = 'matdyn')
        
        if filename != None or configString != None:
            self.syncSetting()        


    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.parse()
        
        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'input', varName, \
                                                self.input, self._path_defaults)

        #self._syncPathInNamelist('flfrc', 'input', 'matdynflfrc')
        #self._syncPathInNamelist('flfrq', 'input', 'matdynflfrq')
        #self._syncPathInNamelist('flvec', 'input', 'matdynflvec')
        #self._syncPathInNamelist('fldos', 'input', 'matdynfldos')


                    
if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 5:19:13 PM$"
