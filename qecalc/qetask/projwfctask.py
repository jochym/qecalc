#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
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

class ProjwfcTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'outdir': './',
        'prefix': 'pwscf',
        }
        
        self.readSetting(filename, configString, sectionName)




    def cmdLine(self):
        return self._getCmdLine('projwfc.x', 'projwfcInput', 'projwfcOutput')


    def name(self):
        return 'projwfc.x'


    def readSetting(self, filename = None, configString = None, sectionName = None, parse = True):
        """
        Initializes Setting, QEInput and QEOutout classes 
        and synchronizes with QEInput object
        """
        
        configDic = {
        'projwfcInput' : 'pdos.in',
        'projwfcOutput': 'pdos.out',
        }     
        
        QETask.readSetting(self, filename, configString)
        
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)
        self.input = QEInput(self.setting.get('projwfcInput'), type = 'projwfc')
        # add pointer to setting for input filenames synchronization 
        self.input._setting = self.setting        
        self.output = QEOutput(self.setting, type='projwfc')
        
        if filename != None or configString != None:
            self.syncSetting()


    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.parse()
        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'inputpp', varName, \
                                                self.input, self._path_defaults)

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Jan 20, 2010 5:24:23 PM$"
