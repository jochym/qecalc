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
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        configDic = {
        'dosInput' : 'dos.in',
        'dosOutput': 'dos.out',
        #'dosfldos' : None
        }

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'fldos': 'fldos.dos',
        'outdir': './',
        'prefix': 'pwscf',
        }
        
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)
        self.input = QEInput(self.setting.get('dosInput'), type = 'dos')
        # add pointer to setting for input filenames synchronization 
        self.input._setting = self.setting        
        self.output = QEOutput(self.setting, type='dos')
#        self._cmdStr = self.setting.paraPrefix + " dos.x " +  \
#                       " -inp " + \
#                       self.setting.dosInput + " > " + \
#                       self.setting.dosOutput + "< /dev/null"


    def cmdLine(self):
        return self._getCmdLine('dos.x', 'dosInput', 'dosOutput')


    def name(self):
        return 'dos.x'


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
__date__ ="$Nov 10, 2009 11:37:18 AM$"
