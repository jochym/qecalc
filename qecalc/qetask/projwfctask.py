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
from qeparser.qesinput import QESInput
from qeparser.qeoutput import QEOutput

class ProjwfcTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        # ****************** Task Specifics ************************************
        self._inputConstructor = 'QESInput'
        # input/output defaults        
        self._configDic = {
        'projwfcInput' : 'pdos.in',
        'projwfcOutput': 'pdos.out',
        }                           
        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'outdir': './',
        'prefix': 'pwscf',
        }
        self._type = 'projwfc'
        # **********************************************************************        
        
        self.readSetting(filename, configString, sectionName)




    def cmdLine(self):
        return self._getCmdLine('projwfc.x', 'projwfcInput', 'projwfcOutput')


    def name(self):
        return 'projwfc.x'


    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.parse()
        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'inputpp', varName, \
                                                self.input, self._path_defaults)

if __name__ == "__main__": pass

__author__="Nikolay Markovskiy"
__date__ ="$Jan 20, 2010 5:24:23 PM$"
