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
from qeparser.qesinput import QESInput
from qeparser.qeoutput import QEOutput

class Q2RTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        # ****************** Task Specifics ************************************
        self._inputConstructor = 'QESInput'
        # input/output defaults
        self._configDic = {
        'q2rInput': 'q2r.in',
        'q2rOutput': 'q2r.out',
        }        
        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'fildyn': None,
        'flfrc': None,
        }        
        self._type = 'q2r'
        # **********************************************************************

        self.readSetting(filename, configString, sectionName)

    def cmdLine(self):
        return self._getCmdLine('q2r.x', 'q2rInput', 'q2rOutput')


    def name(self):
        return 'q2r.x'


    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        
        self.input.parse()

        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'input', varName, \
                                                self.input, self._path_defaults)
        #self._syncPathInNamelist('fildyn', 'input', 'q2rfildyn')
        #self._syncPathInNamelist('flfrc', 'input', 'q2rflfrc')



if __name__ == "__main__": pass

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 3:12:40 PM$"
