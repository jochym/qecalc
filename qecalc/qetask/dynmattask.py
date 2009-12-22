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
from qeparser.qeinput import QEInput
from qeparser.qeoutput import QEOutput

class DynmatTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        configDic = {
        'dynmatInput': 'dyn.in',
        'dynmatOutput': 'dyn.out',
 #       'dynmatfildyn': None,
 #       'dynmatfilout': None,
 #       'dynmatfilmol': None,
 #       'dynmatfilxsf': None,
        }

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'fildyn': 'matdyn',
        'filout': 'dynmat.out',
        'filmol': 'dynmat.mold',
        'filxsf': 'dynmat.axsf',
        }
        
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)
        
        self.input = QEInput(filename = self.setting.get('dynmatInput'), type = 'dynmat')
        self.output = QEOutput(self.setting, type = 'dynmat')
        #self._cmdStr = "dynmat.x < " + self.setting.dynmatInput
        #+ " > " + \
        #               self.setting.dynmatOutput

        
    def cmdLine(self):
        if self.setting.get('dynmatOutput') == self.setting.get('filout'):
            raise NameError('Naming conflict: Change dynmatOutput value')
        return self._getCmdLine('dynmat.x', 'dynmatInput', 'dynmatOutput')


    def name(self):
        return 'dynmat.x'


    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """

        self.input.parse()

        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'input', varName, \
                                                self.input, self._path_defaults)

#        self._syncPathInNamelist('fildyn', 'input', 'dynmatfildyn')
#        self._syncPathInNamelist('filout', 'input', 'dynmatfilout')
#        self._syncPathInNamelist('filmol', 'input', 'dynmatfilmol')
#        self._syncPathInNamelist('filxsf', 'input', 'dynmatfilxsf')


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 3:24:51 PM$"
