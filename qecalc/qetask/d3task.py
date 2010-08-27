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
from qeparser.d3input import D3Input
from qeparser.qeoutput import QEOutput

class D3Task(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setSerial()

        # ****************** Task Specifics ************************************
        self._inputConstructor = 'QESInput'
        # input/output defaults        
        self._configDic = {
        'd3Input': 'd3.in',                           
        'd3Output': 'd3.out',
        }                           
        self._path_defaults = {
        'fildyn': 'd3dyn',
        'fildrho': ' ',
        'fild0rho': ' ',
        'outdir': './',
        'prefix': 'pwscf'
        }
        self._type = 'd3'
        # **********************************************************************        
        
        self.readSetting(filename, configString, sectionName)


    def cmdLine(self):
        return self._getCmdLine('d3.x', 'd3Input', 'd3Output')


    def name(self):
        return 'd3.x'

    
    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        self.input.parse()
        for varName in self._path_defaults.keys():
            self.setting.syncPathInNamelist(varName, 'inputph', varName, \
                                                self.input, self._path_defaults)

        #self._syncPathInNamelist('fildyn', 'inputph', 'd3fildyn')
        #self._syncPathInNamelist('fildrho', 'inputph', 'd3fildrho')
        #self._syncPathInNamelist('fild0rho', 'inputph', 'd3fild0rho')
        #self._syncPathInNamelist('outdir', 'inputph', 'outDir')
