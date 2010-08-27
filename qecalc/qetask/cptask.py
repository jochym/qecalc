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


from qetask import QETask
from qeparser.cpinput import CPInput
from qeparser.qeoutput import QEOutput

class CPTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = False,\
                                                            sectionName = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.setParallel()

        # ****************** Task Specifics ************************************
        self._inputConstructor = 'CPInput'
        # input/output defaults        
        self._configDic = {
        'cpInput': 'cp.in',
        'cpOutput': 'cp.out',
        }        
        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'outdir': './',
        'pseudo_dir': './',
        'prefix': 'cp'
        }
        self._type = 'cp'
        # **********************************************************************        
        
        self.readSetting(filename, configString, sectionName)


    def cmdLine(self):
        return self._getCmdLine('cp.x', 'cpInput', 'cpOutput')


    def name(self):
        return 'cp.x'


    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """

        self.input.parse()

        self.setting.syncPathInNamelist('outdir', 'control', 'outdir', \
                                                self.input, self._path_defaults)
        self.setting.syncPathInNamelist('pseudo_dir', 'control', 'pseudo_dir', \
                                                self.input, self._path_defaults)
        self.setting.syncPathInNamelist('prefix', 'control', 'prefix', \
                                                self.input, self._path_defaults)


        # Solve for pseudopotential locations
        species = self.input.structure.atomicSpecies
        for specie in species.keys():
            setattr(self.setting._paths, 'ps' + specie, self.setting.get('pseudo_dir') + species[specie].potential )

