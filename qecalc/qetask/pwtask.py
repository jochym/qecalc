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

        self.setParallel()

        # QE input file's path containing variables' defaults (will be moved to
        # QE input parser)
        self._path_defaults = {
        'outdir': './',
        'pseudo_dir': './',
        'prefix': 'pwscf'
        }        
        
        
        self.readSetting(filename, configString, sectionName)


    def cmdLine(self):
        return self._getCmdLine('pw.x', 'pwInput', 'pwOutput')


    def name(self):
        return 'pw.x'


    def readSetting(self, filename = None, configString = None, sectionName = None):
        """
        Initializes Setting, QEInput and QEOutout classes 
        and synchronizes with QEInput object
        """
        
        # pw main input and output
        configDic = {
        'pwInput': 'scf.in',
        'pwOutput': 'scf.out',
        }        
        
        QETask.readSetting(self, filename, configString)
        
        #self._setSetiingInputOutput(configDic, input, output,  filename, configString, sectionName)
        
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)
        self.input = PWInput( setting = self.setting ) #self.setting.get('pwInput') )
        # add pointer to setting for input filenames synchronization 
        #self.input._setting = self.setting
        self.output = QEOutput(self.setting, type='pw')
        
        if filename != None or configString != None:
            self.syncSetting()
            

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