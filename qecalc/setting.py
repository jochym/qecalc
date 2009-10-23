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
class Setting:
    def __init__(self, fname= None):
        import os
        import ConfigParser
        
       # Default values, see explanations below:        
        configDic = {
        'isMetallic': 'True',
        'useTorque' : 'False',
        'torqueResourceList': '-l nodes=1:ppn=1',
        'numNodes'  : '1',
        'numProc': '1',
        'paraPrefix': '',
        'paraPostfix': '',
        'pwscfInput': 'scf.in',
        'pwscfOutput': 'scf.out',
        'phInput': 'ph.in',
        'phOutput': 'ph.out',
        'dynmatInput': 'dynmat.in',
        'dynmatOutput': 'dynmat.out',
        'q2rInput': 'q2r.in',
        'q2rOutput': 'q2r.out',
        'matdynInput': 'matdyn.in',
        'matdynOutput': 'matdyn.out',        
        'matdynModes': 'matdyn.modes',
        'matdynFreqs': 'matdyn.freq',
        'matdynfldos' : 'matdyn.phdos'
        }

        try:
          if fname == None:
             raise NameError("Config should be initialized with a filename")
        except NameError:
            raise
        self.configFileName = fname
        self.config = ConfigParser.SafeConfigParser(configDic)
        self.config.read(fname)

        #workDir = os.getcwd() + '/'
        
        # all the relevant input files must be preconfiguered for specific tasks 
        # before using this class

        # parallelization parameters

        self.useTorque =  self.config.getboolean('Setting','useTorque')

        self.numProc = self.config.getint('Setting','numProc')
        self.paraPrefix = self.config.get('Setting', 'paraPrefix')
        self.paraPostfix = self.config.get('Setting', 'paraPostfix')

        self.pwscfInput = self.config.get('Setting', 'pwscfInput')
        # pwscf output file relevant to 'total energy' as well as 'geometry' tasks           
        self.pwscfOutput = self.config.get('Setting', 'pwscfOutput')
        
        self.phInput =  self.config.get('Setting', 'phInput')
        self.phOutput =  self.config.get('Setting', 'phOutput')
        
        # dynmat input/output file relevant to 'single phonon' task        
        self.dynmatInput =  self.config.get('Setting', 'dynmatInput')
        self.dynmatOutput =  self.config.get('Setting', 'dynmatOutput')
        
        # input/output files relevant to 'multiple phonon' task    
        self.q2rInput =  self.config.get('Setting', 'q2rInput')
        self.q2rOutput =  self.config.get('Setting', 'q2rOutput')
        self.matdynInput =  self.config.get('Setting', 'matdynInput')
        self.matdynOutput =  self.config.get('Setting', 'matdynOutput')
        self.matdynModes =  self.config.get('Setting', 'matdynModes')
        self.matdynFreqs =  self.config.get('Setting', 'matdynFreqs')
        self.matdynfldos =  self.config.get('Setting', 'matdynfldos')
        