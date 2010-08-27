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


from qeinput import QEInput

from qecalc.qetask.setting import Setting

class QESInput(QEInput):
    """
    Layer between QEInput and separate input classes to bring additional 
    functionality to QEInput
    QESInput also makes sure QEInput is synchronized with Setting class 
    """
    def __init__(self, filename=None, config=None, type = 'pw', setting = None,\
                                                                  parse = True):       
                        
        self._setting = setting
        if isinstance(filename, Setting):
            self._setting = filename
        if self._setting != None:
            fname = self._setting.get( type + 'Input')
        else:
            fname = filename
            
        QEInput.__init__(self,fname, config, type = type, parse = parse)
        
        if parse and (filename or config  or setting ):
            self.parse()        
    
    def parse(self):
        """
        Parses the configuration file and stores the values in qe dictionary
        """
        if self._setting != None:
            self.readFile( filename = self._setting.get(self.type() + 'Input'))
        else:
            QEInput.parse(self)  
    
    
    def save(self, filename=None):
        if filename == None:
            if self._setting != None:
                self.filename = self._setting.get(self.type()+'Input')
            QEInput.save(self)
        else:
            QEInput.save(self, filename)
            
            
    def readFile(self, filename):
        """
        Reads and parses configuration input from file
            filename: (str) -- File name
        """
        QEInput.readFile(self, filename)
        if self._setting != None:
            self.filename = self._setting.get(self.type() + 'Input')
               