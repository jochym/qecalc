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


from qecalc.qetask.qeparser.pwinput import PWInput
from qecalc.qetask.qeparser.qestructureparser import *

class P_pwinput(QEStructureParser):
    """
        Parser for pw.x input format
    """
    
    def __init__(self, qeConf):
        QEStructureParser.__init__(self, qeConf)
        self.format = "pwinput"
    
    
    def parseStr(self, s):
        """Create Structure instance from a string."""
        input = PWInput(config = s)
        return self.__genStructure(input)
         

    def parse(self, filename):
        """Create Structure instance from a file."""
        self.filename = filename 
        input = PWInput(filename = filename)
        return self.__genStructure(input)
            
            
    def __genStructure( self, input ):
        
        input.parse()
        
        stru = QEStructure(qeConf = self.qeConf)
        
        input.structure.lattice.qeConf = self.qeConf
     
        stru.structure = input.structure.structure
        stru.lattice = input.structure.lattice
        stru.atomicSpecies = input.structure.atomicSpecies
        stru.optConstraints = input.structure.optConstraints
        stru.nat = input.structure.nat
        stru.ntyp = input.structure.ntyp
        stru.atomicPositionsType = input.structure.atomicPositionsType
        
        return stru
    
def getParser(qeConf):
    return P_pwinput(qeConf)
 