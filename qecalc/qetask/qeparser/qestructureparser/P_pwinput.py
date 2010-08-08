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
    
    def __init__(self, qeInput):
        QEStructureParser.__init__(self, qeInput)
        self.format = "pwinput"
    
    
    def parseStr(self, s):
        """Create Structure instance from a string."""
        input = PWInput()
        input.readString(s)
        return self.__genStructure(input)
         

    def parse(self, filename):
        """Create Structure instance from a file."""
        self.filename = filename 
        input = PWInput(filename = filename)
        return self.__genStructure(input)
            
            
    def __genStructure( self, input ):
        #****** new version *******
        input.parse()
        #input.update( qeInput = self._qeInput, forceUpdate = True  )
        stru = QEStructure( input.structure, qeInput = self._qeInput )        
        #print "input:", stru._qeInput.toString()
                
        return stru        
#     ******** Old version        
        input.parse()
        stru = QEStructure(qeInput = self._qeInput)        
        
        input.structure.lattice._qeInput = self._qeInput
        input.structure._qeInput = self._qeInput
        #input.structure._qeInput.update( forceUpdate = True )
        #print 'Input parser', stru._qeInput.toString()         
        # create a shallow copy of all source attributes
        stru.__dict__.update(input.structure.__dict__)   
        stru[:] = input.structure
        #stru._qeInput.update( forceUpdate = True )
        return stru
    #*********************        
                    
    
def getParser(qeInput):
    return P_pwinput(qeInput)
 