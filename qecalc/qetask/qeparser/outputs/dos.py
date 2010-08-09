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


import numpy
from baseoutput import BaseOutput

class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        # dictionary with list of alternative property names,
        # not case and white space sensitive:
        self._propertyNamesDic = { 'electron dos'      : ['dos', 'electrondos', 'electronicdos'],
                              }        
        self.parsers = {
                'electron dos'   : self.getDOS,
                }

    def getDOS(self,setting):
        """
        Extracts electronic DOS from QE output. Return axis and values
        """
        dos = numpy.loadtxt(setting.get('fldos'))
        return [(dos[0:,0], 'eV'), (dos[0:,1], None)]



if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Nov 10, 2009 12:16:38 PM$"
