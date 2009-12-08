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

from baseoutput import BaseOutput

class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        self.parsers = {
                        'qpoints'       : self.getQpoints,
                        }

    def getQpoints(self, setting):
        """
        Extract ireducible qpoints from  ph output
        """
        #read Espresso output into memory:
        file = open(setting.phOutput)
        phOut = file.readlines()
        posList =  \
        [i for i,line in enumerate(phOut) if 'Dynamical matrices for (' in line]
        i = posList[-1] + 3
        qpoints = []
        while len(phOut[i].split()) == 4:
            qpoints.append([ float(w) for w in phOut[i].split()[1:]])
            i = i + 1
        return [(qpoints, None)]


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Nov 3, 2009 6:53:12 PM$"
