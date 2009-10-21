#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
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
#from qestructure import QEStructure

class PWInput(QEInput):
    def __init__(self, filename=None, config=None):
        QEInput.__init__(self,filename, config, type='pw')
        self.structure = None
    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary
            Initializes structure as well"""
        (self.namelists, self.cards) = self.parser.parse()
        #self.structure = QEStructure(self)


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 20, 2009 12:21:58 PM$"
