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

from qesinput import QESInput
from phqpoints import PHQpoints

class PHInput(QESInput):
    def __init__(self, filename=None, config=None, setting = None):
        QESInput.__init__(self,filename, config, type='ph', setting = setting)
        self.qpoints = PHQpoints(self)


    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary
            Initializes structure as well"""
        #(self.header, self.namelists, self.cards, self.attach) = self.parser.parse()
        QESInput.parse(self)
        self.qpoints.parse()


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 8, 2009 12:39:52 PM$"
