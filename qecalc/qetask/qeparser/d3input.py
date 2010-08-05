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
from d3qpoints import D3Qpoints

class D3Input(QEInput):
    def __init__(self, filename=None, config=None):
        QEInput.__init__(self,filename, config, type='d3')
        self.qpoints = D3Qpoints(self)


    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary
            Initializes structure as well"""
        (self.header, self.namelists, self.cards, self.attach) = self.parser.parse()
        self.qpoints.parse()


if __name__ == "__main__": pass

__author__="Nikolay Markovskiy"
__date__ ="$Dec 8, 2009 6:00:01 PM$"
