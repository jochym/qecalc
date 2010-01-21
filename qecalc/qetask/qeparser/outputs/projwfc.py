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
# To change this template, choose Tools | Templates
# and open the template in the editor.

from baseoutput import BaseOutput

class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        self.parsers = {}


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Jan 20, 2010 5:31:02 PM$"
