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
class D3Qpoints(object):
    def __init__(self, qeInput):
        self.qeInput = qeInput
        self.coords = None

    def set(self, qpoint):
        qpoint = numpy.array(qpoint)
        self.coords = qpoint
        self.qeInput.attach = \
                          "%f    %f    %f\n" % (qpoint[0], qpoint[1], qpoint[2])

    def parse(self):
        self.coords = []
        # create a list and get rid of empty lines
        qStrList = [line for line in self.qeInput.attach.split('\n') \
                                    if len(line.split()) == 3]
        for line in qStrList[:]:
            points = [float(w) for w in line.split()]
            if len(points) == 3:
                self.coords.append(points)
            else:
                raise Exception('Wrong format for q-points')
        self.coords = numpy.array(self.coords)


if __name__ == "__main__": pass

__author__="Nikolay Markovskiy"
__date__ ="$Dec 8, 2009 5:57:48 PM$"
