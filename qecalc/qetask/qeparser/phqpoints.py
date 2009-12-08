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
class PHQpoints(object):
    def __init__(self, qeInput):
        self.qeInput = qeInput
        self.isAutomatic = False
        self.grid = None
        self.coords = None

    def set(self, qpoint):
        self.isAutomatic = False
        self.coords = qpoint
        self.grid = None
        self.qeInput.namelist('inputph').remove('nq1')
        self.qeInput.namelist('inputph').remove('nq2')
        self.qeInput.namelist('inputph').remove('nq3')
        self.qeInput.namelist('inputph').remove('ldisp')
        self.qeInput.attach = \
                          "%f    %f    %f\n" % (qpoint[0], qpoint[1], qpoint[2])

    def setAutomatic(self, grid):
        """
        grid - list with dimensions of q-point grid
        """
        self.grid = numpy.array(grid)
        self.isAutomatic = True
        self.coords = None
        self.grid = numpy.array(grid)
        self.qeInput.attach = ''
        self.qeInput.namelist('inputph').add('nq1', str(grid[0]))
        self.qeInput.namelist('inputph').add('nq2', str(grid[1]))
        self.qeInput.namelist('inputph').add('nq3', str(grid[2]))
        self.qeInput.namelist('inputph').add('ldisp', '.true.')


    def parse(self):
        self.coords = []
        if 'nq1' in self.qeInput.namelist('inputph').params:
            self.isAutomatic = True
            self.grid = [int(self.qeInput.namelist('inputph').param('nq1')),
                            int(self.qeInput.namelist('inputph').param('nq2')),
                            int(self.qeInput.namelist('inputph').param('nq3'))]
        else:
            # create a list and get rid of empty lines
            qStrList = [line for line in self.qeInput.attach.split('\n') \
                                                          if not line.isspace()]
            for line in qStrList[1:]:
                points = [float(w) for w in line.split()]
                if len(points) == 3:
                    self.coords.append(points)
                else:
                    raise Exception('Wrong format for q-points')
        self.coords = numpy.array(self.coords)


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 8, 2009 12:22:39 PM$"
