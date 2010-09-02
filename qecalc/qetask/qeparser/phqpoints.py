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

    def set(self, point):
        point = numpy.array(point)
        self.isAutomatic = False
        self.coords = point
        self.grid = None
        self.qeInput.namelist('inputph').remove('nq1')
        self.qeInput.namelist('inputph').remove('nq2')
        self.qeInput.namelist('inputph').remove('nq3')
        self.qeInput.namelist('inputph').remove('ldisp')
        self.qeInput.attach = \
                          "%f    %f    %f\n" % (point[0], point[1], point[2])

    def setAutomatic(self, grid):
        """
        grid - list with dimensions of q-point grid
        """
        self.grid = numpy.array(grid)
        self.isAutomatic = True
        self.coords = None
        self.grid = numpy.array(grid)
        self.qeInput.attach = ''
        self.qeInput.namelist('inputph').set('nq1', str(grid[0]))
        self.qeInput.namelist('inputph').set('nq2', str(grid[1]))
        self.qeInput.namelist('inputph').set('nq3', str(grid[2]))
        self.qeInput.namelist('inputph').set('ldisp', '.true.')


    def parse(self):
        self.coords = []
        if 'nq1' in self.qeInput.namelist('inputph').paramlist():
            self.isAutomatic = True
            self.grid = [int(self.qeInput.namelist('inputph').get('nq1')),
                            int(self.qeInput.namelist('inputph').get('nq2')),
                            int(self.qeInput.namelist('inputph').get('nq3'))]
            self.grid = numpy.array(self.grid)
        else:
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
__date__ ="$Dec 8, 2009 12:22:39 PM$"
