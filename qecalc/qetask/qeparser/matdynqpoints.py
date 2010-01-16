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
class MatdynQpoints(object):
    def __init__(self, qeInput):
        self.qeInput = qeInput
        self.isAutomatic = False
        self.grid = None
        self.coords = None
        self.axis = None

    def set(self, qpoints, axis = None):
        qpoints = numpy.array(qpoints)
        self.isAutomatic = False
        self.qeInput.namelist('input').remove('nk1')
        self.qeInput.namelist('input').remove('nk2')
        self.qeInput.namelist('input').remove('nk3')
        self.qeInput.namelist('input').remove('fldos')
        self.qeInput.namelist('input').remove('dos')
        self.coords = qpoints
        self.grid = None
        # axis == None for custom generated q-point grid
        if axis == None:
            axis = numpy.zeros(qpoints.shape[0])
        self.axis = axis
        string = str(qpoints.shape[0]) + '\n'
        for qpoint, coord in zip(self.coords, self.axis):
            string = string + \
                   "%f    %f    %f    %f\n" % (qpoint[0], qpoint[1], qpoint[2], coord)

        self.qeInput.attach = string

    def setAutomatic(self, grid):
        """
        grid - list with dimensions of q-point grid
        """
        self.isAutomatic = True
        self.coords = None
        self.grid = numpy.array(grid)
        self.qeInput.attach = ''
        self.qeInput.namelist('input').add('nk1', str(grid[0]))
        self.qeInput.namelist('input').add('nk2', str(grid[1]))
        self.qeInput.namelist('input').add('nk3', str(grid[2]))
        self.qeInput.namelist('input').add('dos', '.true.')


    def parse(self):
        self.coords = []
        self.axis = []
        if 'nk1' in self.qeInput.namelist('input').params:
            self.isAutomatic = True
            self.grid = [int(self.qeInput.namelist('input').param('nk1')),
                            int(self.qeInput.namelist('input').param('nk2')),
                            int(self.qeInput.namelist('input').param('nk3'))]
            self.grid = numpy.array(self.grid)
        else:
            # create a list and get rid of empty lines
            qStrList = [line for line in self.qeInput.attach.split('\n') \
                                                          if not line.isspace() and line != '']
            for line in qStrList[1:]:
                points = [float(w) for w in line.split()]
                if len(points) == 3:
                    self.coords.append(points)
                else:
                    if len(points) == 4:
                        self.coords.append(points[0:3])
                        self.axis.append(points[3])
                    else:
                        raise Exception('Wrong format for q-points')
        self.coords = numpy.array(self.coords)
        self.axis = numpy.array(self.axis)

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 22, 2009 12:47:27 PM$"
