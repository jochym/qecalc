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
class PWKpoints(object):
    def __init__(self, qeInput):
        self.qeInput = qeInput
        self.isAutomatic = False

    def parse(self):
        """ Returns array of points. Does not have to be AUTOMATIC """
        kpoints = []
        if 'k_points' not in self.qeInput.cards:
            kpoints = numpy.array(kpoints)
            self.shifts = None
            return
        if self.qeInput.cards['k_points'].arg() == 'automatic':
            self.isAutomatic = True
        for line in self.qeInput.cards['k_points'].lines():
            if self.isAutomatic:
                kpoints.append([int(k) for k in line.split()])
                if len(kpoints[0][:]) > 3:
                    shifts = [int(k) for k in kpoints[0][3:]]
                else:
                    shifts = [int(0) ,int(0), int(0)]
                break
            else:
                kpoints.append([float(k) for k in line.split()])
        
        if self.isAutomatic:
            shifts = numpy.array(shifts)            
        else:
            # first line was number of kpoints
            kpoints.pop(0)
            shifts = None
        kpoints = numpy.array(kpoints)
        if self.isAutomatic:
            self.coords = None
            self.grid = kpoints
            self.shifts = shifts
            self.weights = None
        else:
            if kpoints.shape[1] == 4:
                self.coords = kpoints[:,:3]
                self.weights = kpoints[:,3].flatten()
            else:
                self.coords = kpoints
                self.weights = numpy.zeros(kpoints.shape[0]) + 1.0 
            
            self.shifts = None
        
        
    def set(self, points, weights = None, type = 'triba'):
        """
        Sets k-points from numpy.array of size (n,3)
        Default values of array of weights are 1.0 
        """
        
        typesPossible = ['triba', 'crystal']
        if type not in typesPossible:            
            raise Exception('Wrong type of k-point grid')
        
        points = numpy.array(points)
        self.isAutomatic = False
        self.qeInput.cards['k_points'].removeLines()
        self.qeInput.cards['k_points'].setArg(type)
        
        self.coords = points
        self.grid = None
        # weights == None for custom generated q-point grid
        if weights == None:
            weights = numpy.zeros(points.shape[0])
            weights = weights + 1.0
        self.weights = weights
        string = str(points.shape[0]) + '\n'
        for qpoint, coord in zip(self.coords, self.weights):
            string = string + \
                   "%# .8f %# .8f %# .8f %# .8f\n" % (qpoint[0], qpoint[1], qpoint[2], coord)
        
        self.qeInput.cards['k_points'].addLine(string)


    def setAutomatic(self, points, shifts = None):
        """Takes two lists of values"""
        self.isAutomatic = True
        if shifts == None:
            if len(points) == 3:
                shifts = numpy.array([int(0) ,int(0), int(0)])
            if len(points) == 6:
                shifts = points[3:]
                points = points[:3]
        else:
            if len(shifts) == 3 and len(points) == 3:
                shifts = shifts
                points = points
            else:
                raise Exception('Wrong number of kpoints and/or shifts')

        self.grid = numpy.array(points)
        self.shifts = numpy.array(shifts)
        self.qeInput.cards['k_points'].removeLines()
        self.qeInput.cards['k_points'].setArg('AUTOMATIC')
        string = ""
        for k in self.grid:
            string = string + str(k) + " "
        for s in self.shifts:
            string = string + str(s) + " "
        self.qeInput.cards['k_points'].addLine(string)

if __name__ == "__main__": pass

__author__="Nikolay Markovskiy"
__date__ ="$Oct 20, 2009 12:59:23 PM$"
