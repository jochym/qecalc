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
import numpy
class PWKpoints(object):
    def __init__(self, qeInput):
        self.qeInput = qeInput
        self.isAutomatic = False
        self.getkPointsFromPWSCF()

    def getkPointsFromPWSCF(self):
        """ Returns array of points. Does not have to be AUTOMATIC """
        kpoints = []
        if self.qeInput.cards['k_points'].arg() == 'automatic':
            self.isAutomatic = True
        for line in self.qeInput.cards['k_points'].lines():
            if self.isAutomatic:
                kpoints.append([int(k) for k in line.split()])
                break
            else:
                kpoints.append([float(k) for k in line.split()])
        kpoints = numpy.array(kpoints)
        if self.qeInput.cards['k_points'].arg() == 'automatic':
            self.isAutomatic = True
            if len(kpoints[0,:]) > 3:
                shifts = [int(k) for k in kpoints[0,3:]]
            else:
                self.shifts = numpy.array([0,0,0])

        return numpy.array(kpoints)

    def setAutomatic(self, kpoints, shifts = None):
        if shifts == None:
            shifts = [int(k) for k in self.getkPointsFromPWSCF()[0,3:]]
        self.qeInput.cards['k_points'].removeLines()
        self.qeInput.cards['k_points'].setArg('AUTOMATIC')
        string = ""
        for k in kpoints:
            string = string + str(k) + " "
        for s in shifts:
            string = string + str(s) + " "
        self.qeInput.cards['k_points'].addLine(string)
        self.qeInput.save(self.pwscfInput)

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 20, 2009 12:59:23 PM$"
