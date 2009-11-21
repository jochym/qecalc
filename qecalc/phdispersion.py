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

from qedispersion import QEDispersion

class PHDispersion(QEDispersion):
    def __init__(self, lattice , matdynTask):
        QEDispersion.__init__(self,lattice)
        self.matdynTask = matdynTask

    def launch(self, *pathNPoints):
        #self.path, self.axis =
        self.setPath(*pathNPoints)

        self.matdynTask.input.qpoints.set(self.path, self.axis)
        self.matdynTask.input.save()

        self.matdynTask.launch()
        pol, disp, qpoints =  self.matdynTask.output.property('multi phonon')
        self.dispersion = disp

    def setPhononPath(self, *pathNPoints):
        """
        Obsolete
        """
        self.launch(*pathNPoints)


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 1:18:20 PM$"
