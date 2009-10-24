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

from qecalc import QECalc
import numpy

from qetask.pwtask import PWTask

class PWCalc(QECalc):
    """ Calc for single pw.x use.:
        Task list:
          pw -- PWTask
          Example:
            >>> pwCalc = PWCalc('config.ini')
            >>> pwCalc.launch()
    """
    def __init__(self, filename):
        QECalc.__init__(self, filename)
        self.pw = PWTask(self.setting)
        self.taskList = [self.pw]


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 3:05:10 PM$"
