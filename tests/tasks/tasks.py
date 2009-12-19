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
# To change this template, choose Tools | Templates
# and open the template in the editor.

from qecalc.qetask.pwtask import PWTask
from qecalc.qetask.phtask import PHTask
from qecalc.qetask.matdyntask import MatdynTask
from qecalc.qetask.dynmattask import DynamatTask
from qecalc.qetask.d3task import D3Task
from qecalc.qetask.q2rtask import Q2RTask
from qecalc.qetask.dostask import DOSTask

def testTask(task):
    print 'Testing task ', task.name(), ':'
    task.launch()
    print 'Printing properties:'
    task.output.properties()


if __name__ == "__main__":
    pw      = PWTask()
    ph_sp   = PWTask()
    dynmat  = PWTask()
    q2r     = PWTask()
    matdyn  = PWTask()
    d3      = PWTask()
    dos     = DOSTask()

    taskList = [pw, ph, dynmat]

    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 18, 2009 3:03:02 PM$"
