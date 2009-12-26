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
from qetask.phtask import PHTask
from qetask.d3task import D3Task
from qetask.pwphmerger import PWPHMerger

class D3Calc(QECalc):
    """ Calc for multi phonon calculations:

    Task list:

      pw     -- PWTask

      ph     -- PHTask

      d3     -- D3Task

      pwph   -- PWPHMerger - task merger used for submission of pw.x and ph.x
      commands in a single command string

      taskList = [pwph, d3]

    Example:

      >>> phonCalc = SinglePhononCalc('config.ini')
      >>> phonCalc.launch()
      >>> print phonCalc.dynmat.output('single phonon')

    """
    def __init__(self, filename = None, configString = None, sectionList = None, taskList = None):
        QECalc.__init__(self)

        # tasks definition:
        # specify taskName/taskConstructur pairs
        self._taskSpec = [
                          ['pw', PWTask],
                          ['ph', PHTask],
                          ['d3', D3Task]
                         ]

        self._populateTasks(filename, configString, sectionList, taskList)
        #self.pw = PWTask(filename)
        #self.ph = PHTask(filename)
        #self.d3 = D3Task(filename)
        self.pwph = PWPHMerger(self.pw,self.ph, cleanOutDir = True)
        self.taskList = [self.pwph, self.d3]

        #Hack: make sure d3 task is serial (d3 does not seem to work in parallel)
        for task in [self.pw, self.ph, self.d3]:
            #task.setting.useTorque = 'False'
            task.setSerial()
            task.setting.serialPrefix = ''


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 8, 2009 10:45:30 AM$"
