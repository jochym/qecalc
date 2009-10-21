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

from qetask import QETask

class TaskMerger(QETask):
    def __init__(self, setting, tasks, cleanOutDir = False):
        QETask.__init__(self, setting, cleanOutDir)
        self.tasks = tasks
        self.cmdStr = tasks[0].cmdLine()
        self.name = tasks[0].name
        for task in tasks[1:]:
            self.cmdStr = self.cmdStr + ' ; ' + task.cmdLine()
            self.name = self.name + ' -> ' + self.cmdStr            
    
    def launch(self, cleanOutDir = None):
        if cleanOutDir != None:
            clean = cleanOutDir
        else:
            clean = self.cleanOutDir
        if clean:
            self.cleanOutputDir()
        for task in self.tasks:
            task.input.parse()
        self._run()
        for task in self.tasks:
            task.output.parse(parserList = 'all')