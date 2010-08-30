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

from qetask import QETask

class TaskMerger(QETask):
    """
    Takes list of tasks as an argument
    Merges their names and command lines. Can be manipulated as a single task
    when launching.
    When torque is used, merged task is submitted as a single batch job
    """
    def __init__(self, tasks, cleanOutDir = False, ioTask = None):
        QETask.__init__(self, filename = tasks[0].setting.get('filename'),  \
                        configString = tasks[0].setting.get('configString'),\
                        cleanOutDir = cleanOutDir)

        self.setParallel()
        #self.readSetting(self, filename = tasks[0].setting.get('filename'), \
        #                 configString = tasks[0].setting.get('configString') )
        self._mergedTask = True

        if ioTask == None:
            ioTask = tasks[0]
        self.input = ioTask.input
        self.output = ioTask.output
        self.setting = ioTask.setting

        self.tasks = tasks
#        self._cmdStr = tasks[0].cmdLine()
#        self.name = tasks[0].name
#        for task in tasks[1:]:
#            self._cmdStr = self._cmdStr + ' ; ' + task.cmdLine()
#            self.name = self.name + ' -> ' + self._cmdStr


    def cmdLine(self):
        cmd = self.tasks[0].cmdLine()
        for task in self.tasks[1:]:
            cmd = cmd + ' ; ' + task.cmdLine()
        return cmd


    def name(self):
        name = self.tasks[0].name()
        for task in self.tasks[1:]:
            name = name + ' -> ' + task.name()
        return name

    def launch(self, cleanOutDir = None):

        for task in self.tasks:
            task.input.parse()
            task.syncSetting()
            task.input.save()

#        if cleanOutDir != None:
#            clean = cleanOutDir
#        else:
#            clean = self.cleanOutDir
#        if clean:
#            self.cleanOutputDir()
        self.cleanOutputDir(cleanOutDir)
        self._run()
        
        for task in self.tasks:
            task.output.parse(parserList = 'all')

    def syncSetting(self):
        """
        When this method is called on launch(), the input file is already
        parsed and will be saved before the run...
        """
        for task in self.tasks:
            task.input.parse()
            task.syncSetting()
            