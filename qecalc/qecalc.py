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

class QECalc(object):
    """Base abstract class for 'CalcNameCalc' family of classes.
    All user defined classes should be derived from this
    class. All tasks of user defined Calcs should be appended to the taskList
    for launch() and lookupProperty() methods to work properly

      taskList -- list of all user specified tasks. Should reflect
      their launching order
    """
    def __init__(self):
        #self.setting = Setting(filename)
        self.taskList = []

    def launch(self, taskList = None):
        """Launches all tasks in taskList one after another. With no taskList
        provided internal taskList is used.

        Example:
        
          >>> pwCalc = PWCalc('config.ini')
          >>> pwCalc.launch()
        
        """
        #self.syncInputs()
        if taskList == None:
            taskList = self.taskList
            
        for task in taskList:
            task.launch()

    def getAllTasks(self):

        taskList = []
        for task in self.taskList:
            if task.isMerged():
                for subTask in task.tasks:
                    taskList.append(subTask)
            else:
                taskList.append(task)

        return taskList

    def syncInputs(self): pass

    def lookupProperty(self, propertyName, taskList = None, \
                                                           withUnits = False):
        """
        Will look up a specific output property with propertyName in the
        list of tasks taskList. With no taskList  provided internal taskList
        is used.

        Example:

          >>> phonCalc = SinglePhononCalc('config.ini')
          >>> phonCalc.launch()
          >>> print phonCalc.lookupProperty('total energy', withUnits = True)
          >>> print phonCalc.lookupProperty('single phonon', withUnits = True)
        """
        if taskList == None:
            taskList = self.taskList
            
        value = None
        for task in taskList:
            #try:
            value = task.output.property(propertyName, withUnits = False)
            if value != None:
                value = task.output.property(propertyName, withUnits)
                break
            #except KeyError:
            #    pass
        return value

    def _populateTasks(self, filename = None, configString = None, sectionList = None, taskList = None):
        self.taskList = []
        if sectionList != None:
            if len(sectionList) == len(self._taskSpec):
                for i, task in enumerate(self._taskSpec):
                    setattr(self, task[0],task[1](filename = filename, configString = configString, sectionName = sectionList[i]))
                    self.taskList.append(getattr(self, task[0]))
                return
            else:
                raise NameError('sectionList and _taskSpec dimensions do not match')

        if taskList != None:
            if len(taskList) == len(self._taskSpec):
                for i, task in enumerate(self._taskSpec):
                    if task[0] in taskList[i].name():
                        setattr(self, task[0], taskList[i])
                        self.taskList.append( getattr(self, task[0]))
                    else:
                        raise NameError('task names in taskList and _taskSpec  do not match')
                return
            else:
                raise NameError('taskList and _taskSpec dimensions do not match')

        for i, task in enumerate(self._taskSpec):            
            setattr(self, task[0], task[1](filename = filename, configString = configString))
            self.taskList.append(getattr(self, task[0]))


if __name__ == '__main__':
    pass
