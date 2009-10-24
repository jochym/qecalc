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
from setting import Setting
import numpy

class QECalc(object):
    """Base abstract class for 'CalcNameCalc' family of classes.
    All user defined classes should be derived from this
    class. All tasks of user defined Calcs should be appended to the taskList
    for launch() and lookupProperty() methods to work properly

      setting -- provides access to parallel environment and QE input/output files

      taskList -- list of all user specified tasks. Should reflect
      their launching order
    """
    def __init__(self, filename):
        
        self.setting = Setting(filename)
        self.taskList = []

    def launch(self, taskList = None):
        """Launches all tasks in taskList one after another. With no taskList
        provided internal taskList is used.

        Example:
        
          >>> pwCalc = PWCalc('config.ini')
          >>> pwCalc.launch()
        
        """
        if taskList == None:
            taskList = self.taskList
            
        for task in taskList:
            task.launch()

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
            value = task.output.property(propertyName, withUnits)
            #except KeyError:
            #    pass
        return value

if __name__ == '__main__':
    qe = QECalc('config.ini')
    qe.qeConfig.setNamelistParameter('system', 'ecutwfc', 44)
    qe.qeConfig.save(qe.qeConfig.filename)
    print qe.getkPoints().shape
    qe.setkPointsAutomatic(numpy.array([17, 33, 12, 0, 0, 0]))
    print 'ibrav' in qe.qeConfig.namelists['system'].params
