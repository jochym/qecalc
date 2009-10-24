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
    """Base abstract class for 'CalcNameCalc' family of classes and does not provide
    any functionality. All user defined classes should be derived from this
    class. All tasks of user defined Calcs should be appended to the taskList
    for launch() and lookupProperty methods to work

      setting -- provides access to parallel environment and QE input/output files

      taskList -- list of all user specified tasks. Should reflect
      their launching order
    """
    def __init__(self, fname):
        
        self.setting = Setting(fname)
        self.taskList = []

    def launch(self):
        """Launches all tasks in taskList one after another:
        
          >>> pwCalc = PWCalc('config.ini')
          >>> pwCalc.launch()
        
        """
        for task in self.taskList:
            task.launch()

    def lookupProperty(self, propertyName, taskList = None, \
                                                           withUnits = False):
        """
        Will look up a specific output property with propertyName from the
        list of tasks taskList:
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
