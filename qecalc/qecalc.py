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
from setting import Setting
import numpy

class QECalc(object):
    def __init__(self, fname):
        
        self.setting = Setting(fname)
        self.taskList = []

    def launch(self):
        for task in self.taskList:
            task.launch()

    def lookupProperty(self, propertyName, taskList = None, \
                                                           withUnits = False):
        """Will look up a specific output property from a list if tasks"""
        if taskList == None:
            taskList = self.taskList
        value = None
        for task in taskList:
            try:
                value = task.output.property(propertyName, withUnits)
            except KeyError:
                pass
        return value

if __name__ == '__main__':
    qe = QECalc('config.ini')
    qe.qeConfig.setNamelistParameter('system', 'ecutwfc', 44)
    qe.qeConfig.save(qe.qeConfig.filename)
    print qe.getkPoints().shape
    qe.setkPointsAutomatic(numpy.array([17, 33, 12, 0, 0, 0]))
    print 'ibrav' in qe.qeConfig.namelists['system'].params
