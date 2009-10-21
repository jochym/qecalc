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

if __name__ == '__main__':
    qe = QECalc('config.ini')
    qe.qeConfig.setNamelistParameter('system', 'ecutwfc', 44)
    qe.qeConfig.save(qe.qeConfig.filename)
    print qe.getkPoints().shape
    qe.setkPointsAutomatic(numpy.array([17, 33, 12, 0, 0, 0]))
    print 'ibrav' in qe.qeConfig.namelists['system'].params
