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
from qeparser.qeinput import QEInput
from qeparser.qeoutput import QEOutput

class Q2RTask(QETask):
    def __init__(self, filename = None,configString = None, cleanOutDir = None):
        QETask.__init__(self, filename, configString, cleanOutDir)

        self.name = 'q2r.x'

        configDic = {
        'q2rInput': 'q2r.in',
        'q2rOutput': 'q2r.out'
        }
        self.setting.section(self.name, configDic)

        self.input = QEInput(filename = self.setting.q2rInput, type = 'q2r')
        self.output = QEOutput(self.setting, type = 'q2r')
        self._cmdStr = "q2r.x < " + self.setting.q2rInput + " > " + \
                        self.setting.q2rOutput


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 3:12:40 PM$"
