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

import os
from qetorque import QETorque

class QETask(object):
    def __init__(self, setting, cleanOutDir = False):
        self.setting = setting
        self.cleanOutDir = cleanOutDir
        if self.setting.useTorque:
            self.torque = QETorque(self.setting.configFileName)        
        self.input = None
        self.output = None
        self.cmdStr = None
        self.name = None

    def _check(self, x):
        """Will check the exit status of the program to be executed"""
        signal = x & 0xFF
        exitcode = (x >> 8) & 0xFF
        if exitcode != 0:
            raise Exception("Task " + self.name + " crashed: check your settings" + "Command string:" + self.cmdStr)

    def _run(self):
        os.remove('CRASH')
        if self.setting.paraPrefix != '' and self.setting.paraPrefix in self.cmdStr:
            if self.setting.useTorque:
                self.torque.serial(self.cmdStr)
            else:
                self._check(os.system(self.cmdStr))
        else:
            self._check(os.system(self.cmdStr))
        if os.path.exists('CRASH'):
            raise Exception("Task " + self.name + " crashed: 'CRASH' file was discovered")


    def cleanOutputDir(self):
        """
        Parses pwscfInput file for output dir name and cleans the directory
        """
        from qeparser.qeinput import QEInput
        import shutil
        qeConf = QEInput(self.setting.pwscfInput)
        qeConf.parse()
        outDir = qeConf.namelist('control').param('outdir')[1:-1]
        if self.setting.useTorque:
            os.system('bpsh -a rm -r -f ' + outDir)
            os.system('bpsh -a mkdir ' + outDir)
        else:
            shutil.rmtree(outDir)
            os.mkdir(outDir)        
    
    def cmdLine(self):
        return self.cmdStr

    def launch(self, cleanOutDir = None):
        """
        Parses input. Launches task. Parses output.
        """
        if cleanOutDir != None:
            clean = cleanOutDir
        else:
            clean = self.cleanOutDir
        if clean:
            self.cleanOutputDir()
        self.input.parse()
        self._run()
        self.output.parse(parserList = 'all')

__author__="kolya"
__date__ ="$Oct 18, 2009 5:03:21 PM$"

if __name__ == "__main__":
    print "Hello World";
