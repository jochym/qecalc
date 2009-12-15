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
from setting import Setting
from qetorque import QETorque

class QETask(object):
    """
    Base abstract class for 'TaskNameTask' family of classes. It also provides
    the interface to loumch tasks

    All user defined tasks should be derived from this class

    setting -- provides access to parallel environment and QE input/output files
    """
    def __init__(self, filename = 'config.ini', configString = None, cleanOutDir = None):
       # parallelization parameters
       # Default values, see explanations below:
        #self.name = 'Launcher'
        configDic = {
        'useTorque' : 'False',
        'torqueResourceList': '-l nodes=1:ppn=1',
        'paraPrefix': '',
        'paraPostfix': '',
        'paraRemoteShell': '',
        'outDir': ''
        }
        self.setting = Setting(filename, configString)
        self.setting.section(QETask.name(self), configDic)

        self.cleanOutDir = cleanOutDir

        if self.setting.useTorque == 'True':
            self.setting.useTorque = True
        else:
            self.setting.useTorque = False

        if self.setting.useTorque:
            self._torque = QETorque(self.setting.torqueResourceList)


    def name(self):
        return 'Launcher'


    def _check(self, x):
        """
        Will check the exit status of the program to be executed
        """
        signal = x & 0xFF
        exitcode = (x >> 8) & 0xFF
        if exitcode != 0:
            raise Exception("Task " + self.name() + " crashed: check your settings" + "Command string:" + self.cmdLine())

    def _run(self):
        if os.path.exists('CRASH'):
            os.remove('CRASH')
        os.system(self.setting.paraRemoteShell + ' mkdir -p ' + self.setting.outDir)
        if self.setting.paraPrefix != '' and self.setting.paraPrefix in self.cmdLine():
            if self.setting.useTorque:
                self._torque.serial(self.cmdLine())
            else:
                self._check(os.system(self.cmdLine()))
        else:
            self._check(os.system(self.cmdLine()))
        if os.path.exists('CRASH'):
            raise Exception("Task " + self.name() + " crashed: 'CRASH' file was discovered")


    def _syncSetting(self):
        """
        Will syncronise QE input file with class Setting for given task (QE input
        files may contain multiple input/ouput file names  definitions which
        can be overriden in this method)
        """
        pass


    def cleanOutputDir(self, cleanOutDir = None):
        """
        Cleans the output directory (directory, where large files, used
        for calculation are stored, can be for example  'temp/' or 'scratch/')
        """
        import shutil
        if cleanOutDir == None:
            cleanOutDir = self.cleanOutDir
        if cleanOutDir == None:
            raise Exception('outDir can not be cleaned, since it was not defined')

        os.system(self.setting.paraRemoteShell + ' rm -r -f ' + cleanOutDir)
        os.system(self.setting.paraRemoteShell + ' mkdir ' + cleanOutDir)
#        if self.setting.useTorque:
#
#        else:
#            shutil.rmtree(cleanOutDir)
#            os.mkdir(cleanOutDir)
    
#    def cmdLine(self):
#        """
#        returns command string of a given task
#        """
#        return self._cmdStr

    def launch(self, cleanOutDir = None):
        """
        Parses input. Launches task. Parses output.
        """
        if cleanOutDir == None:
            cleanOutDir = self.cleanOutDir
        if cleanOutDir != None:
            self.cleanOutputDir(cleanOutDir)
        self.input.parse()
        self._syncSetting() # sync setting with QE input file
        self.input.save()
        self._run()
        self.output.parse(parserList = 'all')

__author__="kolya"
__date__ ="$Oct 18, 2009 5:03:21 PM$"

if __name__ == "__main__":
    print "Hello World";
