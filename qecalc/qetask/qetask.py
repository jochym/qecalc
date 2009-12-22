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
    def __init__(self, filename = None, configString = None, cleanOutDir = False):
       # parallelization parameters
       # Default values, see explanations below:
        #self.name = 'Launcher'

        self._isSerial = True

        configDic = {
        'useTorque' : 'False',
        'paraTorqueParams':  '-l nodes=1:ppn=1',
        'serialTorqueParams': '-l nodes=1:ppn=1',        
        'paraPrefix': '',
        'paraPostfix': '',
        'serialPrefix': '',
        'serialPostfix': '',
        'paraRemoteShell': '',
        'outdir': None
        }

        if filename == None and configString == None:
            filename = 'config.ini'

        self.setting = Setting(filename, configString)
        self.setting.section(QETask.name(self), configDic)

        self.cleanOutDir = cleanOutDir

        if self.setting.get('useTorque') == 'True':
            self.setting.set('useTorque', True)
        else:
            self.setting.set('useTorque', False)

        if self.setting.get('useTorque'):
            self._torque = QETorque(self.setting.get('torqueResourceList'))


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

        outdir = self.setting.get('outdir')
        if outdir != None:
            os.system(self.setting.paraRemoteShell + ' mkdir -p ' + outdir)
        if self.setting.paraPrefix != '' and self.setting.paraPrefix in self.cmdLine():
            if self.setting.useTorque:
                self._torque.serial(self.cmdLine())
            else:
                self._check(os.system(self.cmdLine()))
        else:
            self._check(os.system(self.cmdLine()))
        if os.path.exists('CRASH'):
            raise Exception("Task " + self.name() + " crashed: 'CRASH' file was discovered")


    def syncSetting(self):
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
#        import shutil
        if cleanOutDir == None:
            cleanOutDir = self.cleanOutDir
#        if cleanOutDir == None:
#            raise Exception('outDir can not be cleaned, since it was not defined')



        if cleanOutDir:
            outDir = self.setting.get('outdir')
            os.system(self.setting.paraRemoteShell + ' rm -r -f ' + outDir)
            os.system(self.setting.paraRemoteShell + ' mkdir ' + outDir)
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

    def launch(self, cleanOutDir = False):
        """
        Parses input. Launches task. Parses output.
        """
#        if cleanOutDir == None:
#            cleanOutDir = self.cleanOutDir
#        if cleanOutDir != None:
#            self.cleanOutputDir(cleanOutDir)
        self.input.parse()
        self.syncSetting() # sync setting with QE input file
        self.input.save()
        self.cleanOutputDir(cleanOutDir)
        self._run()
        self.output.parse(parserList = 'all')

    def setSerial(self):
        self._isSerial = True


    def setParallel(self):
        self._isSerial = False

    def isSerial(self):
        if self._isSerial:
            return True
        else:
            return False

    def getPrefix(self) :
        if self.isSerial():
            return self.setting.get('serialPrefix')
        else:
            return self.setting.get('paraPrefix')


    def getPostfix(self):
        if self.isSerial():
            return self.setting.get('serialPostfix')
        else:
            return self.setting.get('paraPostfix')

    def _initCmdLineParams(self):
        self._inp = '-inp'
        self._devnul = '< /dev/null'
        if self.isSerial():
            self._inp = '<'
            self._devnul = ''


    def _getCmdLine(self, executable, input, output):
        self._initCmdLineParams()
        return  self.getPrefix() + ' ' + executable + " " + self.getPostfix() + \
                ' ' + self._inp + ' ' +  self.setting.get(input) + ' > ' + \
                self.setting.get(output) + self._devnul

__author__="kolya"
__date__ ="$Oct 18, 2009 5:03:21 PM$"

if __name__ == "__main__":
    print "Hello World";
