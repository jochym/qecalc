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

        self._mergedTask = False
        
        self.cleanOutDir = cleanOutDir

        QETask.readSetting(self, filename, configString)


    def readSetting(self, filename = None, configString = None):
        """
        Initializes setting class
        """
        
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
            configString = ''
       #     filename = 'config.ini'

        self.setting = Setting(filename, configString)
        self.setting.section(QETask.name(self), configDic)
        

        if self.setting.get('useTorque') == 'True':
            self.setting.set('useTorque', True)
        else:
            self.setting.set('useTorque', False)

        if self.setting.get('useTorque'):
            if self.isSerial():
                self._torque = QETorque(self.setting.get('serialTorqueParams'))
            else:
                self._torque = QETorque(self.setting.get('paraTorqueParams'))
        
        #self.syncSetting()      


    def name(self):
        return 'Launcher'


    def parse(self):
        
        self.syncSetting()


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
        #self.input.parse()
        self.syncSetting() # sync setting with QE input file
        self.input.save()
        self.cleanOutputDir(cleanOutDir)
        self._run()
        self.output.parse(parserList = 'all')

    def setSerial(self):
        #self._serialThroughMPI = serialThroughMPI
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
 #           if self._isSerialThroughMPI():
 #               return ''
 #           else:
            return self.setting.get('serialPrefix')
        else:
            return self.setting.get('paraPrefix')


    def getPostfix(self):
        if self.isSerial():
            return self.setting.get('serialPostfix')
        else:
            return self.setting.get('paraPostfix')


    def isMerged(self):
        return self._mergedTask


    def _initCmdLineParams(self):
        self._inp = '-inp'
        self._devnul = '< /dev/null'
        if self.isSerial():
            self._inp = '<'
            self._devnul = ''


#    def _isSerialThroughMPI(self):
#        if self.isSerial():
#            if self._serialThroughMPI:
#                return True
#            else:
#                return False


    def _check(self, x):
        """
        Will check the exit status of the program to be executed
        """
        signal = x & 0xFF
        exitcode = (x >> 8) & 0xFF
        if exitcode != 0:
            raise Exception("Task " + self.name() + " crashed: check your \
settings!\n" + "Command string: " + self.cmdLine())

    def _run(self):
        if os.path.exists('CRASH'):
            os.remove('CRASH')

        outdir = self.setting.get('outdir')
        if outdir != None:
            os.system(self.setting.paraRemoteShell + ' mkdir -p ' + outdir)

        if self.setting.get('useTorque'):
            if self.isSerial():
                if self.setting.get('serialPrefix') == '':
                    self._check(os.system(self.cmdLine()))
                else:
                    self._torque.serial(self.cmdLine(), torqueParams = self.setting.serialTorqueParams)
            else:
                self._torque.serial(self.cmdLine(), torqueParams = self.setting.paraTorqueParams)
        else:
            self._check(os.system(self.cmdLine()))

#        if self.setting.paraPrefix != '' and self.setting.paraPrefix in self.cmdLine():
#            if self.setting.useTorque:
#                self._torque.serial(self.cmdLine())
#            else:
#                self._check(os.system(self.cmdLine()))
#        else:
#            self._check(os.system(self.cmdLine()))
        if os.path.exists('CRASH'):
            raise Exception("Task " + self.name() + " crashed: 'CRASH' file was discovered")

    def _getCmdLine(self, executable, input, output):
        self._initCmdLineParams()
        return  self.getPrefix() + ' ' + executable + " " + self.getPostfix() + \
                ' ' + self._inp + ' ' +  self.setting.get(input) + ' > ' + \
                self.setting.get(output) + self._devnul
                
                
    def _setSetiingInputOutput(self, configDic, filename = None, configString = None, sectionName = None):
        """
        Helper method for readSetting()
        """
        if sectionName == None:
            name = self.name()
        else:
            name = sectionName

        self.setting.section(name, configDic)
        self.input = PWInput( self.setting.pwInput ) #self.setting.get('pwInput') )
        # add pointer to setting for input filenames synchronization 
        self.input._setting = self.setting
        self.output = QEOutput(self.setting, type='pw')
        
        if filename != None or configString != None:
            self.syncSetting()        

__author__="kolya"
__date__ ="$Oct 18, 2009 5:03:21 PM$"

if __name__ == "__main__": pass
