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
import subprocess
import sys
import os
import time
class QETorque:
    def __init__(self, torqueParams):
        # need this in case $HOME includes symbolic links otherwise torque
        # gets confused:
        myHome = os.environ['HOME']
        myRealHome = os.path.realpath(myHome)
        myRealWorkDir = os.getcwd()
        if myRealHome in myRealWorkDir:
            relativeWorkDir = myRealWorkDir.replace(myRealHome, '')
            self._workDir = myHome + relativeWorkDir
        else:
            self._workDir = os.getcwd()
        self._jobID = None

        self.torqueParams = torqueParams

        # use qmgr (not implemented)
    def submit(self, cmdStr, torqueParams = None):
        """Submits job. cmdStr is mpirun + params + program is if one did not
           use torque"""
        if torqueParams == None:
            torqueParams = self.torqueParams
        else:
            self.torqueParams = torqueParams
        if torqueParams == None:
            self.torqueParams = ''

        submitStr = "echo '" + cmdStr + "' | " + 'qsub -V -d ' + \
                     self._workDir + ' ' + self.torqueParams + ' -'
        print "submitStr = " + submitStr
        try:
            p = subprocess.Popen(submitStr, shell=True, stdout = subprocess.PIPE)
            # get stdout:
            self._jobID = p.communicate()[0].rstrip()
            print "jobID = " + self._jobID
            retcode = p.returncode
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode
            else:
                print >>sys.stderr, "Child returned", retcode
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

    def wait(self):
        """Waits till a process with _jobID is finished"""
        while True:
            p = subprocess.Popen('qstat -f ' + self._jobID, shell = True, \
                                  stdout = subprocess.PIPE)
            qstatOut = p.communicate()[0]
            if 'exit_status' not in qstatOut:
                time.sleep(3) # wait three seconds
                continue
            else:
                # check for exit code:
                for line in qstatOut.splitlines():
                    if 'exit_status' in line:
                        exitcode = int(line.split()[2])
                        if exitcode != 0:
                            print >>sys.stderr, 'exitcode = ', exitcode
                            print >>sys.stderr, "Quantum Espresso crashed: " +\
                     "check your settings and/or clean your 'outdir' directory"
                            raise Exception()
                break


    def serial(self, cmdStr, torqueParams = None):
        self.submit(cmdStr, torqueParams)
        self.wait()



if __name__ == "__main__":
    print 'Hello World';

__author__="Nikolay Markovskiy"
__date__ ="$Oct 6, 2009 10:13:47 PM$"