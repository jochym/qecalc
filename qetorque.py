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
import subprocess
import sys
#import launcher
import os
import time
class QETorque:
    def __init__(self, fname):
        self._workDir = os.getcwd()
        self._jobID = None

        configDic = {
        'torqueResourceList': '-l nodes=1:ppn=1'
        }
        self.config = ConfigParser.SafeConfigParser(configDic)
        self.config.read(fname)
        self.torqueResourceList = self.config.get('Setting', 'torqueResourceList')

        # use qmgr (not implemented)
    def submit(self, cmdStr):
        """Submits job. cmdStr is mpirun + params + program is if one did not
           use torque"""
        submitStr = "echo '" + cmdStr + "' | " + 'qsub -V -d ' + \
                     self._workDir + ' ' + self.torqueResourceList + ' -'
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
                for line in qstatOut:
                    if 'exit_status' in line:
                        exitcode = line.split()[2]
                        if exitcode != 0:
                            raise Exception("Quantum Espresso crashed: " +\
                     "check your settings and/or clean your 'outdir' directory")
                break
        print "qstatOut = " + qstatOut


    def serial(self, cmdStr):
        self.submit(cmdStr)
        self.wait()



if __name__ == "__main__":
    torque = QETorque()
    torque.serial('ls')

__author__="Nikolay Markovskiy"
__date__ ="$Oct 6, 2009 10:13:47 PM$"