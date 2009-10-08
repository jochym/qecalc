# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="markovsk"
__date__ ="$Oct 6, 2009 10:13:47 PM$"
import os
import sys
import launcher
class QETorque():    
    def __init__(self):
        self._workDir = os.getcwd()
        self._jobID = None
        # use qmgr 
    def submit(self, cmdStr):
        """Submits job. cmdStr is mpirun + params + program is if one did not
           use torque"""
        scriptName = 'submit'
	file = open(scriptName, 'w')
	file.write('#!/bin/sh\n' + cmdStr)
	file.close()		
        submitStr = "echo '" + cmdStr + "' | " + 'qsub -V -d ' + \
                     self._workDir + ' ' + torqueResourceList + ' -'
        print "submitStr = " + submitStr
	os.system(submitStr)
        
    def serial(self, cmdStr):
        if self._jobID == None:
            self.submit(cmdStr)


if __name__ == "__main__":
    print "Hello World";
