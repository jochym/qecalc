#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import unittest
import os

from qecalc.qetask import *

taskDic = {'pw'       : PWTask,
          'ph'        : PHTask,
          'matdyn'    : MatdynTask,
          'projwfc'   : ProjwfcTask,
          'dynmat'    : DynmatTask,
          'dos'       : DOSTask,
          'd3'        : D3Task,
          'cp'        : CPTask,        
          }


#taskInputs = ['scf.in', 'ph.in', 'matdyn.in', 'dynmat.in', 'q2r.in', 'd3.in', \
#                                                'projwfc.in', 'dos.in', 'cp.in']
# variables
tests_dir = os.path.dirname(os.path.abspath('testTasks.py'))
testdata_dir = os.path.join(tests_dir, 'data')

useStringConfig = True

class TestStructureMethods(unittest.TestCase):
    
    def setUp(self):
        
        self.task = {}        
        filename = os.path.join(testdata_dir, 'tasks_config.ini') 
        for task in taskDic:
            self.task[task] = taskDic[task](filename = filename)
    
    def test_default_constructor(self):
        
        for task in taskDic:
            t = taskDic[task]()
            
    def test_cptask(self):
        
        cmdLine = 'mpiexec -n 2 cp.x -npool 2 -inp data/task_cp.in > cp.out< /dev/null'
        
        cp = self.task['cp']
        #print cp.input.toString()
        self.assertEqual(cmdLine, cp.cmdLine())
        
    
    def test_change_setting(self):
        answer = 'Al   0.000000 0.000000 0.000000 26.980000 Al.pz-vbc.UPF'
        pw = self.task['pw']
        pw.setting.set('pwInput', 'data/al_pw.in')        
        pw.input.parse()    
        self.assertEqual( answer, str(pw.input.structure[0]))
        answer = 'mpiexec -n 2 pw.x -npool 2 -inp data/al_pw.in > scf.out< /dev/null'
        self.assertEqual( pw.cmdLine(), answer )


    def test_configString(self):
        input = os.path.join(testdata_dir, 'al_pw.in') 
        output = os.path.join(testdata_dir, 'mgalb4_pw.out') 
        config = "[pw.x]\npwInput: %s\npwOutput: %s" %(input, output)

        pw = PWTask(configString = config)
        
        pw.output.parse()
        #print pw.output.property('total energy')
                
if __name__ == '__main__':
    unittest.main()               
