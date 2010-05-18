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

from baseoutput import BaseOutput

class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        self.parsers = {
                        'trajectory'       : self.getTrajectory,
                        }
    def getTrajectory(self, setting):

        trajectory = {
                       'pos': [],
                       'vel': [],}

        #read Espresso output into memory:
        file = open(setting.get('cpOutput'))
        cpOut = file.readlines()
        posList =  [i for i,line in enumerate(cpOut) \
                                    if '* Physical Quantities at step:' in line]
        for i, iStart in enumerate(posList):
            if i < len(posList) - 1:
                iEnd = i+1
            else:
                iEnd = len(cpOut) - 1
            for i in range(iStart, iEnd):
                line = cpOut[i]
                if 'ATOMIC_POSITIONS' in line:                
                    atoms = []
                    j = i + 1
                    while cpOut[j].strip() != "":
                        coord = [ float(w) for w in cpOut[j].split()[1:]]
                        atoms.append(coord)
                        j = j + 1
                    trajectory['pos'].append(atoms)
        #print cpOut
        #print trajectory
        return [(trajectory, None)]
        

if __name__ == "__main__":
    print "Hello World";
