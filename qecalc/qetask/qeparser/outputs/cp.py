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
                       'pos':    [],
                       'vel':    [],
                       'forces': [],
                       'etot':   [],
                     }

        #read Espresso output into memory:
        file = open(setting.get('cpOutput'))
        cpOut = file.readlines()
        posList =  [i for i,line in enumerate(cpOut) \
                                    if '* Physical Quantities at step:' in line]
        for i, iStart in enumerate(posList):
            if i < len(posList) - 1:
                iEnd = posList[i+1]
            else:
                iEnd = len(cpOut) - 1
            for iLine in range(iStart, iEnd):
                line = cpOut[iLine]
                if 'ATOMIC_POSITIONS' in line:
                    self._appendAtomicProperty(iLine, cpOut, 'pos', trajectory)
                if 'ATOMIC_VELOCITIES' in line:
                    self._appendAtomicProperty(iLine, cpOut, 'vel', trajectory)
                if 'Forces acting on atoms (au):' in line:
                    self._appendAtomicProperty(iLine, cpOut, 'forces', trajectory)
                if 'total energy =' in line:
                    trajectory['etot'].append(float(line.split()[3]))

        return [(trajectory, None)]

    def _appendAtomicProperty(self, iPos, cpOut, propName, trajectory):
        """
        iPos - position number where the atomic property of interest starts
               (E.g. atomic_positions)
        propName = desired name in trajectory dictionary
        """
        atoms = []
        j = iPos + 1
        while cpOut[j].strip() != "":
            coord = [ float(w) for w in cpOut[j].split()[1:]]
            atoms.append(coord)
            j = j + 1
        trajectory[propName].append(atoms)

if __name__ == "__main__":
    print "Hello World";
