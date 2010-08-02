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

import StringIO
from qecalc.qetask.qeparser.qestructureparser import *

class P_pwoutput(QEStructureParser):
    """
        Loads structure from pw.x output file. If there was geometry
        optimization (relax or vc-relax), the structure will be reinitialized
        from the last step of the optimization. Assumes output is in ALAT units
    """    
    
    def __init__(self, qeConf):
        QEStructureParser.__init__(self, qeConf)
        self.format = "pwoutput"
    
    
    def parseStr(self, s):
        """Create Structure instance from a string."""
        file = StringIO.StringIO(s)
        return self.__genStructure(file)
         

    def parse(self, filename):
        """Create Structure instance from a file."""        
        file = open(filename, 'r')
        self.filename = filename
        return self.__genStructure(file)
            
            
    def __genStructure( self, file ):            
        # silly copy constructor:
        #from qecalc.qetask.qeparser.pwinput import PWInput
        #new_qeConf = PWInput(config = self.qeConf.toString()) 
        #new_qeConf.parse()
        #*************
        stru = QEStructure(qeConf = self.qeConf)     
        pwscfOut = file.readlines()
        pseudoList = []
        atomList = []
        massList = []
        stru.atomicSpecies = OrderedDict()
        stru.atomicPositionsType = 'alat'  
        # parse beginning:
        for i, line in enumerate(pwscfOut):
            if 'lattice parameter (a_0)' in line:
                a_0 = float(line.split()[4])
            if 'bravais-lattice index' in line:
                ibrav = int(line.split('=')[1])
            if 'number of atoms/cell' in line:
                stru.nat = int(line.split('=')[1])
            if 'number of atomic types' in line:
                stru.ntyp = int(line.split('=')[1])
            if 'PseudoPot.' in line:
                pseudoList.append(line.split('read from file')[1].strip())
            if 'atomic species   valence    mass     pseudopotential' in line:
                for j in range(stru.ntyp):
                    atomList.append(pwscfOut[i+j+1].split()[0])
                    massList.append(float(pwscfOut[i+j+1].split()[2]))
            if 'crystal axes: (cart. coord. in units of a_0)' in line:
                latticeVectors = [[float(f)*a_0 for f in pwscfOut[i + 1].split()[3:6] ],
                                  [float(f)*a_0 for f in pwscfOut[i + 2].split()[3:6] ],
                                  [float(f)*a_0 for f in pwscfOut[i + 3].split()[3:6] ]]
                stru.lattice.setLatticeFromQEVectors(ibrav, latticeVectors)
            if 'site n.     atom                  positions (a_0 units)' in line:
                stru.structure = Structure(lattice = stru.lattice.diffpy())
                for n in range(stru.nat):
                    words = pwscfOut[i + n + 1].split()
                    atomSymbol = words[1]
                    coords = [float(w) for w in words[6:9]]
                    constraint = []
                    stru.optConstraints.append(numpy.array(constraint, dtype = int))
                    coords = stru.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
                    stru.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))

        for a, m, p in zip(atomList, massList, pseudoList):
            stru.atomicSpecies[a] = AtomicSpecies(a, m, p)

        #Parse end:
        # Find all geometry optimization steps

        posList =  [i for i,line in enumerate(pwscfOut) if '!    total energy' in line]
        lastSection = pwscfOut[posList[-1]:]        
        # check if geometry is there:
        cellCheck = False
        atomCheck = False
        for i, line in enumerate(lastSection):
            if 'CELL_PARAMETERS (alat)' in line:
                cellCheck = True
            if 'ATOMIC_POSITIONS (alat)' in line:
                atomCheck = True

        if not cellCheck or not atomCheck:
            if len(posList) > 1:
                lastSection = pwscfOut[posList[-2]:]
            else:
                return stru
            
        for i, line in enumerate(lastSection):
            if 'CELL_PARAMETERS (alat)' in line:
                latticeVectors = [[float(f)*a_0 for f in lastSection[i + 1].split() ],
                                  [float(f)*a_0 for f in lastSection[i + 2].split() ],
                                  [float(f)*a_0 for f in lastSection[i + 3].split() ]]
                stru.lattice.setLatticeFromQEVectors(ibrav, latticeVectors)
            if 'ATOMIC_POSITIONS (alat)' in line:
                stru.structure = Structure(lattice = stru.lattice.diffpy())
                for n in range(stru.nat):
                    words = lastSection[i + n + 1].split()
                    atomSymbol = words[0]
                    coords = [float(w) for w in words[1:4]]
                    constraint = []
                    if len(words) > 4:
                        constraint = [int(c) for c in words[4:7]]
                    stru.optConstraints.append(numpy.array(constraint, dtype = int))
                    coords = stru.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
                    stru.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))        
          
        return stru
    
def getParser(qeConf):
    return P_pwoutput(qeConf)    