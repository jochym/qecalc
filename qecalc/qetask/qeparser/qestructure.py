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
from diffpy.Structure.structure import Structure
from diffpy.Structure.atom import Atom

from qelattice import QELattice
import numpy
from qeinput import QEInput
from orderedDict import OrderedDict


class AtomicSpecies():
    def __init__(self, element = 'H', mass = 1.0, pseudopotential = ''):
        self.element = element
        self.pseudopotential = pseudopotential
        self.mass = mass
    def toString(self):
        return '%-3s'%self.element + ' ' + '%.4f'%self.mass + ' ' + self.pseudopotential



class QEStructure():
    
    def __init__(self, qeConf):
        """the structure is initialized from PWSCF config file
           'lattice' and 'structure' are automatically updated"""
        self.filename = qeConf.filename
        self.atomicSpecies = OrderedDict()
        self.lattice = None
        self.formatString = '%# .8f %# .8f %# .8f'
        # optConstraints three 1/0 for each coordinate of each atom
        self.optConstraints = []
        self.qeConf = qeConf
        #self.qeConf.parse()
        #self.setStructureFromQEInput()
        self.lattice = None
        self.structure = None
        self.nat = None
        self.ntyp = None
        
        
    def parseInput(self):
        self.setStructureFromQEInput()
        
    def parseOutput(self, pwscfOutputFile):
        self.setStructureFromPWOutput(pwscfOutputFile)

#    def _parseAtomicPositions(self, pwscfOut):
#        for n in range(self.nat):
#            words = pwscfOut[i + n + 1].split()
#            atomSymbol = words[0]
#            coords = [float(w) for w in words[1:4]]
#            constraint = []
#            if len(words) > 4:
#                constraint = [int(c) for c in words[4:7]]
#            self.optConstraints.append(numpy.array(constraint, dtype = int))
#            print numpy.array(coords[0:3])*a_0
#            coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
#            self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))

    def setStructureFromPWOutput(self, pwscfOutputFile): 
        """
        Loads structure from PWSCF output file. If there was geometry
        optimization (relax or vc-relax), the structure will be reinitialized
        from the last step of the optimization
        """
        file = open(pwscfOutputFile)
        pwscfOut = file.readlines()
        pseudoList = []
        atomList = []
        massList = []
        self.atomicSpecies = OrderedDict()
        self.atomicPositionsType = 'alat'  
        # parse beginning:
        for i, line in enumerate(pwscfOut):
            if 'lattice parameter (a_0)' in line:
                a_0 = float(line.split()[4])
            if 'bravais-lattice index' in line:
                ibrav = int(line.split('=')[1])
            if 'number of atoms/cell' in line:
                self.nat = int(line.split('=')[1])
            if 'number of atomic types' in line:
                self.ntyp = int(line.split('=')[1])
            if 'PseudoPot.' in line:
                pseudoList.append(line.split('read from file')[1].strip())
            if 'atomic species   valence    mass     pseudopotential' in line:
                for j in range(self.ntyp):
                    atomList.append(pwscfOut[i+j+1].split()[0])
                    massList.append(float(pwscfOut[i+j+1].split()[2]))
            if 'crystal axes: (cart. coord. in units of a_0)' in line:
                latticeVectors = [[float(f)*a_0 for f in pwscfOut[i + 1].split()[3:6] ],
                                  [float(f)*a_0 for f in pwscfOut[i + 2].split()[3:6] ],
                                  [float(f)*a_0 for f in pwscfOut[i + 3].split()[3:6] ]]
                self.lattice = QELattice(ibrav = ibrav, base = latticeVectors)
            if 'site n.     atom                  positions (a_0 units)' in line:
                self.structure = Structure(lattice = self.lattice.diffpy())
                for n in range(self.nat):
                    words = pwscfOut[i + n + 1].split()
                    atomSymbol = words[1]
                    coords = [float(w) for w in words[6:9]]
                    constraint = []
                    self.optConstraints.append(numpy.array(constraint, dtype = int))
                    print numpy.array(coords[0:3])*a_0
                    coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
                    self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))

        for a, m, p in zip(atomList, massList, pseudoList):
            self.atomicSpecies[a] = AtomicSpecies(a, m, p)

        print 'Input structure from output file: ', self.toString()
        #Parse end:
        # Find all geometry optimization steps
        posList =  [i for i,line in enumerate(pwscfOut) if '!    total energy' in line]
        lastSection = pwscfOut[posList[-1]:]
        for i, line in enumerate(lastSection):
            if 'CELL_PARAMETERS (alat)' in line:
                latticeVectors = [[float(f)*a_0 for f in lastSection[i + 1].split() ],
                                  [float(f)*a_0 for f in lastSection[i + 2].split() ],
                                  [float(f)*a_0 for f in lastSection[i + 3].split() ]]
                self.lattice = QELattice(ibrav = 0, base = latticeVectors)
                print self.lattice.diffpy().base
            if 'ATOMIC_POSITIONS (alat)' in line:
                self.structure = Structure(lattice = self.lattice.diffpy())
                for n in range(self.nat):
                    words = lastSection[i + n + 1].split()
                    atomSymbol = words[0]
                    coords = [float(w) for w in words[1:4]]
                    constraint = []
                    if len(words) > 4:
                        constraint = [int(c) for c in words[4:7]]
                    self.optConstraints.append(numpy.array(constraint, dtype = int))
                    print numpy.array(coords[0:3])*a_0
                    coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
                    self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))
        self.lattice.ibrav = ibrav
        print 'Output structure from output file: ', self.toString()
    
    
    def setStructureFromQEInput(self):
        """ Loads structure from PWSCF config file"""
        self.atomicSpecies = OrderedDict()
        self.lattice = QELattice(qeConf = self.qeConf)
        self.structure = Structure(lattice = self.lattice.diffpy())
        self.nat  = int(self.qeConf.namelist('system').param('nat'))
        self.ntyp  = int(self.qeConf.namelist('system').param('ntyp'))
        atomicLines = self.qeConf.card('atomic_positions').lines()
        self.atomicPositionsType = self.qeConf.card('atomic_positions').arg()
        if self.atomicPositionsType == 'bohr' or self.atomicPositionsType == 'angstrom':
            raise NotImplementedError
        if self.atomicPositionsType == None:
            self.atomicPositionsType = 'alat'
        for line in atomicLines:
            if '!' not in line:
                words = line.split()
                coords = [float(w) for w in words[1:4]]
                constraint = []
                if len(words) > 4:
                    constraint = [int(c) for c in words[4:7]]
                self.optConstraints.append(numpy.array(constraint, dtype = int))
                atomSymbol = words[0]
                if self.atomicPositionsType == 'alat':
                    coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*self.lattice.a0)
                if self.atomicPositionsType == 'crystal':
                    coords = numpy.array(coords[0:3])
                self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))
                # parse mass ATOMIC_SPECIES section:
                atomicSpeciesLines = self.qeConf.card('atomic_species').lines()
                for line in atomicSpeciesLines:
                    if '!' not in line:
                        atomicSpeciesWords = line.split()
                        element = atomicSpeciesWords[0]
                        mass = float(atomicSpeciesWords[1])
                        ps = atomicSpeciesWords[2]
                        self.atomicSpecies[element] =  AtomicSpecies(element, mass, ps)


    def toString(self):
        s = self.lattice.toString() + '\n'
        if self.atomicPositionsType == 'alat':
            s = s + 'Atomic positions in units of lattice parametr "a":\n'        
        if self.atomicPositionsType == 'crystal':
            s = s + 'Atomic positions in crystal coordinates:\n'
        for atom, constraint in zip(self.structure, self.optConstraints):
            if self.atomicPositionsType == 'alat':      
                coords = self.lattice.diffpy().cartesian(atom.xyz)/self.lattice.a
                coords = self.formatString%(coords[0], coords[1], coords[2])
                #coords = str(coords/self.lattice.a)[1:-1]
            else:
                if self.atomicPositionsType == 'crystal':
                    #coords = str(atom.xyz)[1:-1]
                    coords = self.formatString%(v[0], v[1], v[2])%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                else:
                    raise NonImplementedError
            s = s + '%-3s'%atom.element + '    ' + coords + '  ' \
                    + str(constraint)[1:-1] + '\n'

        s = s + '\n'
        for element, specie in self.atomicSpecies.items():
            s = s + specie.toString() + '\n'

        return s


    def updatePWInput(self, qeConf = None):

        self.lattice.updatePWInput()

        qeConf.namelist('system').remove('ntyp')
        qeConf.namelist('system').remove('nat')
        qeConf.namelist('system').add('ntyp', self.ntyp)
        qeConf.namelist('system').add('nat', self.nat)

        if 'atomic_positions' in qeConf.cards:
            qeConf.removeCard('atomic_positions')
        qeConf.createCard('atomic_positions')
        qeConf.card('atomic_positions').setArg(self.atomicPositionsType)
        for atom, constraint in zip(self.structure, self.optConstraints):
            if self.atomicPositionsType == 'alat':
                coords = self.lattice.diffpy().cartesian(atom.xyz)/self.lattice.a
                coords = self.formatString%(coords[0], coords[1], coords[2])
            else:
                if self.atomicPositionsType == 'crystal':
                    #coords = str(atom.xyz)[1:-1]
                    coords = self.formatString%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                else:
                    raise NonImplementedError
            line = '%-3s'%atom.element + '    ' + coords + '  ' + str(constraint)[1:-1]
#            line = atom.element + ' ' + coords + ' ' + str(constraint)[1:-1]
            qeConf.card('atomic_positions').addLine(line)

        # update ATOMIC_SPECIES card
        if 'atomic_species' in qeConf.cards:
            qeConf.removeCard('atomic_species')
        qeConf.createCard('atomic_species')
        for element, specie in self.atomicSpecies.items():
            qeConf.card('atomic_species').addLine(specie.toString())


    def save(self, fname = None):
        """Writes/updates structure into PW config file,
           if the file does not exist, new one will be created"""
        if fname != None:
            filename = fname
            self.lattice.save(filename)
            qeConf = QEInput(fname)
            qeConf.parse()
        else:
            filename = self.filename
            self.lattice.save(filename)
            qeConf = self.qeConf
            #qeConf.parse()
        self.updatePWInput(qeConf )
            
        qeConf.save(filename)

    def diffpy(self):
        return self.structure

 #       def placeInLattice(self, new_lattice):

#        def getLattice(self):
 #           return self.structure.lattice

if __name__ == '__main__':
    pwInput = QEInput('scf.in', type = 'pw')
    #pwInput.parse()
    myStruct = QEStructure(qeConf = pwInput)
    #myStruct.lattice.ibrav = 4
    #print myStruct.lattice.a
    #print myStruct.lattice.c
    #myStruct.lattice.a = 43
    #pwInput.save()
    #myStruct.saveStructureToPWSCF('scf_2.in')
    #print myStruct.structure
    myStruct.parseOutput('geom.out')
    print myStruct.toString()