from diffpy.Structure.structure import Structure
from diffpy.Structure.atom import Atom
from qelattice import QELattice
import numpy
from parser.configParser import *


class AtomicSpecies():
    def __init__(self, element = 'H', mass = 1.0, pseudopotential = ''):
        self.element = element
        self.pseudopotential = pseudopotential
        self.mass = mass
    def toString(self):
        return self.element + ' ' + str(self.mass) + ' ' + self.pseudopotential



class QEStructure():

    
    def __init__(self, fname):
        """the structure is initialized from PWSCF config file
           'lattice' and 'structure' are automatically updated"""
        self.filename = fname
        self.atomicSpecies = {}
        self.lattice = None
        # optConstraints three 1/0 for each coordinate of each atom
        self.optConstraints = []
        self.qeConf = QEConfig(fname)
        self.qeConf.parse()
        self.setStructureFromPWSCF()

    
    def setStructureFromPWSCF(self):
        """ Loads structure from PWSCF config file"""
        self.lattice = QELattice(fname = self.filename)
        self.structure = Structure(lattice = self.lattice.diffpy())
        self.nat  = int(self.qeConf.namelist('system').param('nat'))
        self.ntyp  = int(self.qeConf.namelist('system').param('ntyp'))
        atomicLines = self.qeConf.card('atomic_positions').getLines()
        self.atomicPositionsType = self.qeConf.card('atomic_positions').argument()
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
                atomicSpeciesLines = self.qeConf.card('atomic_species').getLines()
                for line in atomicSpeciesLines:
                    if '!' not in line:
                        atomicSpeciesWords = line.split()
                        element = atomicSpeciesWords[0]
                        mass = float(atomicSpeciesWords[1])
                        ps = atomicSpeciesWords[2]
                        self.atomicSpecies[element] =  AtomicSpecies(element, mass, ps)


    def saveStructureToPWSCF(self, fname = None):
        """Writes/updates structure into PWSCF config file,
           if the file does not exist, new one will be created"""
        if fname != None:
            filename = fname
            self.lattice.saveLatticeToPWSCF(filename)
            qeConf = QEConfig(fname)
            qeConf.parse()
        else:
            filename = self.filename
            self.lattice.saveLatticeToPWSCF(filename)
            qeConf = self.qeConf
            qeConf.parse()

        qeConf.namelist('system').removeParam('ntyp')
        qeConf.namelist('system').removeParam('nat')
        qeConf.namelist('system').addParam('ntyp', self.ntyp)
        qeConf.namelist('system').addParam('nat', self.nat)

        if 'atomic_positions' in qeConf.cards:
            qeConf.removeCard('atomic_positions')
        qeConf.createCard('atomic_positions')
        qeConf.card('atomic_positions').setArgument(self.atomicPositionsType)
        for atom, constraint in zip(self.structure, self.optConstraints):
            if self.atomicPositionsType == 'alat':
                coords = self.lattice.diffpy().cartesian(atom.xyz)
                coords = str(coords/self.lattice.a)[1:-1]
            else:
                if self.atomicPositionsType == 'crystal':
                    coords = str(atom.xyz)[1:-1]
                else:
                    raise NonImplementedError
            line = atom.element + ' ' + coords + ' ' + str(constraint)[1:-1]
            qeConf.card('atomic_positions').addLine(line)

        # update ATOMIC_SPECIES card
        if 'atomic_species' in qeConf.cards:
            qeConf.removeCard('atomic_species')
        qeConf.createCard('atomic_species')
        for element, specie in self.atomicSpecies.items():
            qeConf.card('atomic_species').addLine(specie.toString())
            
        qeConf.save(filename)

    def diffpy(self):
        return self.structure

 #       def placeInLattice(self, new_lattice):

#        def getLattice(self):
 #           return self.structure.lattice

if __name__ == '__main__':
    myStruct = QEStructure('scf2.in')
    myStruct.lattice.ibrav = 4
    print myStruct.lattice.a
    print myStruct.lattice.c
    myStruct.saveStructureToPWSCF('qwe.in')
    print myStruct.structure
