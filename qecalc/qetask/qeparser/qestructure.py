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
        return self.element + ' ' + str(self.mass) + ' ' + self.pseudopotential



class QEStructure():
    
    def __init__(self, qeConf):
        """the structure is initialized from PWSCF config file
           'lattice' and 'structure' are automatically updated"""
        self.filename = qeConf.filename
        self.atomicSpecies = OrderedDict()
        self.lattice = None
        # optConstraints three 1/0 for each coordinate of each atom
        self.optConstraints = []
        self.qeConf = qeConf
        #self.qeConf.parse()
        self.setStructureFromQEInput()


    def setStructureFromQEInput(self):
        """ Loads structure from PWSCF config file"""
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


    def updatePWInput(self, qeConf = None):
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
                coords = self.lattice.diffpy().cartesian(atom.xyz)
                coords = str(coords/self.lattice.a)[1:-1]
            else:
                if self.atomicPositionsType == 'crystal':
                    #coords = str(atom.xyz)[1:-1]
                    coords = '%f.8  %f.8  %f.8'%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                else:
                    raise NonImplementedError
            line = atom.element + '    ' + coords + '  ' + str(constraint)[1:-1]
#            line = atom.element + ' ' + coords + ' ' + str(constraint)[1:-1]
            qeConf.card('atomic_positions').addLine(line)

        # update ATOMIC_SPECIES card
        if 'atomic_species' in qeConf.cards:
            qeConf.removeCard('atomic_species')
        qeConf.createCard('atomic_species')
        for element, specie in self.atomicSpecies.items():
            qeConf.card('atomic_species').addLine(specie.toString())


    def saveStructureToPWSCF(self, fname = None):
        """Writes/updates structure into PWSCF config file,
           if the file does not exist, new one will be created"""
        if fname != None:
            filename = fname
            self.lattice.saveLatticeToPWSCF(filename)
            qeConf = QEInput(fname)
            qeConf.parse()
        else:
            filename = self.filename
            self.lattice.saveLatticeToPWSCF(filename)
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
    pwInput.parse()
    myStruct = QEStructure(qeConf = pwInput)
    myStruct.lattice.ibrav = 4
    print myStruct.lattice.a
    print myStruct.lattice.c
    myStruct.lattice.a = 43
    #pwInput.save()
    myStruct.saveStructureToPWSCF('scf_2.in')
    print myStruct.structure
