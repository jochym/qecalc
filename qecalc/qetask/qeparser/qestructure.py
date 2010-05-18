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

try:
    from matter import Structure, Atom, Lattice
    from matter.Lattice import cosd
    from matter.SymmetryUtilities import equalPositions
except ImportError:
    from diffpy.Structure.structure import Structure
    from diffpy.Structure.atom import Atom
    from diffpy.Structure.lattice import cosd, Lattice
    from diffpy.Structure.SymmetryUtilities import equalPositions


#from matter import Structure, Atom, Lattice
#from matter.Lattice import cosd
#from matter.SymmetryUtilities import equalPositions

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
        self.formatString = '%# .8f %# .8f %# .8f'
        # optConstraints three 1/0 for each coordinate of each atom
        self.optConstraints = []
        self.qeConf = qeConf
        self.lattice = QELattice()
        self.structure = Structure(lattice = self.lattice.diffpy())
        self.nat = None
        self.ntyp = None
        self.atomicPositionsType = 'crystal'       
        

    def atomLabels(self):
        labels = []
        for l in self.atomicSpecies:
            labels.append(l)
        return labels
        
    def parseInput(self):
        self.setStructureFromQEInput()
        
    def parseOutput(self, pwscfOutputFile):
        self.setStructureFromPWOutput(pwscfOutputFile)

    def setStructureFromPWOutput(self, pwscfOutputFile): 
        """
        Loads structure from PWSCF output file. If there was geometry
        optimization (relax or vc-relax), the structure will be reinitialized
        from the last step of the optimization. Assumes output is in ALAT units
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
                self.lattice.setLatticeFromQEVectors(ibrav, latticeVectors)
            if 'site n.     atom                  positions (a_0 units)' in line:
                self.structure = Structure(lattice = self.lattice.diffpy())
                for n in range(self.nat):
                    words = pwscfOut[i + n + 1].split()
                    atomSymbol = words[1]
                    coords = [float(w) for w in words[6:9]]
                    constraint = []
                    self.optConstraints.append(numpy.array(constraint, dtype = int))
                    #print numpy.array(coords[0:3])*a_0
                    coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
                    self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))

        for a, m, p in zip(atomList, massList, pseudoList):
            self.atomicSpecies[a] = AtomicSpecies(a, m, p)

        #print 'Input structure from output file: ', self.toString()
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
                return


        for i, line in enumerate(lastSection):
            if 'CELL_PARAMETERS (alat)' in line:
                latticeVectors = [[float(f)*a_0 for f in lastSection[i + 1].split() ],
                                  [float(f)*a_0 for f in lastSection[i + 2].split() ],
                                  [float(f)*a_0 for f in lastSection[i + 3].split() ]]
                self.lattice.setLatticeFromQEVectors(ibrav, latticeVectors)
                #self.lattice = QELattice(ibrav = ibrav, base = latticeVectors)
                #print self.lattice.diffpy().base
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
                    coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*a_0)
                    self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))
    
    
    def setStructureFromQEInput(self):
        """ Loads structure from PWSCF config file"""
        self.atomicSpecies = OrderedDict()
        self.lattice.setLatticeFromPWInput(self.qeConf)
        #self.lattice = QELattice(qeConf = self.qeConf)
        self.structure = Structure(lattice = self.lattice.diffpy())
        self.nat = self.ntyp = None
        self.filename = self.qeConf.filename
        self.optConstraints = []
        
        if 'system' in self.qeConf.namelists:        
            self.nat  = int(self.qeConf.namelist('system').param('nat'))
            self.ntyp  = int(self.qeConf.namelist('system').param('ntyp'))
        if 'atomic_positions' in self.qeConf.cards:        
            atomicLines = self.qeConf.card('atomic_positions').lines()
            self.atomicPositionsType = self.qeConf.card('atomic_positions').arg()
            if self.atomicPositionsType == 'angstrom':
                raise NotImplementedError\
         ('atomic positions in bohr and angstrom are not currently supported')
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
                        coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3])*self.lattice.a)
                    if self.atomicPositionsType == 'crystal':
                        coords = numpy.array(coords[0:3])
                    if self.atomicPositionsType == 'bohr':
                        coords = self.lattice.diffpy().fractional(numpy.array(coords[0:3]))
                    self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))
        # parse mass ATOMIC_SPECIES section:
         
        if 'atomic_species' in self.qeConf.cards:
            atomicSpeciesLines = self.qeConf.card('atomic_species').lines()
            for line in atomicSpeciesLines:
                if '!' not in line:
                    if line.strip() != '':                     
                        atomicSpeciesWords = line.split()
                        element = atomicSpeciesWords[0]
                        mass = 0
                        ps = ''
                        if len(atomicSpeciesWords) > 1 :
                            mass = float(atomicSpeciesWords[1])
                        if len(atomicSpeciesWords) > 2:
                            ps = atomicSpeciesWords[2]
                        self.atomicSpecies[element] =  AtomicSpecies(element, mass, ps)


    def load(self, source, **args):
        task = {
            'diffpy': self.setStructureFromDiffpyStructure
        }
        if source == 'diffpy':
            if 'ibrav' in args and args['ibrav'] != 0:
                task['diffpy'] = self.setReducedStructureFromDiffpyStructure

        task[source](**args)
        
        self.updatePWInput(qeConf = self.qeConf)


    def setStructureFromDiffpyStructure(self, structure, massList = [], psList = [], ibrav = 0):
        """
        structure - diffpy.Structure object
        ibrav - Lattice index
        psList - list of strings with pseudopotential names
        diffpyStructure object will be modified with reduced atomic positions
        """      
        diffpyLattice = structure.lattice
        
        self.structure = structure
               
        
        #set lattice and  convert to bohr units
        #qeLattice = QELattice(ibrav = 0, a = 1.889725989, base = diffpyLattice.base)
        qeLattice = QELattice(ibrav = 0, base = diffpyLattice.base)
        qeLattice.a = 1.889725989*qeLattice.a
        
        
        self.lattice = qeLattice
        self.lattice.type = 'generic cubic'

        atomNames = []
        for a in structure:
            if self._element(a) not in atomNames:
                atomNames.append(self._element(a))
        
        #print atomNames
        #print len(massList)
        for i, elem in enumerate(atomNames):
            if len(massList) - 1 < i:
                mass = 0
            else:
                mass = massList[i]
            if len(psList) - 1 < i:
                ps = ''
            else:
                ps = psList[i]               
            self.atomicSpecies[elem] =  AtomicSpecies(elem, mass, ps)
        
        for atom in structure:
            self.optConstraints.append([])        
        
        # for i, atom in enumerate(structure):
            # elem = self._element(atom)
            # if len(massList) - 1 < i:
                # mass = 0
            # else:
                # mass = massList[i]
            # if len(psList) - 1 < i:
                # ps = ''
            # else:
                # ps = psList[i]
            # self.atomicSpecies[elem] =  AtomicSpecies(elem, mass, ps)
            # self.optConstraints.append([])

#        for atom, mass, ps in zip(structure, massList, psList):
#            elem = self._element(atom)
#            self.atomicSpecies[elem] =  AtomicSpecies(elem, mass, ps)
#            self.optConstraints.append([])

        self.nat = len(structure)
        self.ntyp = len(self.atomicSpecies)        
     
                        
    def setReducedStructureFromDiffpyStructure(self, structure, ibrav, massList = [], psList = []):
        """
        structure - diffpy.Structure object
        ibrav - Lattice index
        psList - list of strings with pseudopotential names
        diffpyStructure object will be modified with reduced atomic positions
        """

        #self.atomicSpecies = OrderedDict()
        #self.optConstraints = []
        #self.atomicPositionsType = 'crystal'

        diffpyLattice = structure.lattice

        a = diffpyLattice.a
        b = diffpyLattice.b
        c = diffpyLattice.c
        cAB = cosd(diffpyLattice.gamma)
        cBC = cosd(diffpyLattice.alpha)
        cAC = cosd(diffpyLattice.beta)

        qeLattice = QELattice(ibrav = ibrav, a = a, b = b, c = c,  cBC =  cBC, \
                              cAC = cAC, cAB = cAB)

        self.lattice = qeLattice
        # make a deep copy (does not wok now)
        #reducedStructure = Structure(diffpyStructure)
        reducedStructure = structure

        reducedStructure.placeInLattice(Lattice(base=qeLattice.diffpy().base))

        # collect atoms that are at equivalent position to some previous atom
        duplicates = set([a1
            for i0, a0 in enumerate(reducedStructure) for a1 in reducedStructure[i0+1:]
                if   self._element(a0) == self._element(a1) and equalPositions(a0.xyz, a1.xyz, eps=1e-4)])

        
        # Filter out duplicate atoms.  Use slice assignment so that
        # reducedStructure is not replaced with a list.
        reducedStructure[:] = [a for a in reducedStructure if not a in duplicates]

        self.structure = reducedStructure

        atomNames = []
        for a in reducedStructure:
            if self._element(a) not in atomNames:
                atomNames.append(self._element(a))
        
        #print atomNames
        #print len(massList)
        for i, elem in enumerate(atomNames):
            if len(massList) - 1 < i:
                mass = 0
            else:
                mass = massList[i]
            if len(psList) - 1 < i:
                ps = ''
            else:
                ps = psList[i]      
            #print mass, ps
            # atomDict[a] = 
        # for i, atom in enumerate(reducedStructure):
            # elem = self._element(atom)
            # if len(massList) - 1 < i:
                # mass = 0
            # else:
                # mass = massList[i]
            # if len(psList) - 1 < i:
                # ps = ''
            # else:
                # ps = psList[i]            
            self.atomicSpecies[elem] =  AtomicSpecies(elem, mass, ps)
        
        for atom in reducedStructure:
            self.optConstraints.append([])

        # convert to bohr units
        self.lattice.setLattice(ibrav, self.lattice.a*1.889725989, \
                                 self.lattice.b*1.889725989,
                                 self.lattice.c*1.889725989)

        self.nat = len(reducedStructure)
        self.ntyp = len(self.atomicSpecies)

        # use rstrip to avoid duplicate line feed
        #print reducedStructure.writeStr(format='discus').rstrip()
        #print reducedStructure.writeStr(format='discus')
        #print reducedStructure.writeStr().rstrip()
        #print reducedStructure
        #self.lattice = setLatticeFromDiffpyLattice(structure.lattice, ibrav)


    # def toString(self):
        # s = self.lattice.toString() + '\n'
        # if self.atomicPositionsType == 'alat':
            # s = s + 'Atomic positions in units of lattice parametr "a":\n'        
        # if self.atomicPositionsType == 'crystal':
            # s = s + 'Atomic positions in crystal coordinates:\n'
        # for atom, constraint in zip(self.structure, self.optConstraints):
            # if self.atomicPositionsType == 'alat':      
                # coords = self.lattice.diffpy().cartesian(atom.xyz)/self.lattice.a
                # coords = self.formatString%(coords[0], coords[1], coords[2])
            # else:
                # if self.atomicPositionsType == 'crystal':
                    # coords = self.formatString%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                # else:
                    # raise NonImplementedError
            # s = s + '%-3s'%self._element(atom) + '    ' + coords + '  ' \
                    # + str(constraint)[1:-1] + '\n'

        # s = s + '\n'
        # for element, specie in self.atomicSpecies.items():
            # s = s + specie.toString() + '\n'

        # return s

    def toString(self, string = None):
        if string != None:
            string = self.lattice.toString(string = string)
            qeConf = QEInput(config = string)
            qeConf.parse()
        else:
            if self.qeConf != None:
                qeConf = self.qeConf
            else:
                qeConf = QEInput(config = '')

        self.updatePWInput(qeConf)
        return qeConf.toString()

        
    def updatePWInput(self, qeConf = None):

        if qeConf == None:
            qeConf = self.qeConf

        self.lattice.updatePWInput(qeConf)

        if 'system' not in qeConf.namelists:
            qeConf.addNamelist('system')            
        qeConf.namelist('system').remove('ntyp')
        qeConf.namelist('system').remove('nat')
        if self.ntyp != None:
            qeConf.namelist('system').add('ntyp', self.ntyp)
        if self.nat != None:
            qeConf.namelist('system').add('nat', self.nat)
        
        if len(qeConf.namelist('system').params) == 0:
            qeConf.removeNamelist('system')  

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
                    coords = self.formatString%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                else:
                    if self.atomicPositionsType == 'bohr':
                        coords = self.lattice.diffpy().cartesian(atom.xyz)
                        coords = self.formatString%(coords[0], coords[1], coords[2])
                    else:
                        raise NonImplementedError
            line = '%-3s'%self._element(atom) + '    ' + coords + '  ' + str(constraint)[1:-1]
            qeConf.card('atomic_positions').addLine(line)
        
        if len(qeConf.card('atomic_positions').lines()) == 0:
            qeConf.removeCard('atomic_positions')

        # update ATOMIC_SPECIES card
        if 'atomic_species' in qeConf.cards:
            qeConf.removeCard('atomic_species')
        qeConf.createCard('atomic_species')
        for element, specie in self.atomicSpecies.items():
            qeConf.card('atomic_species').addLine(specie.toString())
        
        if len(qeConf.card('atomic_species').lines()) == 0:
            qeConf.removeCard('atomic_species')        


    def save(self, fname = None):
        """Writes/updates structure into PW config file,
           if the file does not exist, new one will be created"""
        filename = fname
        if fname != None:
            self.lattice.save(filename)
            qeConf = QEInput(fname)
            qeConf.parse()
        else:
            filename = self.filename
            self.lattice.save(filename)
            qeConf = self.qeConf
            
        self.updatePWInput(qeConf)
        qeConf.save(filename)

        

    def diffpy(self):
        return self.structure


    def _element(self, atom):
        """
        Is needed for suport both diffpy and matter classess
        """
        if 'element' in dir(atom):
            return atom.element
        else:
            if 'symbol' in dir(atom):
                return atom.symbol
            else:
                raise

if __name__ == '__main__':
    print "Hello";