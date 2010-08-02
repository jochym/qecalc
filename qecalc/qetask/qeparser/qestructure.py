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
    from diffpy.Structure.structure import Structure
    from diffpy.Structure.atom import Atom
    from diffpy.Structure.lattice import cosd, Lattice
    from diffpy.Structure.SymmetryUtilities import equalPositions
except ImportError:
    from matter import Structure, Atom, Lattice
    from matter.Lattice import cosd
    from matter.SymmetryUtilities import equalPositions

from qelattice import QELattice
import numpy
from qeinput import QEInput
from orderedDict import OrderedDict


class AtomicSpecies():
    def __init__(self, element = 'H', mass = 1.0, pseudopotential = ''):
        self._element = element
        self.pseudopotential = pseudopotential
        self.mass = mass
    def __str__(self):
        return '%-3s'%self.element + ' ' + '%.4f'%self.mass + ' ' + self.pseudopotential
    def toString(self):
        return str(self)
        #return '%-3s'%self.element + ' ' + '%.4f'%self.mass + ' ' + self.pseudopotential
    def _get_element(self):
        return self._element

    def _set_element(self, value):
        self._element = value

    element = property(_get_element, _set_element, doc ="element")        


class QEStructure( list ):
    
    def __init__(self, qeInput = None):
        """the structure is initialized from PWSCF config file
           'lattice' and 'structure' are automatically updated"""
        #self.filename = qeInput.filename
        self.atomicSpecies = OrderedDict()
        self.formatString = '%# .8f %# .8f %# .8f'
        # optConstraints three 1/0 for each coordinate of each atom
        self._optConstraints = []
        self.lattice = QELattice()
        self.lattice._qeInput = qeInput
        self._qeInput = qeInput
        if qeInput != None:
            self.lattice._qeInput.structure = self
            self.lattice._qeInput.structure.lattice = self.lattice
            
        self.structure = Structure(lattice = self.lattice.diffpy())
        self._nat = None
        self._ntyp = None
        self._atomicPositionsType = 'crystal'                           
         
        
    def _get_nat(self):
        return self._nat

    def _set_nat(self, value):
        self._nat = value
        self.lattice._qeInput.update()

    nat = property(_get_nat, _set_nat, doc ="number of atoms")
    
 
    def _get_ntyp(self):
        return self._ntyp

    def _set_ntyp(self, value):
        self._ntyp = value
        self.lattice._qeInput.update()

    ntyp = property(_get_ntyp, _set_ntyp, doc ="number of types")
    
    
    def _get_atomicPositionsType(self):
        return self._atomicPositionsType

    def _set_atomicPositionsType(self, value):
        self._atomicPositionsType = value
        self.lattice._qeInput.update()

    atomicPositionsType = property(_get_atomicPositionsType, \
                                   _set_atomicPositionsType, \
                           doc ="type of atomic positions (crystal, alat ...)")    


    def _get_optConstraints(self):
        return self._optConstraints

    def _set_optConstraints(self, value):
        self._optConstraints = value
        self.lattice._qeInput.update()

    optConstraints = property(_get_optConstraints, _set_optConstraints, \
                               doc ="optimization constraints list")
                

    def __str__(self):
        """simple string representation"""        
        s = str(self.lattice) + '\n'
        if self.atomicPositionsType == 'alat':
            s = s + 'Atomic positions in units of lattice parametr "a":\n'        
        if self.atomicPositionsType == 'crystal':
            s = s + 'Atomic positions in crystal coordinates:\n'
        for atom, constraint in zip(self.structure, self.optConstraints):
            if self.atomicPositionsType == 'alat':
                coords = self.lattice.diffpy().cartesian(atom.xyz)/self.lattice.a
                coords = self.formatString%(coords[0], coords[1], coords[2])
            else:
                if self.atomicPositionsType == 'crystal':
                    coords = self.formatString%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                else:
                    raise NonImplementedError
            s = s + '%-3s'%self._element(atom) + '    ' + coords + '  ' \
                    + str(constraint)[1:-1] + '\n'
        s = s + '\n'
        
        for element, specie in self.atomicSpecies.items():
            s = s + specie.toString() + '\n'
            
        return s


    def atomLabels(self):
        labels = []
        for l in self.atomicSpecies:
            labels.append(l)
        return labels
        
    def parseInput(self, qeInput):
        #self._qeInput = qeInput
        from qestructureparser.qestructureparser import QEStructureParser
        new_structure = QEStructureParser(qeInput).parseqeInput()
        self.__Init(new_structure)
#        print self._qeInput.toString()
        #self.lattice._qeInput.update()
        #self._setStructureFromQEInput(qeInput)
    
    
    def __Init(self, structure):
        QEStructure.__init__(self)
        if structure is not None:
            self.__dict__.update(structure.__dict__)
            self.lattice.__dict__.update(structure.lattice.__dict__)
            self = structure        
    
    
#    def parseOutput(self, pwscfOutputFile):
#        self._setStructureFromPWOutput(pwscfOutputFile)


    def read(self, filename, format = 'pwinput'):
        """Load structure from a file, any original data become lost.

        filename -- file to be loaded
        format   -- structure formats
                    'pwinput'  - pw.x input
                    'pwoutput' - pw.x output

        Return instance of data Parser used to process file.  This
        can be inspected for information related to particular format.
        """        
        from  qecalc.qetask.qeparser.qestructureparser import parser_index
        
        if self._qeInput == None:                
            self._qeInput = PWInput()
            self._qeInput.parse()
        
        if format in parser_index:             
            module = __import__("qestructureparser.P_" + format, globals(), \
                                locals(), ['P_' + format], -1)
            parser = module.getParser(self._qeInput)
            new_structure = parser.parse(filename)
        else:            
            diffpyStruct = Structure()
            parser = diffpyStruct.read(filename, format = format)
            new_structure = QEStructure(qeInput = self._qeInput)
            new_structure._setStructureFromDiffpyStructure(diffpyStruct, \
                                        massList = [], psList = [], ibrav = 0)

        new_structure.lattice._qeInput.update()
        self.__Init(new_structure)
        return parser

        
    def readStr(self): pass


    def load(self, source, **args):
        #from qecalc.qetask.qeparser.qesrtuctureparser import *
        task = {
            'diffpy': self._setStructureFromDiffpyStructure,
 #           'pwinput': self._setStructureFromPWConfig,
 #           'pwoutput': self._setStructureFromPWOutput,
        }
        if source == 'diffpy':
            if 'ibrav' in args and args['ibrav'] != 0:
                task['diffpy'] = self._setReducedStructureFromDiffpyStructure

        task[source](**args)
        
        self.lattice._qeInput.update()


    def _setStructureFromDiffpyStructure(self, structure, massList = [], psList = [], ibrav = 0):
        """
        structure - diffpy.Structure object
        ibrav - Lattice index
        psList - list of strings with pseudopotential names
        diffpyStructure object will be modified with reduced atomic positions
        """      
        diffpyLattice = structure.lattice
        
        self.structure = structure
        
        self.atomicSpecies = OrderedDict()
               
        
        #set lattice and  convert to bohr units
        #qeLattice = QELattice(ibrav = 0, a = 1.889725989, base = diffpyLattice.base)
        qeLattice = QELattice(ibrav = 0, base = diffpyLattice.base)
        qeLattice.a = 1.889725989*qeLattice.a
        qeLattice._qeInput = self._qeInput
        
        self.lattice = qeLattice
        self.lattice.type = 'generic cubic'

        atomNames = []
        for a in structure:
            if self._element(a) not in atomNames:
                atomNames.append(self._element(a))
        
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
     
                        
    def _setReducedStructureFromDiffpyStructure(self, structure, ibrav, massList = [], psList = []):
        """
        structure - diffpy.Structure object
        ibrav - Lattice index
        psList - list of strings with pseudopotential names
        diffpyStructure object will be modified with reduced atomic positions
        """
        import copy

        #self.atomicSpecies = OrderedDict()
        #self.optConstraints = []
        #self.atomicPositionsType = 'crystal'

        diffpyLattice = copy.deepcopy(structure.lattice)

        self.atomicSpecies = OrderedDict()
        
        a = diffpyLattice.a
        b = diffpyLattice.b
        c = diffpyLattice.c
        cAB = cosd(diffpyLattice.gamma)
        cBC = cosd(diffpyLattice.alpha)
        cAC = cosd(diffpyLattice.beta)

        qeLattice = QELattice(ibrav = ibrav, a = a, b = b, c = c,  cBC =  cBC, \
                              cAC = cAC, cAB = cAB)

        qeLattice._qeInput = self._qeInput
        self.lattice = qeLattice
        # make a deep copy:
        reducedStructure = Structure(atoms = structure)
        #reducedStructure = structure

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


    def toString(self, string = None):
        if string != None:
            string = self.lattice.toString(string = string)
            qeInput = QEInput(config = string)
            qeInput.parse()
        else:
            if self.lattice._qeInput != None:
                qeInput = self.lattice._qeInput
            else:
                qeInput = QEInput(config = '')

        #self.updatePWInput(qeInput)
        return qeInput.toString()

 
    def save(self, fname = None):
        """Writes/updates structure into PW config file,
           if the file does not exist, new one will be created"""
        filename = fname
        if fname != None:
            self.lattice.save(filename)
            qeInput = QEInput(fname)
            qeInput.parse()
        else:
            filename = self.filename
            self.lattice.save(filename)
            qeInput = self.lattice._qeInput
            
        #self.updatePWInput(qeInput)
        qeInput.save(filename)
     
    
    def updatePWInput(self, qeInput = None):
        """
        Deprecated
        """
        self._qeInput.update()

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

if __name__ == '__main__': pass
