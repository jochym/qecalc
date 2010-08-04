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

from qeatom import QEAtom

try:
    from diffpy.Structure.structure import Structure
    from diffpy.Structure.lattice import cosd, Lattice
    from diffpy.Structure.SymmetryUtilities import equalPositions
except ImportError:
    from matter import Structure, Lattice
    from matter.Lattice import cosd
    from matter.SymmetryUtilities import equalPositions

from qelattice import QELattice
import numpy
from qeinput import QEInput
from orderedDict import OrderedDict


class QEStructure( Structure ):
    
    def __init__(self, qeInput = None, filename = None):
        """the structure is initialized from PWSCF config file
           'lattice' and 'structure' are automatically updated"""
           
        Structure.__init__(self)
        self.formatString = '%# .8f %# .8f %# .8f'
        # optConstraints three 1/0 for each coordinate of each atom
        self._optConstraints = []
        self.lattice = QELattice()
        
        if filename != None:
            from qecalc.qetask.qeparser.pwinput import PWInput
            qeInput = PWInput(filename = filename)
            qeInput.parse()
            self.parseInput(qeInput)
            return
        
        self.lattice._qeInput = qeInput
        self._qeInput = qeInput
        if qeInput != None:
            self.lattice._qeInput.structure = self
            self.lattice._qeInput.structure.lattice = self.lattice            
        self._atomicPositionsType = 'crystal'                           
         
         
    def addNewAtom(self, *args, **kwargs):
        """Add new Atom instance to the end of this Structure.

        All arguments are forwarded to Atom constructor.

        No return value.
        """
        kwargs['lattice'] = self.lattice
        a = QEAtom(*args, **kwargs)
        list.append(self, a)
        self._uncache('labels')
        return
    

    def placeInLattice(self, new_lattice):
        """place structure into new_lattice coordinate system

        sets lattice to new_lattice and recalculate fractional coordinates
        of all atoms so their absolute positionls remain the same

        return self
        """
        Tx = numpy.dot(self.lattice.diffpy().base, new_lattice.diffpy().recbase)
        Tu = numpy.dot(self.lattice.diffpy().normbase, new_lattice.diffpy().recnormbase)
        for a in self:
            a.xyz = numpy.dot(a.xyz, Tx)
        tmpInput = self.latice._qeInput
        self.lattice = new_lattice
        self.lattice._qeInput = tmpInput
        self.lattice._qeInput.structure = self        
        return self
        

    ##############################################################################
    # overloaded list methods - taken from diffpy.Structure
    ##############################################################################

    def append(self, a, copy=True):
        """Append atom to a structure and update its lattice attribute.

        a    -- instance of QEAtom
        copy -- flag for appending a copy of a.
                When False, append a and update a.owner.

        No return value.
        """
        self._uncache('labels')
        adup = copy and QEAtom(a) or a
        adup.lattice = self.lattice
        list.append(self, adup)
        return

    def insert(self, idx, a, copy=True):
        """Insert atom a before position idx in this Structure.

        idx  -- position in atom list
        a    -- instance of QEAtom
        copy -- flag for inserting a copy of a.
                When False, append a and update a.lattice.

        No return value.
        """
        self._uncache('labels')
        adup = copy and QEAtom(a) or a
        adup.lattice = self.lattice
        list.insert(self, idx, adup)
        return

    def extend(self, atoms, copy=True):
        """Extend Structure by appending copies from a list of atoms.

        atoms -- list of QEAtom instances
        copy  -- flag for extending with copies of QEAtom instances.
                 When False extend with atoms and update their lattice
                 attributes.

        No return value.
        """
        self._uncache('labels')
        if copy:    adups = [QEAtom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.extend(self, adups)
        return

    def __setitem__(self, idx, a, copy=True):
        """Set idx-th atom to a.

        idx  -- index of atom in this Structure
        a    -- instance of Atom
        copy -- flag for setting to a copy of a.
                When False, set to a and update a.lattice.

        No return value.
        """
        self._uncache('labels')
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.__setitem__(self, idx, adup)
        return

    def __setslice__(self, lo, hi, atoms, copy=True):
        """Set Structure slice from lo to hi-1 to the sequence of atoms.

        lo    -- low index for the slice
        hi    -- high index of the slice
        atoms -- sequence of Atom instances
        copy  -- flag for using copies of Atom instances.  When False, set
                 to existing instances and update their lattice attributes.

        No return value.
        """
        self._uncache('labels')
        if copy:    adups = [QEAtom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.__setslice__(self, lo, hi, adups)
        return
    
            
    def _get_nat(self):
        return len(self)
    nat = property(_get_nat, doc ="number of atoms")
    
 
    def _get_ntyp(self):
        return len(self.atomicSpecies)
    
    ntyp = property(_get_ntyp, doc ="number of types")
    
    
    def _get_atomicPositionsType(self):
        return self._atomicPositionsType

    def _set_atomicPositionsType(self, value):
        self._atomicPositionsType = value
        self.lattice._qeInput.update()

    atomicPositionsType = property(_get_atomicPositionsType, \
                                   _set_atomicPositionsType, \
                           doc ="type of atomic positions (crystal, alat ...)")    


    def _get_atomicSpecies(self):
        atomicSpecies = OrderedDict()
        for a in self:
            atomicSpecies[a.element] = AtomicSpecies(element = a.element, \
                                        mass = a.mass, potential = a.potential)
        
        return atomicSpecies

    atomicSpecies = property(_get_atomicSpecies, \
              doc ="returns an ordered dictionary with atomic species' objects")   


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
        for atom, constraint in zip(self, self.optConstraints):
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
        from qestructureparser.qestructureparser import QEStructureParser
        new_structure = QEStructureParser(qeInput).parseqeInput()
        self.__Init(new_structure)
    
    
    def __Init(self, structure):
        QEStructure.__init__(self)
        if structure is not None:
            self.__dict__.update(structure.__dict__)
            self.lattice.__dict__.update(structure.lattice.__dict__)
            self[:] = structure     


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

        
    def readStr(self, s, format = 'pwinput'):
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
            new_structure = parser.parseStr(s)
        else:            
            diffpyStruct = Structure()
            parser = diffpyStruct.readStr(s, format = format)
            new_structure = QEStructure(qeInput = self._qeInput)
            new_structure._setStructureFromDiffpyStructure(diffpyStruct, \
                                        massList = [], psList = [], ibrav = 0)

        new_structure.lattice._qeInput.update()
        self.__Init(new_structure)
        return parser


    def load(self, source, **args):
        task = {
            'diffpy': self._setStructureFromDiffpyStructure,
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
        psList - list of strings with potential names
        diffpyStructure object will be modified with reduced atomic positions
        """      
        diffpyLattice = structure.lattice
            
        qeLattice = QELattice(ibrav = 0, base = diffpyLattice.base)
        qeLattice.a = 1.889725989*qeLattice.a
        qeLattice._qeInput = self._qeInput
        
        self.lattice = qeLattice
        self.lattice.type = 'generic cubic'

        atomNames = []        
        for a in structure:
            if self._element(a) not in atomNames:
                atomNames.append(self._element(a))
        
        atomicSpecies = {}
        for i, elem in enumerate(atomNames):
            if len(massList) - 1 < i:
                mass = 0
            else:
                mass = massList[i]
            if len(psList) - 1 < i:
                ps = ''
            else:
                ps = psList[i]               
            atomicSpecies[elem] =  (mass, ps)
        
        for atom in structure:
            self.optConstraints.append([])       
        
        self[:] = []
        for atom in structure:
            elem = self._element(atom)
            self.addNewAtom(atype = elem, xyz = atom.xyz, \
                            mass = atomicSpecies[elem][0], \
                            potential = atomicSpecies[elem][1],\
                            lattice = self.lattice)             
             
                        
    def _setReducedStructureFromDiffpyStructure(self, structure, ibrav, massList = [], psList = []):
        """
        structure - diffpy.Structure object
        ibrav - Lattice index
        psList - list of strings with potential names
        diffpyStructure object will be modified with reduced atomic positions
        """
        import copy

        diffpyLattice = copy.deepcopy(structure.lattice)

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

        reducedStructure.placeInLattice(Lattice(base=qeLattice.diffpy().base))

        # collect atoms that are at equivalent position to some previous atom
        duplicates = set([a1
            for i0, a0 in enumerate(reducedStructure) for a1 in reducedStructure[i0+1:]
                if   self._element(a0) == self._element(a1) and equalPositions(a0.xyz, a1.xyz, eps=1e-4)])

        
        # Filter out duplicate atoms.  Use slice assignment so that
        # reducedStructure is not replaced with a list.
        reducedStructure[:] = [a for a in reducedStructure if not a in duplicates]

        atomNames = []
        for a in reducedStructure:
            if self._element(a) not in atomNames:
                atomNames.append(self._element(a))
        
        atomicSpecies = {}
        for i, elem in enumerate(atomNames):
            if len(massList) - 1 < i:
                mass = 0
            else:
                mass = massList[i]
            if len(psList) - 1 < i:
                ps = ''
            else:
                ps = psList[i]
            atomicSpecies[elem] =  (mass, ps)
               
        self[:] = []
        
        # convert to bohr units
        self.lattice.setLattice(ibrav, self.lattice.a*1.889725989, \
                                 self.lattice.b*1.889725989,
                                 self.lattice.c*1.889725989)

        for atom in reducedStructure:
            self.optConstraints.append([])
        for atom in reducedStructure:
            elem = self._element(atom)
            self.addNewAtom(atype = elem, xyz = atom.xyz, \
                            mass = atomicSpecies[elem][0], \
                            potential = atomicSpecies[elem][1],\
                            lattice = self.lattice)


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
            
        qeInput.save(filename)
     
    
    def updatePWInput(self, qeInput = None):
        """
        Deprecated
        """
        self._qeInput.update()

    def diffpy(self):
        stru = Structure(lattice = self.lattice.diffpy())
        for atom in self:
            stru.addNewAtom(atype = atom.element, xyz = atom.xyz, \
                                              lattice = self.lattice.diffpy() )
        return stru


    def _element(self, atom):
        """
        Is needed for support both diffpy and matter classes
        """
        if 'element' in dir(atom):
            return atom.element
        else:
            if 'symbol' in dir(atom):
                return atom.symbol
            else:
                raise
            
class AtomicSpecies():
    def __init__(self, element = 'H', mass = 1.0, potential = ''):
        self.element = element
        self.potential = potential
        self.mass = mass
    def __str__(self):
        return '%-3s'%self.element + ' ' + '%.4f'%self.mass + ' ' + self.potential
    def toString(self):
        return str(self)              

if __name__ == '__main__': pass
