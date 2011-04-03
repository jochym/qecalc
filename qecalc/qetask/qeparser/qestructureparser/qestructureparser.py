#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from qecalc.qetask.qeparser.qestructureparser import *
from qecalc.qetask.qeparser.orderedDict import OrderedDict

class QEStructureParser():
    def __init__(self, qeInput):
        self._qeInput = qeInput
        self.format = None
        self.filename = None
    
    def parseStr(self, s):
        """Create Structure instance from a string."""
        raise NotImplementedError, \
                "parseStr not defined for '%s' format" % self.format
         

    def parse(self, filename):
        """Create Structure instance from a file."""        
        raise NotImplementedError, \
                "parse not defined for '%s' format" % self.format
 
 
    def parseqeInput(self, qeInput = None):
        """ Loads structure from PWSCF config class"""
        
        if qeInput != None:
            self._qeInput = qeInput
        autoUpdate = self._qeInput.autoUpdate
        self._qeInput.autoUpdate = False
        stru = QEStructure(qeInput = self._qeInput)
        
        # make qeInput consistent with the current instance of the structure
        stru._qeInput.structure = stru        
                
        stru.lattice = self.__getLattice(self._qeInput)     
        
        noncolin = False
        if 'noncolin' in stru._qeInput.namelists['system'].paramlist():
            noncolin =  stru._qeInput.namelist('system').get('noncolin').lower()
            if noncolin == '.true.':
                noncolin = True
            else:
                noncolin = False  
        
        lda_plus_u = False
        if 'lda_plus_u' in stru._qeInput.namelists['system'].paramlist():
            lda_plus_u =  stru._qeInput.namelist('system').get('lda_plus_u').lower()            
            if lda_plus_u == '.true.':
                lda_plus_u = True
            else:
                lda_plus_u = False
                
        nspin = 1
        if 'nspin' in stru._qeInput.namelists['system'].paramlist():
            nspin =  int(stru._qeInput.namelist('system').get('nspin').lower())                               
        
        if 'atomic_positions' in stru.lattice._qeInput.cards:        
            atomicLines = stru.lattice._qeInput.card('atomic_positions').lines()
            stru.atomicPositionsType = stru.lattice._qeInput.card('atomic_positions').arg()
            if stru.atomicPositionsType == None:
                stru.atomicPositionsType = 'alat'
            for line in atomicLines:
                if '!' not in line:
                    words = line.split()
                    coords = [float(w) for w in words[1:4]]
                    constraint = []
                    if len(words) > 4:
                        constraint = [int(c) for c in words[4:7]]
                    atomSymbol = words[0]
                    if stru.atomicPositionsType == 'alat':
                        coords = stru.lattice.matter().fractional(numpy.array(coords[0:3])*stru.lattice.a)
                    if stru.atomicPositionsType == 'crystal':
                        coords = numpy.array(coords[0:3])
                    if stru.atomicPositionsType == 'bohr' or stru.atomicPositionsType == 'angstrom':
                        coords = stru.lattice.matter().fractional(numpy.array(coords[0:3]))
                    stru.addNewAtom(atype = atomSymbol, xyz = numpy.array(coords[0:3]), \
                                    optConstraint = numpy.array(constraint, dtype = int))
        # parse mass ATOMIC_SPECIES section:
        atomicSpecies = OrderedDict()
        # default values:
        for a in stru:
            atomicSpecies[a.symbol] = (0, '')
        if 'atomic_species' in stru.lattice._qeInput.cards:
            atomicSpeciesLines = stru.lattice._qeInput.card('atomic_species').lines()
            for line in atomicSpeciesLines:
                if '!' not in line:
                    if line.strip() != '':                     
                        atomicSpeciesWords = line.split()
                        symbol = atomicSpeciesWords[0]
                        mass = 0
                        ps = ''
                        if len(atomicSpeciesWords) > 1 :
                            mass = float(atomicSpeciesWords[1])
                        if len(atomicSpeciesWords) > 2:
                            ps = atomicSpeciesWords[2]
                        atomicSpecies[symbol] =  (float(mass), ps)
        
        for a in stru:
            mass = atomicSpecies[a.symbol][0]
            ps  = atomicSpecies[a.symbol][1]
            a.mass = mass
            a.potential = ps
        self._qeInput.autoUpdate = autoUpdate                
        return stru
        
        
    def __getNamelistParams(self, input, namelist, param):
        """
        Extracts data from parameters with parentheses
        namelist - name of namelist
        param - parameter name without parentheses        
        """
        paraDic = {}
        for par in input.namelists[namelist].paramlist():
            if param in par and '(' in par:
                index = int(par.strip(')')[0].strip('(')[1])
                paraDic[index] = input.namelist(namelist).get(par)
        return paraDic

    def __getLattice(self, qeInput ):

        if qeInput == None:
            raise QELatticeError("__getLattice: qeInput was not properly initialized")      
          
        lat = QELattice()
        lat._qeInput = qeInput
        
        # make qeInput consistent with the current instance of the lattice
        lat._qeInput.structure.lattice = lat
           
             
        if 'ibrav' in lat._qeInput.namelists['system'].paramlist():
            ibrav  = int(lat._qeInput.namelist('system').get('ibrav'))            
        else:
            raise QELatticeError("config file should have ibrav defined")
        if ibrav < 0:            
            raise QELatticeError("ibrav should be integer >= 0")
                    
        #************************************************
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        if 'celldm(1)' in qeInput.namelists['system'].paramlist():
            lat._type = 'celldm'
            a = float(qeInput.namelist('system').get('celldm(1)'))

            if ibrav == 0:
                # lattice is set in the units of celldm(1)
                # need to parse CELL_PARAMETERS
                #if 'cell_parameters' not in qeInput.cards:
                #    return  #qeInput.createCard('cell_parameters')
                cellParLines = qeInput.card('cell_parameters').lines()
                cellParType = qeInput.card('cell_parameters').arg()
                if cellParType == 'cubic' or cellParType == None:
                    lat._type = 'generic cubic'
                else:
                    if cellParType == 'hexagonal':
                        lat._type = 'generic hexagonal'
                # convert card into list
                base = []
                for line in cellParLines:
                    if '!' not in line:
                        words = line.split()
                        base.append([float(w) for w in words])
                latPar =  [a, None, None, None, None, None, numpy.array(base)*a]
            if ibrav > 0 and ibrav < 4:
                latPar = [a, a, a, cBC, cAC, cAB, None]

            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c_a = float(qeInput.namelist('system').get('celldm(3)'))
                latPar = [a, a, c_a*a, cBC, cAC, cAB, None]
            if ibrav == 5:
                cAB = float(qeInput.namelist('system').get('celldm(4)'))
                latPar = [a, a, a, cAB, cAB, cAB, None]
            if ibrav > 7 and ibrav < 12:
                b_a = float(qeInput.namelist('system').get('celldm(2)'))
                c_a = float(qeInput.namelist('system').get('celldm(3)'))
                latPar = [a, b_a*a, c_a*a, cBC, cAC, cAB, None]
            if ibrav == 12 or ibrav == 13:
                b_a = float(qeInput.namelist('system').get('celldm(2)'))
                c_a = float(qeInput.namelist('system').get('celldm(3)'))
                cAB = float(qeInput.namelist('system').get('celldm(4)'))
                latPar = [a, b_a*a, c_a*a, cBC, cAC, cAB, None]
            if ibrav == 14:
                b_a = float(qeInput.namelist('system').get('celldm(2)'))
                c_a = float(qeInput.namelist('system').get('celldm(3)'))
                cBC = float(qeInput.namelist('system').get('celldm(4)'))
                cAC = float(qeInput.namelist('system').get('celldm(5)'))
                cAB = float(qeInput.namelist('system').get('celldm(6)'))
                latPar = [a, b_a*a, c_a*a, cBC, cAC, cAB, None]
        else:
            if ibrav == 0:
                raise QELatticeError("Should specify celldm(1) if use 'generic' lattice")
            a = float(qeInput.namelist('system').get('A'))
            lat._type = 'traditional'   # A, B, C, cosAB, cosAC, cosBC
            if ibrav > 0 and ibrav < 4:
                latPar = [a, a, a, cBC, cAC, cAB, None]
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c = float(qeInput.namelist('system').get('C'))
                latPar = [a, a, c, cBC, cAC, cAB, None]
            if ibrav == 5:
                cAB = float(qeInput.namelist('system').get('cosAB'))
                latPar = [a, a, a, cAB, cAB, cAB, None]
            if ibrav > 7 and ibrav < 12:
                b = float(qeInput.namelist('system').get('B'))
                c = float(qeInput.namelist('system').get('C'))
                latPar = [a, b, c, cBC, cAC, cAB, None]
            if ibrav == 12 or ibrav == 13:
                b = float(qeInput.namelist('system').get('B'))
                c = float(qeInput.namelist('system').get('C'))
                cAB = float(qeInput.namelist('system').get('cosAB'))
                latPar = [a, b, c, cBC, cAC, cAB, None]
            if ibrav == 14:
                b = float(qeInput.namelist('system').get('B'))
                c = float(qeInput.namelist('system').get('C'))
                cBC = float(qeInput.namelist('system').get('cosBC'))
                cAC = float(qeInput.namelist('system').get('cosAC'))
                cAB = float(qeInput.namelist('system').get('cosAB'))
                latPar = [a, b, c, cBC, cAC, cAB, None]               
        
        a, b, c, cBC, cAC, cAB, base = latPar
        lat.setLattice(ibrav, a, b, c, cBC, cAC, cAB, base, updateInput = False)
        
        return lat
                                           