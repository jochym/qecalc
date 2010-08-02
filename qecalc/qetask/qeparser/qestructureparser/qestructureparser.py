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

class QEStructureParser():
    def __init__(self, qeConf):
        self.qeConf = qeConf
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
 
 
    def parseQEConf(self, qeConf = None):
        """ Loads structure from PWSCF config class"""
        
        if qeConf != None:
            self.qeConf = qeConf
        self.qeConf.autoUpdate = False
        stru = QEStructure(qeConf = self.qeConf)
        
        # make qeConf consistent with the current instance of the structure
        stru._qeConf.structure = stru        
                
        stru.atomicSpecies = OrderedDict()
        stru.lattice = self.__getLattice(self.qeConf)               
        
        stru.structure = Structure(lattice = stru.lattice.diffpy())       
        stru.nat = stru.ntyp = None
        #self.filename = self.qeConf.filename
        stru.optConstraints = []
        if 'system' in stru.lattice.qeConf.namelists:
            stru.nat  = int(stru.lattice.qeConf.namelist('system').param('nat'))
            stru.ntyp  = int(stru.lattice.qeConf.namelist('system').param('ntyp'))
        if 'atomic_positions' in stru.lattice.qeConf.cards:        
            atomicLines = stru.lattice.qeConf.card('atomic_positions').lines()
            stru.atomicPositionsType = stru.lattice.qeConf.card('atomic_positions').arg()
#            if self.atomicPositionsType == 'angstrom':
#                raise NotImplementedError\
#         ('atomic positions in bohr and angstrom are not currently supported')
            if stru.atomicPositionsType == None:
                stru.atomicPositionsType = 'alat'
            for line in atomicLines:
                if '!' not in line:
                    words = line.split()
                    coords = [float(w) for w in words[1:4]]
                    constraint = []
                    if len(words) > 4:
                        constraint = [int(c) for c in words[4:7]]
                    stru.optConstraints.append(numpy.array(constraint, dtype = int))
                    atomSymbol = words[0]
                    if stru.atomicPositionsType == 'alat':
                        coords = stru.lattice.diffpy().fractional(numpy.array(coords[0:3])*stru.lattice.a)
                    if stru.atomicPositionsType == 'crystal':
                        coords = numpy.array(coords[0:3])
                    if stru.atomicPositionsType == 'bohr' or stru.atomicPositionsType == 'angstrom':
                        coords = stru.lattice.diffpy().fractional(numpy.array(coords[0:3]))
                    stru.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))
        # parse mass ATOMIC_SPECIES section:
         
        if 'atomic_species' in stru.lattice.qeConf.cards:
            atomicSpeciesLines = stru.lattice.qeConf.card('atomic_species').lines()
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
                        stru.atomicSpecies[element] =  AtomicSpecies(element, mass, ps)
        self.qeConf.autoUpdate = True                
        return stru
                        

    def __getLattice(self, qeConf ):

        if qeConf == None:
            raise QELatticeError("__getLattice: qeConf was not properly initialized")      
          
        lat = QELattice()
        lat.qeConf = qeConf
        
        # make qeConf consistent with the current instance of the lattice
        lat.qeConf.structure.lattice = lat
           
             
        if 'ibrav' in lat.qeConf.namelists['system'].params:
            ibrav  = int(lat.qeConf.namelist('system').param('ibrav'))            
        else:
            raise QELatticeError("config file should have ibrav defined")
        if ibrav < 0:            
            raise QELatticeError("ibrav should be integer >= 0")
                    
        #************************************************
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        if 'celldm(1)' in qeConf.namelists['system'].params:
            lat._type = 'celldm'
            a = float(qeConf.namelist('system').param('celldm(1)'))

            if ibrav == 0:
                # lattice is set in the units of celldm(1)
                # need to parse CELL_PARAMETERS
                #if 'cell_parameters' not in qeConf.cards:
                #    return  #qeConf.createCard('cell_parameters')
                cellParLines = qeConf.card('cell_parameters').lines()
                #print cellParLines
                cellParType = qeConf.card('cell_parameters').arg()
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
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                latPar = [a, a, c_a*a, cBC, cAC, cAB, None]
            if ibrav == 5:
                cAB = float(qeConf.namelist('system').param('celldm(4)'))
                latPar = [a, a, a, cAB, cAB, cAB, None]
            if ibrav > 7 and ibrav < 12:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                latPar = [a, b_a*a, c_a*a, cBC, cAC, cAB, None]
            if ibrav == 12 or ibrav == 13:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                cAB = float(qeConf.namelist('system').param('celldm(4)'))
                latPar = [a, b_a*a, c_a*a, cBC, cAC, cAB, None]
            if ibrav == 14:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                cBC = float(qeConf.namelist('system').param('celldm(4)'))
                cAC = float(qeConf.namelist('system').param('celldm(5)'))
                cAB = float(qeConf.namelist('system').param('celldm(6)'))
                latPar = [a, b_a*a, c_a*a, cBC, cAC, cAB, None]
        else:
            if ibrav == 0:
                raise QELatticeError("Should specify celldm(1) if use 'generic' lattice")
            a = float(qeConf.namelist('system').param('A'))
            lat._type = 'traditional'   # A, B, C, cosAB, cosAC, cosBC
            if ibrav > 0 and ibrav < 4:
                latPar = [a, a, a, cBC, cAC, cAB, None]
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c = float(qeConf.namelist('system').param('C'))
                latPar = [a, a, c, cBC, cAC, cAB, None]
            if ibrav == 5:
                cAB = float(qeConf.namelist('system').param('cosAB'))
                latPar = [a, a, a, cAB, cAB, cAB, None]
            if ibrav > 7 and ibrav < 12:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                latPar = [a, b, c, cBC, cAC, cAB, None]
            if ibrav == 12 or ibrav == 13:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                cAB = float(qeConf.namelist('system').param('cosAB'))
                latPar = [a, b, c, cBC, cAC, cAB, None]
            if ibrav == 14:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                cBC = float(qeConf.namelist('system').param('cosBC'))
                cAC = float(qeConf.namelist('system').param('cosAC'))
                cAB = float(qeConf.namelist('system').param('cosAB'))
                latPar = [a, b, c, cBC, cAC, cAB, None]               
        
        a, b, c, cBC, cAC, cAB, base = latPar
        lat.setLattice(ibrav, a, b, c, cBC, cAC, cAB, base, updateInput = False)
        
        return lat
                                           