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

from qeinput import QEInput
from qestructure import QEStructure
from pwkpoints import PWKpoints

class PWInput(QEInput):
    def __init__(self, filename=None, config=None):
        QEInput.__init__(self,filename, config, type='pw')
        self.structure = QEStructure(self)
        self.kpoints = PWKpoints(self)
    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary
            Initializes structure as well"""
        (self.header, self.namelists, self.cards, self.attach) = self.parser.parse()
        self.structure.parseInput(self)
        self.kpoints.parse()

#    def toString(self):       
#        s = ''
#        namelistOrder = {
#        0 : 'control',
#        1 : 'system',
#        2 : 'electrons',
#        3 : 'ions',
#        4 : 'cell',
#        5 : 'phonon',
#        6 : 'ee'
#        }
#
#        cardOrder = {
#        0 : 'atomic_species',
#        1 : 'atomic_positions',
#        2 : 'k_points',
#        3 : 'cell_parameters',
#        4:  'occupations',
#        5 : 'climbing_images',
#        6 : 'constraints',
#        7 : 'collective_vars'
#        }
#        for i in range(len(namelistOrder)):
#            if namelistOrder[i] in self.namelists:
#                s += self.namelists[namelistOrder[i]].toString()
#                s += '\n'
#
#        for i in range(len(cardOrder)):
#            if cardOrder[i] in self.cards:
#                s += self.cards[cardOrder[i]].toString()
#                s += '\n'
#        return s

    def _update(self, qeConf = None):
        if qeConf == None:
            qeConf = self        

        lattice = self.structure.lattice
        structure = self.structure   

        #************************* updateLattice:*******************************
        if 'system' not in qeConf.namelists:
            qeConf.createNamelist('system')
        # clear geometry from qeConf:
        qeConf.namelist('system').remove('a')
        qeConf.namelist('system').remove('b')
        qeConf.namelist('system').remove('c')
        qeConf.namelist('system').remove('cosab')
        qeConf.namelist('system').remove('cosbc')
        qeConf.namelist('system').remove('cosac')
        qeConf.namelist('system').remove('celldm(1)')
        qeConf.namelist('system').remove('celldm(2)')
        qeConf.namelist('system').remove('celldm(3)')
        qeConf.namelist('system').remove('celldm(4)')
        qeConf.namelist('system').remove('celldm(5)')
        qeConf.namelist('system').remove('celldm(6)')
        if 'cell_parameters' in qeConf.cards:
            qeConf.removeCard('cell_parameters')
        if lattice._type == 'celldm':
            qeConf.namelist('system').add('ibrav', lattice._ibrav)
            qeConf.namelist('system').add('celldm(1)', lattice._a)
            qeConf.namelist('system').add('celldm(2)', float(lattice._b)/lattice._a)
            qeConf.namelist('system').add('celldm(3)', float(lattice._c)/lattice._a)
            if lattice._ibrav < 14:
                qeConf.namelist('system').add('celldm(4)', lattice._cAB)
            else:
                qeConf.namelist('system').add('celldm(4)', lattice._cBC)
                qeConf.namelist('system').add('celldm(5)', lattice._cAC)
                qeConf.namelist('system').add('celldm(6)', lattice._cAB)
        else:
            if self._type == 'traditional':
                qeConf.namelist('system').add('ibrav', lattice._ibrav)
                qeConf.namelist('system').add('A', lattice._a)
                qeConf.namelist('system').add('B', lattice._b)
                qeConf.namelist('system').add('C', lattice._c)
                qeConf.namelist('system').add('cosAB', lattice._cAB)
                qeConf.namelist('system').add('cosAC', lattice._cAC)
                qeConf.namelist('system').add('cosBC', lattice._cBC)
            else:
                if 'generic' in lattice._type:
                    qeConf.namelist('system').add('celldm(1)', lattice._a)
                    lattice._ibrav = 0
                    qeConf.namelist('system').add('ibrav', lattice._ibrav)
                    if lattice._type == 'generic hexagonal':
                        cardArg = 'hexagonal'
                    if lattice._type == 'generic cubic' or lattice._type == None:
                        cardArg = 'cubic'
                    qeConf.createCard('cell_parameters')
                    qeConf.card('cell_parameters').setArg(cardArg)
                    qeConf.card('cell_parameters').removeLines()
                    for i in range(3):
                        v = lattice._primitiveLattice.base[i,:]/float(lattice._a)
                        qeConf.card('cell_parameters').addLine(\
                                       lattice.formatString%(v[0], v[1], v[2]))
        
        #************************* updateStructure:*****************************
        
        if 'system' not in qeConf.namelists:
            qeConf.addNamelist('system')            
        qeConf.namelist('system').remove('ntyp')
        qeConf.namelist('system').remove('nat')
        if structure.ntyp != None:
            qeConf.namelist('system').add('ntyp', structure.ntyp)
        if structure.nat != None:
            qeConf.namelist('system').add('nat', structure.nat)
        
        if len(qeConf.namelist('system').params) == 0:
            qeConf.removeNamelist('system')  

        if 'atomic_positions' in qeConf.cards:
            qeConf.removeCard('atomic_positions')
        qeConf.createCard('atomic_positions')
        qeConf.card('atomic_positions').setArg(structure.atomicPositionsType)
        for atom, constraint in zip(structure.structure, structure.optConstraints):
            if structure.atomicPositionsType == 'alat':
                coords = structure.lattice.diffpy().cartesian(atom.xyz)/structure.lattice.a
                coords = structure.formatString%(coords[0], coords[1], coords[2])
            else:
                if structure.atomicPositionsType == 'crystal':
                    coords = structure.formatString%(atom.xyz[0], atom.xyz[1], atom.xyz[2])
                else:
                    if structure.atomicPositionsType == 'bohr' or structure.atomicPositionsType == 'angstrom':
                        coords = structure.lattice.diffpy().cartesian(atom.xyz)
                        coords = structure.formatString%(coords[0], coords[1], coords[2])
                    else:
                        raise NonImplementedError
            line = '%-3s'%structure._element(atom) + '    ' + coords + '  ' + str(constraint)[1:-1]
            qeConf.card('atomic_positions').addLine(line)
        
        if len(qeConf.card('atomic_positions').lines()) == 0:
            qeConf.removeCard('atomic_positions')

        # update ATOMIC_SPECIES card
        if 'atomic_species' in qeConf.cards:
            qeConf.removeCard('atomic_species')
        qeConf.createCard('atomic_species')
        for element, specie in structure.atomicSpecies.items():
            qeConf.card('atomic_species').addLine(specie.toString())
        
        if len(qeConf.card('atomic_species').lines()) == 0:
            qeConf.removeCard('atomic_species')                


    def outDir(self):
        self.parse()
        return self.namelist('control').param('outdir', quotes = False)

textA = """
&control
calculation='scf'
! restart_mode='from_scratch',
wf_collect = .true.,
tstress = .true. ,
tprnfor = .true. ,
verbosity = 'high',
prefix='mgb2',
pseudo_dir = '/home/markovsk/projects/pslib/espresso/mgb2/',
lkpoint_dir = .false. ,
outdir='temp/'
/
&system
ibrav=4,
celldm(1) = 5.78739785,
celldm(2) = 5.78739785,
celldm(3) = 1.135794331,
! celldm(1) = 5.8260,
! celldm(2) = 5.8260,
! celldm(3) = 1.1420,
nat = 3,
ntyp = 2,
nspin = 1,
nbnd = 12,
occupations='smearing',
! degauss=0.025
degauss=0.025,
smearing = 'methfessel-paxton' ,
ecutwfc =32.0,
ecutrho =256.0,
la2f = .false.
/
&electrons
conv_thr = 1.0d-12
diago_full_acc=.TRUE.
/
ATOMIC_SPECIES
Mg 24.305 mg_6.ncpp
B 11.000 B.pbe-n-van_ak.UPF
!b_rc_1.4_pcc.ncpp
! B 10.811 B.pbe-tmnc.UPF


ATOMIC_POSITIONS alat
Mg 0.000000000 0.0000000000000000 0.000000000
B 0.500000000 0.2886751345948129 0.5678971655
B 0.000000000 0.5773502691896257 0.5678971655


K_POINTS AUTOMATIC
24 24 24 0 0 0
"""

if __name__ == "__main__":
    inp = PWInput(config=textA)
    inp.parse()
    print inp.toString()

__author__="Nikolay Markovskiy"
__date__ ="$Oct 20, 2009 12:21:58 PM$"
