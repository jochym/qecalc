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
        # Boolean flag, if True, QEInput is updated on change of any property in
        # structure, lattice, or atom
        self.autoUpdate = True         
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

    def update(self, qeInput = None):
        """
        Loads current mathematical representation of Structure and Lattice 
        into  QEInput parsing object. 
        
        Will be deprecated soon. Currently it is already automatically called \
        for most  of Lattice  and structure operations.
        E.g. 'input.structure.lattice.a = new_a'   will automatically update
        input       
        """
        if qeInput == None:
            qeInput = self        
        if self.autoUpdate is False:
            return
        lattice = self.structure.lattice
        structure = self.structure   

        #************************* updateLattice:*******************************
        if 'system' not in qeInput.namelists:
            qeInput.createNamelist('system')
        # clear geometry from qeInput:
        qeInput.namelist('system').remove('a')
        qeInput.namelist('system').remove('b')
        qeInput.namelist('system').remove('c')
        qeInput.namelist('system').remove('cosab')
        qeInput.namelist('system').remove('cosbc')
        qeInput.namelist('system').remove('cosac')
        qeInput.namelist('system').remove('celldm(1)')
        qeInput.namelist('system').remove('celldm(2)')
        qeInput.namelist('system').remove('celldm(3)')
        qeInput.namelist('system').remove('celldm(4)')
        qeInput.namelist('system').remove('celldm(5)')
        qeInput.namelist('system').remove('celldm(6)')
        if 'cell_parameters' in qeInput.cards:
            qeInput.removeCard('cell_parameters')
        if lattice._type == 'celldm':
            qeInput.namelist('system').add('ibrav', lattice._ibrav)
            qeInput.namelist('system').add('celldm(1)', lattice._a)
            qeInput.namelist('system').add('celldm(2)', float(lattice._b)/lattice._a)
            qeInput.namelist('system').add('celldm(3)', float(lattice._c)/lattice._a)
            if lattice._ibrav < 14:
                qeInput.namelist('system').add('celldm(4)', lattice._cAB)
            else:
                qeInput.namelist('system').add('celldm(4)', lattice._cBC)
                qeInput.namelist('system').add('celldm(5)', lattice._cAC)
                qeInput.namelist('system').add('celldm(6)', lattice._cAB)
        else:
            if lattice._type == 'traditional':
                qeInput.namelist('system').add('ibrav', lattice._ibrav)
                qeInput.namelist('system').add('A', lattice._a)
                qeInput.namelist('system').add('B', lattice._b)
                qeInput.namelist('system').add('C', lattice._c)
                qeInput.namelist('system').add('cosAB', lattice._cAB)
                qeInput.namelist('system').add('cosAC', lattice._cAC)
                qeInput.namelist('system').add('cosBC', lattice._cBC)
            else:
                if 'generic' in lattice._type:
                    qeInput.namelist('system').add('celldm(1)', lattice._a)
                    lattice._ibrav = 0
                    qeInput.namelist('system').add('ibrav', lattice._ibrav)
                    if lattice._type == 'generic hexagonal':
                        cardArg = 'hexagonal'
                    if lattice._type == 'generic cubic' or lattice._type == None:
                        cardArg = 'cubic'
                    qeInput.createCard('cell_parameters')
                    qeInput.card('cell_parameters').setArg(cardArg)
                    qeInput.card('cell_parameters').removeLines()
                    for i in range(3):
                        v = lattice.diffpy().base[i,:]/float(lattice._a)
                        qeInput.card('cell_parameters').addLine(\
                                       lattice.formatString%(v[0], v[1], v[2]))
        
        #************************* updateStructure:*****************************
        
        if 'system' not in qeInput.namelists:
            qeInput.addNamelist('system')            
        qeInput.namelist('system').remove('ntyp')
        qeInput.namelist('system').remove('nat')
        if structure.ntyp != None:
            qeInput.namelist('system').add('ntyp', structure.ntyp)
        if structure.nat != None:
            qeInput.namelist('system').add('nat', structure.nat)
        
        if len(qeInput.namelist('system').params) == 0:
            qeInput.removeNamelist('system')  

        if 'atomic_positions' in qeInput.cards:
            qeInput.removeCard('atomic_positions')
        qeInput.createCard('atomic_positions')
        qeInput.card('atomic_positions').setArg(structure.atomicPositionsType)
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
            qeInput.card('atomic_positions').addLine(line)
        
        if len(qeInput.card('atomic_positions').lines()) == 0:
            qeInput.removeCard('atomic_positions')

        # update ATOMIC_SPECIES card
        if 'atomic_species' in qeInput.cards:
            qeInput.removeCard('atomic_species')
        qeInput.createCard('atomic_species')
        for element, specie in structure.atomicSpecies.items():
            qeInput.card('atomic_species').addLine(specie.toString())
        
        if len(qeInput.card('atomic_species').lines()) == 0:
            qeInput.removeCard('atomic_species')                


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
