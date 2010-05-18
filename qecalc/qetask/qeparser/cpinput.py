#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
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

class CPInput(QEInput):
    def __init__(self, filename=None, config=None):
        QEInput.__init__(self,filename, config, type='cp')
        self.structure = QEStructure(self)
        self.kpoints = PWKpoints(self)
    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary
            Initializes structure as well"""
        (self.header, self.namelists, self.cards, self.attach) = self.parser.parse()
        self.structure.parseInput()
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

    def outDir(self):
        self.parse()
        return self.namelist('control').param('outdir', quotes = False)

textA = """
  &control
    title = 'Ammonia',
    calculation = 'fpmd',
    restart_mode = 'from_scratch', ! 'restart',
    ndr = 51,
    ndw = 51,
    nstep  = 100,
    iprint = 10,
    isave  = 100,
    tstress = .TRUE.,
    tprnfor = .TRUE.,
    dt    = 5.0d0,
    etot_conv_thr = 1.d-9,
    ekin_conv_thr = 1.d-4,
    prefix = 'nh3_mol'
    pseudo_dir='./',
    outdir='temp/',
 /
 &system
    ibrav = 14,
    celldm(1) = 12.0,
    celldm(2) = 1.0,
    celldm(3) = 1.0,
    celldm(4) = 0.0,
    celldm(5) = 0.0,
    celldm(6) = 0.0,
    nat  = 4,
    ntyp = 2,
    nbnd = 4,
    nelec = 8,
    ecutwfc = 80.0,
!    ecfixed = 68.0,
!    qcutz = 68.0,
!    q2sigma = 8.0,
    input_dft = 'BLYP'
 /
 &electrons
    emass = 400.d0,
    emass_cutoff = 2.5d0,
    orthogonalization = 'ortho',
    ortho_eps = 5.d-8,
    ortho_max = 15,
    electron_dynamics = 'sd',
!    electron_damping = 0.3,
    electron_velocities = 'zero',
    electron_temperature = 'not_controlled',
 /
 &ions
    ion_dynamics = 'damp',
    ion_damping = 0.2,
    ion_radius(1) = 0.8d0,
    ion_radius(2) = 0.8d0,
    ion_velocities = 'zero',
    ion_temperature = 'not_controlled',
    ion_nstepe = 10
 /
 &cell
    cell_dynamics = 'none',
    cell_velocities = 'zero',
    press = 0.0d0,
 /
ATOMIC_SPECIES
 N 16.0 N.BLYP.UPF 4
 H  1.0 H.fpmd.UPF 4
ATOMIC_POSITIONS (bohr)
   N     0.000825    0.000825   0.0000
   H     0.159883333   -0.020358333   -0.0184
   H    -0.019208333    0.159883333   -0.017866667
   H    -0.014958333   -0.015058333   0.159883333
"""

if __name__ == "__main__":

    inp = CPInput(config=textA)
    inp.parse()
    print inp.toString()

__author__="Nikolay Markovskiy"
__date__ ="$May 17, 2010 4:06:15 PM$"
