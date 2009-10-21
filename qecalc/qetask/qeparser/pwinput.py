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

from qeinput import QEInput
from qestructure import QEStructure

class PWInput(QEInput):
    def __init__(self, filename=None, config=None):
        QEInput.__init__(self,filename, config, type='pw')
        self.structure = None
    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary
            Initializes structure as well"""
        (self.namelists, self.cards) = self.parser.parse()
        self.structure = QEStructure(self)

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
