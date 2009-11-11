#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Configuration parameters from INPUT_PWCOND.html file
# First value is default value

# Namelist: INPUTCOND

namelists = ('inputcond',)

cards     = ('k_and_energy_points',)

# Namelist: INPUTCOND
namelist_inputcond = ('outdir',
                      'prefixt',
                      'prefixl',
                      'prefixs',
                      'prefixr',
                      'band_file',
                      'tran_file',
                      'save_file',
                      'fil_loc',
                      'lwrite_cond',
                      'lread_cond',
                      'lwrite_loc',
                      'lread_loc',
                      'ikind',
                      'iofspin',
                      'llocal',
                      'bdl',
                      'bds',
                      'bdr',
                      'nz1',
                      'energy0',
                      'denergy',
                      'nenergy',
                      'ecut2d',
                      'ewind',
                      'epsproj',
                      'orbj_in',
                      'orbj_fin ')

# Card: K_and_Energy_Points
card_k_and_energy_points = ('nkpts',
                            'kx',
                            'ky',
                            'weight',
                            'nenergy')

__date__ = "$Sep 2, 2009 12:40:06 PM$"


