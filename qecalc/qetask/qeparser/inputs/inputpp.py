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

# Configuration parameters from INPUT_PP.html file
# First value is default value

# Namelist: INPUTPP

namelists   = ('inputpp',
               'plot')
cards       = ()

namelist_inputpp = ('prefix',
                    'outdir',
                    'filplot',
                    'plot_num',
                    'spin_component',
                    'spin_component',
                    'sample_bias',
                    'stm_wfc_matching',
                    'z',
                    'dz',
                    'kpoint',
                    'kband',
                    'lsign',
                    'spin_component',
                    'emin',
                    'emax',
                    'spin_component',
                    'spin_component')

namelist_plot = ('nfile',
                 'filepp',
                 'weight',
                 'iflag',
                 'output_format',
                 'fileout',
                 'e1',
                 'x0',
                 'nx',
                 'e1',
                 'e2',
                 'x0',
                 'nx',
                 'ny',
                 'e1',
                 'e2',
                 'e3',
                 'x0',
                 'nx',
                 'ny',
                 'nz',
                 'radius',
                 'nx',
                 'ny')


__date__ = "$Sep 2, 2009 11:15:18 AM$"


