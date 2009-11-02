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

# Configuration parameters from INPUT_GIPAW.html file
# First value is default value

# Namelist: INPUTGIPAW

namelists   = ('inputgipaw',)
cards       = ()

namelist_inputgipaw = ('job',
                    'prefix',
                    'tmp_dir',
                    'conv_threshold',
                    'isolve',
                    'q_gipaw',
                    'iverbosity',
                    'filcurr',
                    'filfield',
                    'read_recon_in_paratec_fmt',
                    'file_reconstruction',
                    'use_nmr_macroscopic_shape',
                    'nmr_macroscopic_shape(3,3)',
                    'spline_ps')


__date__ = "$Sep 2, 2009 11:51:30 AM$"


