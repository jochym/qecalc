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

# Configuration parameters from INPUT_PH.html file
# First value is default value

# Namelist: INPUTPH

namelists   = ('inputph',)
cards       = ()

namelist_inputph = ('amass',
                    'outdir',
                    'prefix',
                    'niter_ph',
                    'tr2_ph',
                    'alpha_mix(niter)',
                    'nmix_ph',
                    'iverbosity',
                    'reduce_io',
                    'max_seconds',
                    'fildyn',
                    'fildrho',
                    'fildvscf',
                    'epsil',
                    'lrpa',
                    'lnoloc',
                    'trans',
                    'lraman',
                    'eth_rps',
                    'eth_ns',
                    'dek',
                    'recover',
                    'elph',
                    'zue',
                    'elop',
                    'fpol',
                    'lnscf',
                    'ldisp',
                    'nq1',
                    'nq2',
                    'nq3',
                    'iq1',
                    'iq2',
                    'iq3',
                    'nrapp',
                    'maxirr',
                    'nat_todo ')

#Line-of-input:  xq(1) xq(2) xq(3)

#Line-of-input: irrep(1) irrep(2) ... irrep(nrapp)

#Line-of-input: atom(1) atom(2) ... atom(nat_todo)

__date__ = "$Sep 2, 2009 10:49:46 AM$"


