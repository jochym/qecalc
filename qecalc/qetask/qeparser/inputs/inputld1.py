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

# Configuration parameters from INPUT_LD1.html file
# First value is default value

namelists = ('input',
             'inputp',
             'test')
cards       = ()

cards     = ('allelectroncards',
             'pseudopotentialgenerationcards',
             'pseudopotentialtestcards')

# Namelist: INPUT
namelist_input = ('title',
                  'zed',
                  'atom',
                  'xmin',
                  'dx',
                  'rmax',
                  'beta',
                  'tr2',
                  'iswitch',
                  'nld',
                  'rlderiv',
                  'eminld',
                  'emaxld',
                  'deld',
                  'rel',
                  'lsmall',
                  'lsd',
                  'dft',
                  'latt',
                  'isic',
                  'rytoev_fact',
                  'cau_fact',
                  'vdw',
                  'prefix',
                  'verbosity',
                  'config',
                  'rel_dist',
                  'write_coulomb')

# Namelist: INPUTP
namelist_inputp = ('zval',
                   'pseudotype',
                   'file_pseudopw',
                   'file_recon',
                   'lloc',
                   'rcloc',
                   'nlcc',
                   'new_core_ps',
                   'rcore',
                   'tm',
                   'rho0',
                   'lpaw',
                   'which_augfun',
                   'rmatch_augfun',
                   'lsave_wfc',
                   'author',
                   'file_chi',
                   'file_beta',
                   'file_qvan',
                   'file_screen',
                   'file_core',
                   'file_wfcaegen',
                   'file_wfcncgen',
                   'file_wfcusgen')

# Namelist: TEST
namelist_test = ('nconf',
                 'file_pseudo',
                 'ecutmin',
                 'ecutmax',
                 'decut',
                 'rm',
                 'configts',
                 'lsdts',
                 'frozen_core')

# Card: AllElectronCards
card_allelectroncards = ('nwf',
                         'nl',
                         'n',
                         'l',
                         'oc',
                         'isw',
                         'jj')

# Card: PseudoPotentialGenerationCards
card_pseudopotentialgenerationcards = ('nwfs',
                                       'nls',
                                       'nns',
                                       'lls',
                                       'ocs',
                                       'ener',
                                       'rcut',
                                       'rcutus',
                                       'jjs')

# Card: PseudoPotentialTestCards
card_pseudopotentialtestcards = ('nwfts',
                                 'elts',
                                 'nnts',
                                 'llts',
                                 'octs',
                                 'enerts',
                                 'rcutts',
                                 'rcutusts',
                                 'iswts',
                                 'jjts')


__date__ = "$Sep 2, 2009 12:02:23 PM$"


