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

# XXX: Borrowed from inputpw.py. Modify fields according to the documentation

namelists = ('control',
             'system',
             'electrons',
             'ions',
             'cell',
             'phonon')

cards =   ('atomic_species',
           'atomic_positions',
           'k_points',
           'cell_parameters',
           'climbing_images',
           'constraints',
           'collective_vars',
           'occupations')

# Namelist: CONTROL
namelist_control = ('calculation',  #: ('scf', 'nscf', 'bands', 'phonon', 'relax', 'md', 'vc-relax', 'vc-md', 'neb', 'smd', 'metadyn'),
                    'title',        #: (''),
                    'verbosity',    #: ('high', 'default', 'low', 'minimal'),
                    'restart_mode', #: ('from_scratch', 'restart'),
                    'wf_collect',   #: ('.false.', '.true.'),
                    'nstep',        #: ('???'),
                    'iprint',
                    'tstress',
                    'tprnfor',
                    'dt',
                    'outdir',
                    'wfcdir',
                    'prefix',
                    'lkpoint_dir',
                    'max_seconds',
                    'etot_conv_thr',
                    'forc_conv_thr',
                    'disk_io',
                    'pseudo_dir',
                    'tefield',
                    'dipfield',
                    'lelfield',
                    'lberry',
                    'gdir',
                    'nppstr',
                    'nberrycyc')

system_namelist = ('ibrav',
                   'celldm',
                   'A',
                   'B',
                   'C',
                   'cosAB',
                   'cosAC',
                   'cosBC',
                   'nat',
                   'ntyp',
                   'nbnd',
                   'nelec',
                   'tot_charge',
                   'ecutwfc',
                   'ecutrho',
                   'nr1',
                   'nr2',
                   'nr3',
                   'nr1s',
                   'nr2s',
                   'nr3s',
                   'nosym',
                   'noinv',
                   'occupations',
                   'degauss',
                   'smearing',
                   'nspin',
                   'noncolin',
                   'starting_magnetization',
                   'nelup',
                   'neldw',
                   'multiplicity',
                   'tot_magnetization',
                   'ecfixed',
                   'qcutz',
                   'q2sigma',
                   'xc_type',
                   'lda_plus_u',
                   'Hubbard_alpha',
                   'Hubbard_U',
                   'starting_ns_eigenvalue(m,ispin,I)',
                   'U_projection_type',
                   'edir',
                   'emaxpos',
                   'eopreg',
                   'eamp',
                   'angle1',
                   'angle2',
                   'constrained_magnetization',
                   'fixed_magnetization',
                   'B_field',
                   'lambda',
                   'report',
                   'lspinorb',
                   'assume_isolated')

namelist_electrons = ('electron_maxstep',
                      'conv_thr',
                      'mixing_mode',
                      'mixing_beta',
                      'mixing_ndim',
                      'mixing_fixed_ns',
                      'diagonalization',
                      'ortho_para',
                      'diago_thr_init',
                      'diago_cg_maxiter',
                      'diago_david_ndim',
                      'diago_full_acc',
                      'efield',
                      'startingpot',
                      'startingwfc',
                      'tqr')

namelist_ions = ('ion_dynamics',
                 'phase_space',
                 'pot_extrapolation',
                 'wfc_extrapolation',
                 'remove_rigid_rot',
                 'ion_temperature',
                 'tempw',
                 'tolp',
                 'delta_t',
                 'nraise',
                 'refold_pos',
                 'upscale',
                 'bfgs_ndim',
                 'trust_radius_max',
                 'trust_radius_min',
                 'trust_radius_ini',
                 'w_1',
                 'w_2',
                 'num_of_images',
                 'opt_scheme',
                 'CI_scheme',
                 'first_last_opt',
                 'temp_req',
                 'ds',
                 'k_max',
                 'k_min',
                 'path_thr',
                 'use_masses',
                 'use_freezing',
                 'fe_step',
                 'g_amplitude',
                 'fe_nstep',
                 'sw_nstep')

namelist_cell = ('cell_dynamics | press | wmass | cell_factor | press_conv_thr | cell_dofree',)

namelist_phonon = ('modenum | xqq',)

card_atomic_species = ('atomic_species',)        # X | Mass_X | PseudoPot_X

card_atomic_positions = ('atomic_positions',)    # X | x | y | z | if_pos(1) | if_pos(2) | if_pos(3) |

card_k_points = ('k_points',)                    # nks | xk_x | xk_y | xk_z | wk | nk1 | nk2 | nk3 | sk1 | sk2 | sk3

card_cell_parameters = ('cell_parameters',)      # v1 | v2 | v3

card_climbing_images = ('climbing_images',)      # index1, index2, ... indexN

card_constraints = ('constraints',)              # nconstr | constr_tol | constr_type | constr(1) | constr(2) | constr(3) | constr(4) | constr_target

card_collective_vars = ('collective_vars',)      # ncolvar | colvar_tol | colvar_type | colvar(1) | colvar(2) | colvar(3) | colvar(4)

card_occupations = ('occupations',)              # f_inp1 | f_inp2


__date__ = "$May 17, 2010 4:04:06 PM$"


