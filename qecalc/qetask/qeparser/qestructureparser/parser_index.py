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

parser_index = {

    # automatic format detection - tries all parsers one by one
#    'auto' : {
#        'module' : 'P_auto',
#        'file_extension' : '',
#        'file_pattern' : '*.*',
#        },

    # PW input format format
    'pwinput' : {
        'module' : 'P_pwinput',
        'file_extension' : '.in',
        'file_pattern' : '*.in',
        },

    # PW output format
    'pwoutput' : {
        'module' : 'P_pwoutput',
        'file_extension' : '.out',
        'file_pattern' : '*.out',
        },        
}