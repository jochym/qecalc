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



class QEStructureError(Exception):
    """Exception for inappropriate structure operation
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value
    
    
class QELatticeError(Exception):
    """Exception for inappropriate lattice operation
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value    