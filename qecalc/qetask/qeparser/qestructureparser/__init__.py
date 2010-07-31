#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


try:
    from diffpy.Structure.structure import Structure
    from diffpy.Structure.atom import Atom
    from diffpy.Structure.lattice  import cosd, Lattice
except ImportError:
    from matter import Structure, Atom, Lattice
    from matter.Lattice import cosd

from qecalc.qetask.qeparser.qestructure import  QEStructure, AtomicSpecies
from qecalc.qetask.qeparser.qestructure import QELattice
from qecalc.qetask.qeparser.orderedDict import OrderedDict

import numpy

from qecalc.qetask.qeparser.ParserErrors import *

from qecalc.qetask.qeparser.qestructureparser.parser_index import parser_index  

from qecalc.qetask.qeparser.qestructureparser.qestructureparser \
                                                        import QEStructureParser