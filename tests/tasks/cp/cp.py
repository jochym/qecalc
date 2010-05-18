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
from qecalc.qetask.cptask import CPTask
configString = """
# all the relevant input files must be preconfigured for specific tasks
# before using this class

[Launcher]
# parallelization parameters
# if this section is empty - serial mode is used
paraPrefix:   mpiexec -n 8

outdir: temp/

[cp.x]
# pw input/output files
cpInput:  cp.in
cpOutput: cp.out

"""

if __name__ == "__main__":
    cp = CPTask(configString = configString)
    cp.launch()

__author__="Nikolay Markovskiy"
__date__ ="$May 17, 2010 4:21:30 PM$"
