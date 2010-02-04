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

from qeutils.converger import Converger
configString = """
# all the relevant input files must be preconfigured for specific tasks
# before using this class

[Launcher]
# parallelization parameters
# if this section is empty - serial mode is used
paraPrefix:   mpiexec -n 8
paraPostfix: -npool 8


outdir: temp/


[pw.x]
# pw input/output files
pwfInput:  scf.in
pwOutput: scf.out

[Converger]
# taskName can be 'total energy' or 'single phonon'
taskName:   total energy
tolerance:  0.1
what:       conv_thr
startValue: 1e-6
#step      : 5
multiply: 0.1
"""

if __name__ == "__main__":
    task = Converger(configString = configString)
    opt_value = task.converge()
