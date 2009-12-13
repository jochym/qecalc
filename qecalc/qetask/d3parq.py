#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from mpi4py import MPI
from d3task import D3Task
import sys

def D3ParQ(qpoints, d3):
    outDir = d3.setting.outDir
    myrank = MPI.COMM_WORLD.Get_rank()



if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        raise
    
    configName = sys.argv[1]
    d3 = D3Task(configName)
    outDir = d3.setting.outDir
    myrank = MPI.COMM_WORLD.Get_rank()
    outDir = outDir + '/d3task_' + str(myrank)

    d3.setting.outdir = outDir

__author__="Nikolay Markovskiy"
__date__ ="$Dec 11, 2009 2:20:08 PM$"
