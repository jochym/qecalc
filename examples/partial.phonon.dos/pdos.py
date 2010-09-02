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

from qecalc.qetask.matdyntask import MatdynTask
from qecalc.qetask.pwtask import PWTask
from qeutils.phdos import PhononDOS

if __name__ == "__main__":
    matdyn = MatdynTask(configString = "")
    pw = PWTask(configString = "")
    dos = PhononDOS(matdynTask = matdyn, structure = pw.input.structure)

    # grid has to be huge since a simple histogramming will be used
    qgrid = [64, 64, 64]

    # generate qpoints and launch matdyn task
    dos.launch(nqpoints = qgrid, partialDOS = True )

    #dos.DOS()
    # partDOS will calculate DOS for each atomic site. It will also generate
    # total DOS and DOS by element type
    dos.partDOS()

    # save all available DOSes to disk
    dos.save()