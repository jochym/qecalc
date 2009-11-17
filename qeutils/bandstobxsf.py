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
# To change this template, choose Tools | Templates
# and open the template in the editor.

import numpy

def bandstobxsf(energy, recipLatticeVectors, bands, fname):

    file = open(fname,'w')

    kpoints = [bands.shape[1],bands.shape[2],bands.shape[3]]
    s =  """ BEGIN_INFO
   #
   # this is a Band-XCRYSDEN-Structure-File
   # aimed at Visualization of Fermi Surface
   #"""
    s = s + '\n' + '   Fermi Energy:        ' + str(energy) + '\n'
    s = s + """ END_INFO
 BEGIN_BLOCK_BANDGRID_3D
 band_energies
 BANDGRID_3D_BANDS\n"""
    s = s + '%5u\n'%bands.shape[0] + \
             '%5u %5u %5u\n'%(kpoints[0], kpoints[1],kpoints[2]) + \
             '  %#.8f %#.8f %#.8f\n'%(0.0, 0.0, 0.0)
    for i in range(recipLatticeVectors.shape[0]):
        s = s + '  %#.8f %#.8f %#.8f\n'%(recipLatticeVectors[i,0], \
                            recipLatticeVectors[i,1], recipLatticeVectors[i,2])
    for iband in range(bands.shape[0]):
        s = s + 'BAND:%5u\n'%(iband+1)

        for kz in range(kpoints[2]):
            for ky in range(kpoints[1]):
                for kx in range(kpoints[0]):
                    s = s + '%#.4f '%(bands[iband, kx, ky, kz])
                s = s + '\n'
            s = s + '\n'

    s = s + ' END_BANDGRID_3D\n END_BLOCK_BANDGRID_3D'
    file.write(s)
    file.close()



if __name__ == "__main__":
    from qecalc.qetask.pwtask import PWTask
    from qeutils import kmesh
    pw = PWTask('config.ini')
    pw.input.parse()
    pw.output.parse()
    #print pw.output.property('fermi energy')
    kpoints = [17,17,17]

    bands = kmesh.packBands(kpoints,pw.output.property('bands')[1])
#    numpy.zeros((pw.output.property('bands')[1].shape[1], kpoints[0], kpoints[1], kpoints[2]))
#    print bands
#    for iband in range(pw.output.property('bands')[1].shape[1]):
#        for kz in range(kpoints[2]):
#            for ky in range(kpoints[1]):
#                for kx in range(kpoints[0]):
#                    bands[iband,kx,ky,kz] = \
#                      pw.output.property('bands')[1][\
#                      kpoints[2]*kpoints[1]*kz + kpoints[1]*ky + kx,iband]
    #bands = pw.output.property('bands')[1].transpose().reshape(\
    #                             -1, kpoints[0],kpoints[1],kpoints[2])
    print bands
    print bands.shape
    #pw.output.property('fermi energy')[0]
    bandstobxsf(7.8071, \
                pw.input.structure.lattice.reciprocalBase(), \
                bands, 'test.bxsf')
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Nov 15, 2009 5:12:00 PM$"
