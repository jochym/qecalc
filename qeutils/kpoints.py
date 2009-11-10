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
import numpy
class Kpoints:
    def _kMesh(self, nkx, nky, nkz):
        """Will generate k-point mesh in crystal coordinates in [0,1] box
        Mesh is uniform. I.e. lattice symmetry is not taken into account."""
        j = numpy.complex(0,1)
        data = numpy.mgrid[0.0:1.0:nkx*j,0:1.0:nky*j,0:1.0:nkz*j]
        data = numpy.transpose(data.reshape(3,nkx*nky*nkz))
        return data

    def kMeshCart(self, nq1, nq2, nq3):
        """Will generate k-point mesh in cartesian coordinates
        Lattice symmetry is not taken into account."""
        kpts = self._kMesh(nq1, nq2, nq3)
        kpts_cart = numpy.zeros(kpts.shape)
        for k in range(kpts.shape[0]):
            # convert into cartesian coordinates in units of inv lattice vector a
            kpts_cart[k,:] = self.structure.lattice.recipCartesian(kpts[k,:])
#            self.structure.lattice.diffpy().cartesian(kpt)/ \
#                        self.structure.lattice.a

        return kpts_cart

    def recipCartesian(self, kPoint):
        """Conversts vector on fractional coordinates in reciprocal space into
           a vector in cartesian coordinates"""
        recip_base = self.diffpy().reciprocal().base*self._a
        return numpy.dot( kPoint, recip_base)