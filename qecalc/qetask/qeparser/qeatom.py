#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File refactored  by:      Nikolay Markovskiy
#      based on Atom.py from diffpyStructure package coded by Pavol Juhas
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""class QEAtom for storing properties of a single atom"""
import numpy

class CartesianCoordinatesArray(numpy.ndarray):
    """Helper array for accessing Cartesian coordinates.
    Converts and updates related array of corresponding fractional
    coordinates.

    Data members:
        lattice -- instance of Lattice defining fractional coordinates
        xyz     -- instance of numpy.array storing fractional coordinates
    """

    def __new__(self, lattice, xyz):
        return numpy.zeros(3, dtype=float).view(self)

    def __init__(self, lattice, xyz):
        self.lattice = lattice
        self.xyz = xyz
        self[:] = self.lattice.cartesian(self.xyz)
        pass

    def __setitem__(self, idx, value):
        """Set idx-th coordinate and update linked self.xyz

        idx     -- index in xyz array
        value   -- new value of x, y or z
        """
        numpy.ndarray.__setitem__(self, idx, value)
        self.xyz[:] = self.lattice.fractional(self)
        return

# End of CartesianCoordinatesArray


class QEAtom(object):
    """Atom --> class for storing atom information

    Data members:
        element     -- type of the atom
        xyz         -- fractional coordinates
        name        -- atom label
        xyz_cartn   -- absolute Cartesian coordinates, property synced with xyz
        lattice     -- coordinate system for fractional coordinates,
                       an instance of Lattice or None for Cartesian system

    Private data:
    """


    def __init__(self, atype=None, xyz=None, name=None, mass = None,  \
                 potential = None, lattice=None):
        """Create atom of a specified type at given lattice coordinates.
        Atom(a) creates a copy of Atom instance a.

        atype       -- element symbol string or Atom instance
        xyz         -- fractional(crystal) coordinates
        name        -- atom label
        mass        -- atom mass
        potential   -- pseudopotential file name
        lattice     -- QE coordinate system for fractional coordinates
        """
        # declare data members
        self.element = None
        self.xyz = numpy.zeros(3, dtype=float)
        self.name = ''
        self.mass = None
        self.potential = None
        self.lattice = None
        # assign them as needed
        if isinstance(atype, QEAtom):
            atype_dup = atype.__copy__()
            self.__dict__.update(atype_dup.__dict__)
        else:
            self.element = atype
        # take care of remaining arguments
        if xyz is not None:         self.xyz[:] = xyz
        if name is not None:        self.name = name
        if mass is not None:        self.mass = mass
        if potential is not None:   self.potential = potential
        if lattice is not None:     self.lattice = lattice
        return


    def __repr__(self):
        """simple string representation"""
        xyz = self.xyz
        s = "%-4s %8.6f %8.6f %8.6f %8.6f %s" % \
                (self.element, xyz[0], xyz[1], xyz[2], self.mass, self.potential)
        return s

    def __copy__(self):
        """Return a copy of this instance.
        """
        adup = QEAtom(self.element)
        adup.__dict__.update(self.__dict__)
        # create copies for what should be copied
        adup.xyz = numpy.array(self.xyz)
        return adup

    ####################################################################
    # property handlers
    ####################################################################

    # xyz_cartn

    def _get_xyz_cartn(self):
        if not self.lattice:
            rv = self.xyz
        else:
            rv = CartesianCoordinatesArray(self.lattice, self.xyz)
        return rv

    def _set_xyz_cartn(self, value):
        if not self.lattice:
            self.xyz[:] = value
        else:
            self.xyz = self.lattice.fractional(value)
            self.lattice._qeInput.update()
        return

    xyz_cartn = property(_get_xyz_cartn, _set_xyz_cartn, doc =
        """absolute Cartesian coordinates of an atom
        """ )

# End of class QEAtom
