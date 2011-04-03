#!/usr/bin/env python
"""class QEAtom for storing properties of a single atom"""
import numpy
from qelattice import QELattice

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

class QEAtom(object):
    """Atom --> class for storing atom information
    Data members:
        atype     -- type of the atom
        xyz         -- fractional coordinates
        name        -- atom label
        xyz_cartn   -- absolute Cartesian coordinates, property synced with xyz
        lattice     -- coordinate system for fractional coordinates,
                       an instance of Lattice or None for Cartesian system
    """

    def __init__(self, atype=None, xyz=None, mass = None,  \
                 potential = None, lattice=None, optConstraint = None, name=None):
        """Create atom of a specified type at given lattice coordinates.
        Atom(a) creates a copy of Atom instance a.
        atype         -- symbol string or Atom instance
        xyz           -- fractional(crystal) coordinates
        name          -- atom label
        mass          -- atom mass
        potential     -- pseudopotential file name
        lattice       -- QE coordinate system for fractional coordinates
        optConstraint -- list of up to three constraints for each coordinate 
                         for QE geometry optimization (0 or 1)
        """
        # declare data members
        self.symbol = None
        self._xyz = numpy.zeros(3, dtype=float)
        self.name = ''
        self._mass = 0
        self._potential = ''
        self.lattice = QELattice()
        from pwinput import PWInput
        self.lattice._qeInput = PWInput()
        self._optConstraint = numpy.array([], dtype = int)
        # assign them as needed
        if isinstance(atype, QEAtom):
            atype_dup = atype.__copy__()
            self.__dict__.update(atype_dup.__dict__)
        else:
            self.symbol = atype
        # take care of remaining arguments
        if xyz is not None:             self._xyz[:] = xyz
        if name is not None:            self.name = name
        if mass is not None:            self._mass = mass
        if potential is not None:       self._potential = potential
        if lattice is not None:         self.lattice = lattice
        if optConstraint is not None:  self._optConstraint = \
                                        numpy.array(optConstraint, dtype = int)
        return


    def __repr__(self):
        """simple string representation"""
        xyz = self.xyz
        s = "%-4s %8.6f %8.6f %8.6f %8.6f %s" % \
                (self.symbol, xyz[0], xyz[1], xyz[2], self.mass, self.potential)
        return s

    def __copy__(self):
        """Return a copy of this instance.
        """
        adup = QEAtom(self.symbol)
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

    
    def _get_xyz(self):
        return self._xyz
    
    def _set_xyz(self, value):
        self._xyz = value
        self.lattice._qeInput.update()

    xyz = property(_get_xyz, _set_xyz, doc =
        """fractional coordinates of an atom
        """ )
    
    
    def _get_mass(self):
        return self._mass
    
    def _set_mass(self, value):
        self._mass = value
        self.lattice._qeInput.update()
        
    mass = property(_get_mass, _set_mass, doc =
        """mass of an atom """ ) 
    

    def _get_potential(self):
        return self._potential
    
    def _set_potential(self, value):
        self._potential = value
        self.lattice._qeInput.update()
        
    potential = property(_get_potential, _set_potential, doc =
        """property. Name of paseudopotential of an atom """)
    
    
    def _get_optConstraint(self):
        return self._optConstraint
    
    def _set_optConstraint(self, value):
        self._optConstraints = value
        self.lattice._qeInput.update()
        
    optConstraint = property(_get_optConstraint, _set_optConstraint, doc =
        """optimization constraint, e.g. [1, 0, 1] of an atom for QE geometry  
optimization""")
# End of class QEAtom
