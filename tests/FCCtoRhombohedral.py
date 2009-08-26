#!/usr/bin/env python

"""Express Ni FCC structure in its rhombohedral primitive cell.
"""

from diffpy.Structure import Structure
from diffpy.Structure import Lattice
from diffpy.Structure.SymmetryUtilities import equalPositions

# load Ni fcc structure, Ni.stru is in tests/testdata/
nifcc = Structure(filename='Ni.stru')
# base vectors for rhombohedral lattice in Cartesian coordinates
rhombbase = nifcc.lattice.cartesian(
    [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0, 0.5]])

# make a deep copy of nifcc
nirhomb = Structure(nifcc)
nirhomb.placeInLattice(Lattice(base=rhombbase))

# collect atoms that are at equivalent position to some previous atom
duplicates = set([a1
    for i0, a0 in enumerate(nirhomb) for a1 in nirhomb[i0+1:]
        if equalPositions(a0.xyz, a1.xyz, eps=1e-4)])

# Filter out duplicate atoms.  Use slice assignment so that
# nirhomb is not replaced with a list.
nirhomb[:] = [a for a in nirhomb if not a in duplicates]

# use rstrip to avoid duplicate line feed
print nirhomb.writeStr(format='discus').rstrip()
