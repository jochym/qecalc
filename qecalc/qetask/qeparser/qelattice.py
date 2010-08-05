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
try:
    from diffpy.Structure import Structure
    from diffpy.Structure.lattice import Lattice, cosd
except ImportError:
    from matter import Structure, Lattice
    from matter.Lattice import cosd


from math import sqrt, degrees
import numpy
from qeinput import QEInput


class QELattice(object):
    """Class QELattice for working with crystal lattices in QE notation
       Uses diffpy.Lattice class for storage
       
       Following properties are dynamically linked to other properties  
       (E.g. lattice vectors) and QEInput(if QEInput.autoUpdate = True(default)):
       ibrav                  -- lattice type, setting it into a different value 
                                 will automatically update lattice vectors, 
                                 QE parsing object and structure
                                 if ibrav = 0, only 'a' parameter is relevant
       a, b, c, cBC, cAC ,cAB -- lattice parameters, setting any of them will 
                                 dynamically update the QE parsing object
                                 (e.g. pw.input) and structure(if relevant). 
                                 They are also ibrav sensitive.  E.g. if 
                                 ibrav = 1 (simple cubic). Setting 'a' to a \
                                 different value will also modify  b and c. 
                                  
       type                   -- lattice type to save into PWSCF cfg file 
                                 'celldm'  (default)
                                 'traditional' - using a,b,c,cosAC, cosAB, cosBC;
                                 'generic cubic', 'generic hexagonal' - assume 
                                  that section 'CELL_PARAMETERS' exists
                                 'generic' types also  assume/set ibrav = 0
      setLattice()            -- will set everything at once
    """


    def __init__(self, ibrav = 1,a = 1. ,b = 1.,c = 1.,
                 cBC = 0.,cAC = 0. ,cAB = 0., base = None ):
        self.formatString = '%# .8f %# .8f %# .8f'
        self._qeInput = None
        self._type = 'celldm'
        
        # diffpyStructure container class, used for lattice operations
        self.__primitiveLattice = Lattice()
        
        # Lattice vectors in bohr or angstrom:
        self._base = None

        # initialize the lattice if there is enough information 
        if ibrav > 0 and base != None:
            self.setLatticeFromQEVectors(ibrav, base)
        else:
            self.setLattice(ibrav ,a ,b , c, cBC ,cAC ,cAB, base)


    def __str__(self):
        """simple string representation"""
        st = ''
        if self._ibrav == 0:
            st = st + '"generic" cell:\n'
        else:
            qeBaseTuple = self._getQEBaseFromParCos(self._ibrav, self._a, self._b, \
                                       self._c, self._cBC, self._cAC, self._cAB)
            qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]
            st = st + '"' + qeBaseTuple[2] + '" cell:\n'
        
        for i in range(3):
            v = self.__primitiveLattice.base[i,:]
            st = st + self.formatString%(v[0], v[1], v[2])
            st = st + '\n'
      
        return st    
                
            
    def cartesian(self, u):
        """return cartesian coordinates of a lattice vector"""
        return self.__primitiveLattice.cartesian(u)
    
    def fractional(self, rc):
        """return fractional coordinates of a cartesian vector"""
        return self.__primitiveLattice.fractional(u)

    def dot(self, u, v):
        """return dot product of 2 lattice vectors"""
        return self.__primitiveLattice.dot(u, v)

    def norm(self, u):
        """return norm of a lattice vector"""
        return self.__primitiveLattice.norm(u)

    def dist(self, u, v):
        """Return distance of 2 points in lattice coordinates.
        """
        return __primitiveLattice.dist(u, v)

    def angle(self, u, v):
        """Return angle(u, v) --> angle of 2 lattice vectors in degrees.
        """
        return __primitiveLattice.angle(u, v)

    def recipCartesian(self, kPoint):
        """Conversts a vector in fractional coordinates in reciprocal space into
           a vector in cartesian coordinates"""
        recip_base = self.diffpy().reciprocal().base*self._a
        return numpy.dot( kPoint, recip_base)

    def reciprocalBase(self):
        """
        Get reciprocal lattice vectors in units of 2*pi/a
        """
        return self.diffpy().reciprocal().base*self._a


    def latticeParams(self):
        """Returns tuple of six lattice parameters:
            a, b, c, cos(BC), cos(AC), cos(AB)
        """
        return [self._a, self._b,self._c, self._cBC, self._cAC, self._cAB]        
        
        
    def setLattice(self, ibrav, a = None, b = None, c = None,
                   cBC = None, cAC = None, cAB = None, base = None, updateInput = True):
        """ 'base', numpy array of lattice vectors, and 'a'  will only be used
            if ibrav == 0. Otherwise, ibrav + lattice parameters will be used"""
        from math import acos
#        if [ibrav, a,b,c,cBC,cAC,cAB, base] == 8*[None]:
#            return None
        if ibrav == None:
            raise NonImplementedError('ibrav should be specified')
        self._ibrav = ibrav
        if self._ibrav == 0:
#            print 'Found "generic" cell:'
            if base == None:
                raise NonImplementedError('base must be specified')
            if a == None:
				raise # a and the base must be provided for ibrav = 0
            qeBase = numpy.array(base, dtype = float)#*self._a/a
            self._a = a
#            print qeBase
            #self._a = 1.0
            if 'generic' not in self._type:
                self._type = 'generic cubic'
            self.__primitiveLattice.setLatBase(qeBase)
            #self._standardLattice.setLatBase(qeBase)
        else:
            # Make sure all lattice parameters are mutually consistent 
            # according to ibrav. base array is not used:
            if ibrav < 4 or ibrav == 5: 
                if a is not None: 
                    self._a = self._b = self._c = a
                else:
                    if b is not None: 
                        self._a = self._b = self._c = b
                    else:
                        if c is not None: 
                            self._a = self._b = self._c = c                        
                        else:
                            self._b = self._c = self._a
                if ibrav < 4:                            
                    self._cBC = self._cAC = self._cAB = 0.0
            if ibrav == 4 or ibrav == 6 or ibrav == 7: 
                if a is not None: 
                    self._a = self._b = a
                else:
                    if b is not None: 
                        self._a = self._b = b
                    else:
                        self._b = self._a
                if c is not None: 
                    self._c = c
                if ibrav == 4:
                    self._cAB = cosd(120.0)
                    self._cBC = self._cAC = 0.0
            if ibrav == 5:
                if cBC is not None:
                    self._cBC = self._cAC = self._cAB = cBC
                else:
                    if cAC is not None:
                        self._cBC = self._cAC = self._cAB = cAC
                    else:
                        if cAB is not None:
                            self._cBC = self._cAC = self._cAB = cAB
                        else:
                            self._cAC = self._cAB = self._cBC
            if  ibrav == 6 or ibrav == 7:
                self._cBC = self._cAC = self._cAB = 0.0
            if ibrav > 7 and ibrav <= 14:
                if a is not None: 
                    self._a = a
                if b is not None: 
                    self._b = b
                if c is not None: 
                    self._c = c
                if ibrav > 7 and ibrav < 12:
                    self._cBC = self._cAC = self._cAB = 0.0
            if ibrav == 12 or ibrav == 13:
                if cAB is not None:
                    self._cAB = cAB
                else:
                    if cBC is not None or cAC is not None:
                        raise Exception("Should specify cos(AB) only for" + \
                                         " ibrav = 12 or ibrav = 13" )
                self._cBC = self._cAC = 0.0
            if ibrav == 14:
                if cBC is not None:
                    self._cBC = cBC
                if cAC is not None:
                    self._cAC = cAC
                if cAB is not None:
                    self._cAB = cAB
            qeBaseTuple = self._getQEBaseFromParCos(self._ibrav, self._a, self._b,
                                               self._c, self._cBC, self._cAC, self._cAB)
            qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]            
            self.__primitiveLattice.setLatBase(qeBase)
            alpha = degrees(acos(self._cBC))
            beta = degrees(acos(self._cAC))
            gamma = degrees(acos(self._cAB))
            #self._standardLattice.setLatPar(self._a,self._b,self._c,alpha,beta,gamma)
        #self._a0 = a
        self._base = qeBase
        
        if self._qeInput != None and updateInput == True:
            #self.updatePWInput()
            self._qeInput.update()       


    def diffpy(self):
        '''Returns diffpy.Lattice object. Do not use it for reading  QE
        (standard cell) lattice parameters. Use latticeParams, or a, b, c , ...
        instead'''
        return self.__primitiveLattice


    def updatePWInput(self, qeInput = None):
        """
        Deprecated
        """
        self._qeInput.update()
    

    def setLatticeFromQEVectors(self, ibrav, vectors):
        """ Will extract conventional lattice parameters from primitive vectors.
            'vectors' - is a list with primitive vectors (in QE notation),
            including lattice parameters. For example from PWSCF output"""
        from numpy import dot
        # default values:
        a = b = c = 1.0
        cBC = cAC = cAB = 0.0
        v = numpy.array(vectors, dtype = float)
        if ibrav == 0:
            self.setLattice(ibrav = 0, a =  1.0, base = vectors )
            return
            #raise NotImplementedError
        # sc simple cubic:
        if ibrav == 1:
            a = v[0,0]

        if ibrav == 2:
            a = b = c = sqrt( 2.0*dot(v[0,:],v[0,:]) )

        if ibrav == 3:
            a = b = c = 2.0*sqrt( dot(v[0,:],v[0,:])/3.0)

        if ibrav == 4:
            a = b = sqrt( dot(v[0,:],v[0,:]))
            c = sqrt( dot(v[2,:],v[2,:]))
            cAB = cosd(120.0)

        if ibrav == 5:
            a = b = c = sqrt( dot(v[0,:],v[0,:]))
            cBC = cAC = cAB = dot(v[0,:],v[2,:])/a**2

        if ibrav == 6:
            a = b = sqrt( dot(v[0,:],v[0,:]))
            c = sqrt( dot(v[2,:],v[2,:]))

        if ibrav == 7:
            a = b = v[1,0] - v[2,0]
            c = v[1,2] + v[2,2]

        if ibrav == 8:
            a = v[0,0]
            b = v[1,1]
            c = v[2,2]

        if ibrav == 9:
            a = v[0,0] - v[1,0]
            b = v[0,1] + v[1,1]
            c = v[2,2]

        if ibrav == 10:
            a = v[2,0] - v[0,0] - v[1,0]
            b = v[2,1] - v[0,1] + v[1,1]
            c = v[0,2] - v[1,2] + v[2,2]

        if ibrav == 11:
            a = v[0,0] - v[1,0]
            b = v[1,1] - v[2,1]
            c = v[0,2] - v[2,2]

        if ibrav == 12:
            a = v[0,0]
            b = sqrt( dot(v[1,:],v[1,:]))
            cAB = v[1,0]/b
            c = v[2,2]

        if ibrav == 13:
            a = v[0,0] + v[2,0]
            b = sqrt( dot(v[1,:],v[1,:]))
            c = v[2,2] - v[0,2]
            cAB = v[1,0]/b

        if ibrav == 14:
            a = v[0,0]
            b = sqrt( dot(v[1,:],v[1,:]))
            cAB = v[1,0]/b
            c = sqrt( dot(v[2,:],v[2,:]))
            cAC = v[2,0]/c
            cBC = v[2,1]*sqrt(1.0 - cAB**2)/c + cAC*cAB        
        self.setLattice(ibrav, a, b, c, cBC, cAC, cAB)


    ####################################################################
    # property handlers
    ####################################################################

    # lattice parameters


 #   def _get_a0(self):
 #       if self._a0 != None:
 #           return self._a0
 #       else:
 #           return self._a
 #   a0 = property(_get_a0, doc ="old lattice parameter a0")

    def _get_a(self):
        return self._a

    def _set_a(self, value):
        #self._a = value
        self.setLattice(ibrav = self._ibrav, a = float(value), base = self._base/self._a*value)

    a = property(_get_a, _set_a, doc ="lattice parameter a")


    def _get_b(self):
        return self._b

    def _set_b(self, value):
        #self._b = value
        self.setLattice(ibrav = self._ibrav, b = float(value))

    b = property(_get_b, _set_b, doc ="""lattice parameter b""")


    def _get_c(self):
        return self._c

    def _set_c(self, value):
        #self._c = value
        self.setLattice(ibrav = self._ibrav, c = float(value))

    c = property(_get_c, _set_c, doc ="""lattice parameter c""")


    def _get_cBC(self):
        return self._cBC

    def _set_cBC(self, value):
        #self._cBC = value
        self.setLattice(ibrav = self._ibrav, cBC = float(value))

    cBC = property(_get_cBC, _set_cBC, doc ="""lattice parameter cBC""")


    def _get_cAC(self):
        return self._cAC

    def _set_cAC(self, value):
        #self._cAC = value
        self.setLattice(ibrav = self._ibrav, cAC = float(value))

    cAC = property(_get_cAC, _set_cAC, doc ="""lattice parameter cAC""")


    def _get_cAB(self):
        return self._cAB

    def _set_cAB(self, value):
        #self._cAB = value
        self.setLattice(ibrav = self._ibrav, cAB = float(value))

    cAB = property(_get_cAB, _set_cAB, doc ="""lattice parameter cAB""")


    def _get_ibrav(self):
        return self._ibrav

    def _set_ibrav(self, value):
        if value < 0: value = 0
        ibravOld = self._ibrav
        self._ibrav = value
        if value == 0:
            base = self._base#/self._a
            if ibravOld != 4:
                self._type = 'generic cubic'
            else:
                self._type = 'generic hexagonal'
            self.setLattice(ibrav = self._ibrav, a = self._a, base = base)
        else:
            if 'generic' in self._type:
                self._type = 'celldm'
            base = self._getQEBaseFromParCos( ibrav = self._ibrav, \
                   a = self._a, b = self._b, c = self._c, cBC = self._cBC, \
                   cAC = self._cAC, cAB = self._cAB )
            qebase = base[0]*numpy.array(base[1])
            self.setLatticeFromQEVectors(self._ibrav, qebase)
#            self.setLattice(self._ibrav)

    ibrav = property(_get_ibrav, _set_ibrav, doc ="""Lattice symmetry parameter
ibrav. Changing it will update the lattice, but leave atomic fractional 
coordinates unchanged""")


    def _get_type(self):
        return self._type

    def _set_type(self, value):
        if 'generic' in value:
            self._type = value
            self._ibrav = 0
            base = self._base #/self._a
            self.setLattice(ibrav = self._ibrav, a = self._a, base = base)
        else:
            if self._ibrav == 0:
                pass
            else:
                self._type = value

    type = property(_get_type, _set_type, doc ="""QE lattice type: 'celldm',
    'traditional' or 'generic cubic', 'generic hexagonal'(implies ibrav = 0""")

    def _getQEBaseFromParCos( self, ibrav = 1, a = 1., b = 1., c = 1.,
                                    cBC = 0.,cAC = 0. ,cAB = 0.):
        c_a = float(c)/a
        # description dictionary of QE base vectors:
        QEBase = {
        # sc simple cubic:
        1 : (a, [[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]],
                 'Simple Cubic'),
        # fcc face centered cubic:
        2 : (a/2., [[-1, 0, 1],
                    [0, 1, 1],
                    [-1, 1, 0]],
                    'Face Centered Cubic'),
        # bcc body entered cubic:
        3 : (a/2., [[1, 1, 1],
                    [-1, 1, 1],
                    [-1, -1, 1]],
                    'Body Centered Cubic'), 
        # simple hexagonal and trigonal(p):
        4 : (a, [[1, 0, 0],
                 [-0.5, sqrt(3.0)/2.0, 0.],
                 [0,    0,          c_a]],
                 'Simple Hexagonal or Trigonal(P)'),
        # trigonal(r):
        5 : (a, [[sqrt((1.-cAB)/2.),-sqrt((1.-cAB)/6.), sqrt((1.+2.*cAB)/3.)],
                 [0, 2.*sqrt((1.-cAB)/6.),  sqrt((1.+2.*cAB)/3.)],
                 [-sqrt((1.-cAB)/2.), -sqrt((1.-cAB)/6.), sqrt((1.+2.*cAB)/3.)]],
                 'Trigonal(R)'),
        # simple tetragonal (p):
        6 : (a, [[1, 0, 0],
                 [0, 1, 0.],
                 [0, 0, c_a]],
                 'Simple Tetragonal(P)'),
        # body centered tetragonal (i):
        7 : (a/2., [[1, -1, c_a],
                    [1,  1, c_a],
                    [-1, -1, c_a]],
                    'Body Centered Tetragonal (I)'),
        # simple orthorhombic (p):
        8 : (1.0, [[a, 0., 0.],
                    [0., b, 0.],
                    [0., 0., c]],
                    'Simple Orthorhombic (P)'),  
        # bco base centered orthorhombic:
        9:  (1.0,  [[a/2., b/2., 0.],
                    [-a/2., b/2., 0.],
                    [0., 0., c]],
                    'Base Centered Orthorhombic'),
        # face centered orthorhombic:
        10: (1.0,  [[a/2., 0., c/2.],
                    [a/2., b/2., 0.],
                    [0., b/2., c/2.]],
                    'Face Centered Orthorhombic' ),
        # body centered orthorhombic:
        11: (1.0,  [[a/2., b/2., c/2.],
                    [-a/2., b/2., c/2.],
                    [-a/2., -b/2., c/2.]],
                    'Body Centered Orthorhombic'),
        # monoclinic (p):
        12: (1.0,  [[a, 0, 0],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [0, 0, c]],
                    'Monoclinic (P)'),
        # base centered monoclinic:
        13: (1.0,  [[a/2., 0, -c/2.],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [a/2., 0, c/2.]],
                    'Base Centered Monoclinic'),
        # triclinic:
        14: (1.0,  [[a, 0, 0],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [c*cAC, c*( cBC-cAC*cAB )/sqrt(1.-cAB**2), c*sqrt( 1. + 
                    2.*cBC*cAC*cAB - cBC**2 - cAC**2 - cAB**2)/sqrt(1.-cAB**2)]],
                    'Triclinic')
                    
        }
        return QEBase[ibrav]
    

if __name__ == '__main__':
    pass ;
