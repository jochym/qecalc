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
from diffpy.Structure.lattice import Lattice, cosd

#from matter import Lattice
#from matter.Lattice import cosd

from math import sqrt, degrees
import numpy
from qeinput import QEInput


class QELattice(object):
    """Class QELattice for working with crystal lattices in QE notation
       Uses diffpy.Lattice class for storage
       Following parameters are dynamically linked to other properties 
       (E.g. lattice vectors):
       ibrav - lattice type
                if ibrav = 0, only 'a' prameter is relevant
       a, b, c, cBC, cAC ,cAB - lattice parameters
       type - lattice type to save into PWSCF cfg file ('celldm');
              'traditional' - using a,b,c,cosAC, cosAB, cosBC;
              'generic cubic', 'generic hexagonal' - assume exiasting
              section 'CELL_PARAMETERS', 'generic' types also 
              assume/set ibrav = 0
        setLattice() will set everything at once
       """

    def __init__(self, ibrav = 1,a = 1. ,b = 1.,c = 1.,
                 cBC = 0.,cAC = 0. ,cAB = 0., qeConf = None, base = None ):
#        Lattice.__init__(self)
        self.formatString = '%# .8f %# .8f %# .8f'
        self.qeConf = qeConf
        self._type = 'celldm'
        self._primitiveLattice = Lattice()
        self._standardLattice = Lattice()
        self._base = None
        self._a0 = None
        if self.qeConf != None:
            self.setLatticeFromPWInput(self.qeConf)
        else:
            if ibrav > 0 and base != None:
                self.setLatticeFromQEVectors(ibrav, base)
            else:
                self.setLattice(ibrav ,a ,b , c, cBC ,cAC ,cAB, base)

    # def _setDefaults(self):
        # self.formatString = '%# .8f %# .8f %# .8f'
        # self._type = 'celldm'
        # self._primitiveLattice = Lattice()
        # self._standardLattice = Lattice()        
        # self._base = None
        # self._a0 = None
        # self.setLattice(1. ,a ,b , c, cBC ,cAC ,cAB, base)
                
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


    def setLatticeFromPWInput(self, qeConf = None):
            
        if qeConf == None:
            qeConf = self.qeConf
        else:
            self.qeConf = qeConf
        if qeConf == None:
            raise NotImplementedError("writeLatticeToPWSCF: qeConf was not properly initialized")
    
    
        #self.setLattice(ibrav = 1,a = 1,b = 1,c = 1, cBC = 0,cAC = 0 ,cAB = 0)
        #if 'system' not in self.qeConf.namelists:
        #    if 'cell_parameters' not in self.qeConf.cards:
        #        return
        #    else:
        #jr        self.ibrav = 0
    
        if 'ibrav' in self.qeConf.namelists['system'].params:
            ibrav  = int(self.qeConf.namelist('system').param('ibrav'))
            if ibrav >= 0:
                a, b, c, cBC, cAC, cAB, base = self.parsePWInput(ibrav, qeConf)
            else:
                raise NotImplementedError("ibrav should be integer >= 0")
        else:
            raise NotImplementedError("config file should have ibrav defined")
        self.setLattice(ibrav, a, b, c, cBC, cAC, cAB, base)


    def setLatticeFromPWSCF(self, fname):
        self.qeConf = QEInput(fname)
        self.qeConf.parse()
        setLatticeFromPWInput(self, qeConf)

        
    def setLattice(self, ibrav, a = None, b = None, c = None,
                   cBC = None, cAC = None, cAB = None, base = None):
        """ 'base', numpy array of lattice vectors, and 'a'  will only be used
            if ibrav == 0. Otherwise, ibrav + lattice parameters will be used"""
        from math import acos
#        if [ibrav, a,b,c,cBC,cAC,cAB, base] == 8*[None]:
#            return None
        if ibrav == None:
            raise NonImplementedError('ibrav should be specified')
        self._ibrav = ibrav
        self._a0 = a
        if self._ibrav == 0:
#            print 'Found "generic" cell:'
            if base == None:
                raise NonImplementedError('base must be specified')
            if a == None: a = 1.0
            qeBase = numpy.array(base, dtype = float)*a
#            print qeBase
            self._a = 1.0
            if 'generic' not in self._type:
                self._type = 'generic cubic'
            self._primitiveLattice.setLatBase(qeBase)
            self._standardLattice.setLatBase(qeBase)
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
            self._primitiveLattice.setLatBase(qeBase)
            alpha = degrees(acos(self._cBC))
            beta = degrees(acos(self._cAC))
            gamma = degrees(acos(self._cAB))
            self._standardLattice.setLatPar(self._a,self._b,self._c,alpha,beta,gamma)            
        self._base = qeBase

    # def toString(self):
        # st = ''
        # if self._ibrav == 0:
            # st = st + '"generic" cell:\n'
        # else:
            # qeBaseTuple = self._getQEBaseFromParCos(self._ibrav, self._a, self._b,
                                               # self._c, self._cBC, self._cAC, self._cAB)
            # qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]
            # st = st + '"' + qeBaseTuple[2] + '" cell:\n'
        
        # for i in range(3):
            # v = self._primitiveLattice.base[i,:]
            # st = st + self.formatString%(v[0], v[1], v[2])
            # st = st + '\n'
      
        # return st

    def toString(self, string = None):        
        if string != None:        
            qeConf = QEInput(config = string)
            qeConf.parse()                
        else:
            if self.qeConf != None:
                qeConf = self.qeConf
            else:
                qeConf = QEInput(config = '')
        self.updatePWInput(qeConf)    
        return qeConf.toString()

    def latticeParams(self):
        return [self._a, self._b,self._c, self._cBC, self._cAC, self._cAB]


    def diffpy(self):
        '''Returns diffpy.Lattice object. Do not use it for reading  QE
        (standard cell) lattice parameters. Use latticeParams, or a, b, c , ...
        instead'''
        return self._primitiveLattice

    def parsePWInput(self, ibrav, qeConf = None):
        if qeConf == None:
            qeConf = self.qeConf
        if qeConf == None:
            raise NotImplementedError("writeLatticeToPWSCF: qeConf was not properly initialized")
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        # reset Lattice:
        #self.setLattice(1. ,1 , 1, cBC ,cAC ,cAB)
        if 'celldm(1)' in qeConf.namelists['system'].params:
            self._type = 'celldm' # celldm(i), i=1,6
            a = float(qeConf.namelist('system').param('celldm(1)'))

            if ibrav == 0:
                # lattice is set in the units of celldm(1)
                # need to parse CELL_PARAMETERS
                #if 'cell_parameters' not in qeConf.cards:
                #    return  #qeConf.createCard('cell_parameters')
                cellParLines = qeConf.card('cell_parameters').lines()
                #print cellParLines
                cellParType = qeConf.card('cell_parameters').arg()
                if cellParType == 'cubic' or cellParType == None:
                    self._type = 'generic cubic'
                else:
                    if cellParType == 'hexagonal':
                        self._type = 'generic hexagonal'
                # convert card into list
                base = []
                for line in cellParLines:
                    if '!' not in line:
                        words = line.split()
                        base.append([float(w) for w in words])
                return a, None, None, None, None, None, numpy.array(base)
            if ibrav > 0 and ibrav < 4:
                return a, a, a, cBC, cAC, cAB, None
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                return a, a, c_a*a, cBC, cAC, cAB, None
            if ibrav == 5:
                cAB = float(qeConf.namelist('system').param('celldm(4)'))
                return a, a, a, cAB, cAB, cAB, None
            if ibrav > 7 and ibrav < 12:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                return a, b_a*a, c_a*a, cBC, cAC, cAB, None
            if ibrav == 12 or ibrav == 13:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                cAB = float(qeConf.namelist('system').param('celldm(4)'))
                return a, b_a*a, c_a*a, cBC, cAC, cAB, None
            if ibrav == 14:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                cBC = float(qeConf.namelist('system').param('celldm(4)'))
                cAC = float(qeConf.namelist('system').param('celldm(5)'))
                cAB = float(qeConf.namelist('system').param('celldm(6)'))
                return a, b_a*a, c_a*a, cBC, cAC, cAB, None
        else:
            if ibrav == 0:
                print "Should specify celldm(1) if use 'generic' lattice"
                raise NotImplementedError
            a = float(qeConf.namelist('system').param('A'))
            self._type = 'traditional'   # A, B, C, cosAB, cosAC, cosBC
            if ibrav > 0 and ibrav < 4:
                return a, a, a, cBC, cAC, cAB, None
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c = float(qeConf.namelist('system').param('C'))
                return a, a, c, cBC, cAC, cAB, None
            if ibrav == 5:
                cAB = float(qeConf.namelist('system').param('cosAB'))
                return a, a, a, cAB, cAB, cAB, None
            if ibrav > 7 and ibrav < 12:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                return a, b, c, cBC, cAC, cAB, None
            if ibrav == 12 or ibrav == 13:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                cAB = float(qeConf.namelist('system').param('cosAB'))
                return a, b, c, cBC, cAC, cAB, None
            if ibrav == 14:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                cBC = float(qeConf.namelist('system').param('cosBC'))
                cAC = float(qeConf.namelist('system').param('cosAC'))
                cAB = float(qeConf.namelist('system').param('cosAB'))
                return a, b, c, cBC, cAC, cAB, None



    def getLatticeParamsFromPWSCF(self, ibrav, fname):
        qeConf = QEInput(fname)
        qeConf.parse()
        return self.parsePWInput(ibrav, qeConf)


    def updatePWInput(self, qeConf = None):
        if qeConf == None:
            qeConf = self.qeConf
        if qeConf == None:
            raise NotImplementedError("writeLatticeToPWSCF: qeConf was not properly initialized")
        if 'system' not in qeConf.namelists:
            qeConf.createNamelist('system')
        # clear geometry from qeConf:
        qeConf.namelist('system').remove('a')
        qeConf.namelist('system').remove('b')
        qeConf.namelist('system').remove('c')
        qeConf.namelist('system').remove('cosab')
        qeConf.namelist('system').remove('cosbc')
        qeConf.namelist('system').remove('cosac')
        qeConf.namelist('system').remove('celldm(1)')
        qeConf.namelist('system').remove('celldm(2)')
        qeConf.namelist('system').remove('celldm(3)')
        qeConf.namelist('system').remove('celldm(4)')
        qeConf.namelist('system').remove('celldm(5)')
        qeConf.namelist('system').remove('celldm(6)')
        if 'cell_parameters' in qeConf.cards:
            qeConf.removeCard('cell_parameters')
        if self._type == 'celldm':
            qeConf.namelist('system').add('ibrav', self._ibrav)
            qeConf.namelist('system').add('celldm(1)', self._a)
            qeConf.namelist('system').add('celldm(2)', self._b/self._a)
            qeConf.namelist('system').add('celldm(3)', self._c/self._a)
            if self._ibrav < 14:
                qeConf.namelist('system').add('celldm(4)', self._cAB)
            else:
                qeConf.namelist('system').add('celldm(4)', self._cBC)
                qeConf.namelist('system').add('celldm(5)', self._cAC)
                qeConf.namelist('system').add('celldm(6)', self._cAB)
        else:
            if self._type == 'traditional':
                qeConf.namelist('system').add('ibrav', self._ibrav)
                qeConf.namelist('system').add('A', self._a)
                qeConf.namelist('system').add('B', self._b)
                qeConf.namelist('system').add('C', self._c)
                qeConf.namelist('system').add('cosAB', self._cAB)
                qeConf.namelist('system').add('cosAC', self._cAC)
                qeConf.namelist('system').add('cosBC', self._cBC)
            else:
                if 'generic' in self._type:
                    qeConf.namelist('system').add('celldm(1)', 1.0)
                    self._ibrav = 0
                    qeConf.namelist('system').add('ibrav', self._ibrav)
                    if self._type == 'generic hexagonal':
                        cardArg = 'hexagonal'
                    if self._type == 'generic cubic' or self._type == None:
                        cardArg = 'cubic'
                    qeConf.createCard('cell_parameters')
                    qeConf.card('cell_parameters').setArg(cardArg)
                    qeConf.card('cell_parameters').removeLines()
                    for i in range(3):
                        v = self._primitiveLattice.base[i,:]
                        qeConf.card('cell_parameters').addLine(\
                                       self.formatString%(v[0], v[1], v[2]))



    def save(self, fname = None):
        """Will save the lattice either into its own file or into supplied with fname.
           It will also create all relevant sections/cards"""
        from os.path import exists
        filename = fname
        if fname != None:
            if not exists(filename):
                f = open(filename, 'w')
            qeConf = QEInput(fname)
            qeConf.parse()               
        else:
            qeConf = self.qeConf
            filename = qeConf.filename
            
        self.updatePWInput(qeConf)
        
        qeConf.save(filename)
        
        
    def recipCartesian(self, kPoint):
        """Conversts vector on fractional coordinates in reciprocal space into
           a vector in cartesian coordinates"""
        recip_base = self.diffpy().reciprocal().base*self._a
        return numpy.dot( kPoint, recip_base)

    def reciprocalBase(self):
        """
        Get reciprocal lattice vectors in units of 2*pi/a
        """
        return self.diffpy().reciprocal().base*self._a


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


    ####################################################################
    # property handlers
    ####################################################################

    # lattice parameters


    def _get_a0(self):
        if self._a0 != None:
            return self._a0
        else:
            return self._a
    a0 = property(_get_a0, doc ="old lattice parameter a0")

    def _get_a(self):
        return self._a

    def _set_a(self, value):
        self._a = value
        self.setLattice(ibrav = self._ibrav, a = self._a)

    a = property(_get_a, _set_a, doc ="lattice parameter a")


    def _get_b(self):
        return self._b

    def _set_b(self, value):
        self._b = value
        self.setLattice(ibrav = self._ibrav, b = self._b)

    b = property(_get_b, _set_b, doc ="""lattice parameter b""")


    def _get_c(self):
        return self._c

    def _set_c(self, value):
        self._c = value
        self.setLattice(ibrav = self._ibrav, c = self._c)

    c = property(_get_c, _set_c, doc ="""lattice parameter c""")


    def _get_cBC(self):
        return self._cBC

    def _set_cBC(self, value):
        self._cBC = value
        self.setLattice(ibrav = self._ibrav, cBC = self._cBC)

    cBC = property(_get_cBC, _set_cBC, doc ="""lattice parameter cBC""")


    def _get_cAC(self):
        return self._cAC

    def _set_cAC(self, value):
        self._cAC = value
        self.setLattice(ibrav = self._ibrav, cAC = self._cAC)

    cAC = property(_get_cAC, _set_cAC, doc ="""lattice parameter cAC""")


    def _get_cAB(self):
        return self._cAB

    def _set_cAB(self, value):
        self._cAB = value
        self.setLattice(ibrav = self._ibrav, cAB = self._cAB)

    cAB = property(_get_cAB, _set_cAB, doc ="""lattice parameter cAB""")


    def _get_ibrav(self):
        return self._ibrav

    def _set_ibrav(self, value):
        if value < 0: value = 0
        ibravOld = self._ibrav
        self._ibrav = value
        if value == 0:
            base = self._base/self._a
            if ibravOld != 4:
                self._type = 'generic cubic'
            else:
                self._type = 'generic hexagonal'
            self.setLattice(ibrav = self._ibrav, a = self._a, base = base)
        else:
            if 'generic' in self._type:
                self._type = 'celldm'
            self.setLatticeFromQEVectors(self._ibrav, self.diffpy().base)
#            self.setLattice(self._ibrav)

    ibrav = property(_get_ibrav, _set_ibrav, doc ="""Lattice symmetry parameter
                    ibrav""")


    def _get_type(self):
        return self._type

    def _set_type(self, value):
        if 'generic' in value:
            self._type = value
            self._ibrav = 0
            base = self._base/self._a
            self.setLattice(ibrav = self._ibrav, a = self._a, base = base)
        else:
            if self._ibrav == 0:
                pass
            else:
                self._type = value

    type = property(_get_type, _set_type, doc ="""QE lattice type: 'celldm',
    'traditional' or 'generic cubic', 'generic hexagonal'(implies ibrav = 0""")

if __name__ == '__main__':
    print "Hello";
