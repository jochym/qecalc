from diffpy.Structure.lattice import Lattice, cosd
from math import sqrt, degrees
import numpy
from parser.configParser import *


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

    def __init__(self, ibrav = 1,a = 1 ,b = 1,c = 1,
                 cBC = 0.,cAC = 0. ,cAB = 0., fname = None, base = None ):
#        Lattice.__init__(self)
        self.filename = fname
        self.__type = 'celldm'
        self.qeConf = None   # should be none if nothing to parse
        self.__primitiveLattice = Lattice()
        self.__standardLattice = Lattice()
        self.__base = None
        if self.filename != None:
            self.setLatticeFromPWSCF(self.filename)
        else:
            if ibrav > 0 and base != None:
                self.setLatticeFromQEVectors(ibrav, base)
            else:
                self.setLattice(ibrav ,a ,b , c, cBC ,cAC ,cAB, base)


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
            raise NotImplementedError
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


    def setLatticeFromPWSCF(self, fname = None):
        if fname != None:
            self.filename = fname
            self.qeConf = QEConfig(fname)
            self.qeConf.parse()
        if 'ibrav' in self.qeConf.namelists['system'].params:
            ibrav  = int(self.qeConf.namelist('system').param('ibrav'))
            if ibrav >= 0:
                a, b, c, cBC, cAC, cAB, base = self.getLatticeParamsFromPWSCF(ibrav, fname)
            else:
                raise NotImplementedError("ibrav should be integer >= 0")
        else:
            raise NotImplementedError("config file should have ibrav defined")
        self.setLattice(ibrav, a, b, c, cBC, cAC, cAB, base)


    def setLattice(self, ibrav, a = None, b = None, c = None,
                   cBC = None, cAC = None, cAB = None, base = None):
        """ 'base', numpy array of lattice vectors, and 'a'  will only be used
            if ibrav == 0. Otherwise, ibrav + lattice parameters will be used"""
        from math import acos
#        if [ibrav, a,b,c,cBC,cAC,cAB, base] == 8*[None]:
#            return None
        if ibrav == None:
            raise NonImplementedError('ibrav should be specified')
        self.__ibrav = ibrav
        if self.__ibrav == 0:
            print 'Found "generic" cell:'
            if base == None:
                raise NonImplementedError('base must be specified')
            if a == None: a = 1.0
            qeBase = numpy.array(base, dtype = float)*a
            print qeBase
            self.__a = 1.0
            self.__primitiveLattice.setLatBase(qeBase)
            self.__standardLattice.setLatBase(qeBase)
        else:
            if a is not None: self.__a = a
            if b is not None: self.__b = b
            if c is not None: self.__c = c
            if cBC is not None: self.__cBC = cBC
            if cAC is not None: self.__cAC = cAC
            if cAB is not None: self.__cAB = cAB
            qeBaseTuple = self.__getQEBaseFromParCos(self.__ibrav, self.__a, self.__b,
                                               self.__c, self.__cBC, self.__cAC, self.__cAB)
            qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]            
            print 'Found "' + qeBaseTuple[2] + '" cell'
            print 'Setting the base vectors according to QE conventions:'
            print qeBase
            self.__primitiveLattice.setLatBase(qeBase)
            alpha = degrees(acos(self.__cBC))
            beta = degrees(acos(self.__cAC))
            gamma = degrees(acos(self.__cAB))
            self.__standardLattice.setLatPar(self.__a,self.__b,self.__c,alpha,beta,gamma)
#            print "Standard Lattice:"
#            print self.__standardLattice.base
        self.__base = qeBase


    def latticeParams(self):
        return [self.__a, self.__b,self.__c, self.__cBC, self.__cAC, self.__cAB]


    def diffpy(self):
        '''Returns diffpy.Lattice object. Do not use it for reading  QE
        (standard cell) lattice parameters. Use latticeParams instead'''
        return self.__primitiveLattice


    def getLatticeParamsFromPWSCF(self, ibrav, fname):
        qeConf = QEConfig(fname)
        qeConf.parse()
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        if 'celldm(1)' in qeConf.namelists['system'].params:
            self.__type = 'celldm' # celldm(i), i=1,6
            a = float(qeConf.namelist('system').param('celldm(1)'))
            
            if ibrav == 0:
                # lattice is set in the units of celldm(1)
                # need to parse CELL_PARAMETERS
                cellParLines = qeConf.card('cell_parameters').getLines()                
                cellParType = qeConf.card('cell_parameters').argument()
                if cellParType == 'cubic' or cellParType == None:
                    self.__type = 'generic cubic'
                else:
                    if cellParType == 'hexagonal':
                        self.__type = 'generic hexagonal'
                # convert card into list
                base = []
                for line in cellParLines:
                    if '!' not in line:
                        words = line.split()
                        base.append([float(w) for w in words])
                return 1.0, None, None, None, None, None, base
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
            self.__type = 'traditional'   # A, B, C, cosAB, cosAC, cosBC
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


    def saveLatticeToPWSCF(self, fname = None):
        """Will save the lattice either into its own file or into supplied with fname.
           It will also create all relevant sections/cards"""
        if fname != None:
            filename = fname
            qeConf = QEConfig(fname)
            qeConf.parse()
        else:
            if self.filename != None:
                filename = self.filename
                qeConf = self.qeConf
            else:
                raise  # saveLatticeToPWSCF, filename was not supplied
            
        if qeConf == None:
            raise NotImplementedError("writeLatticeToPWSCF: qeConf was not properly initialized")
        if 'system' not in qeConf.namelists:
            qeConf.createNamelist('system')
        # clear geometry from qeConf:
        qeConf.namelist('system').removeParam('a')
        qeConf.namelist('system').removeParam('b')
        qeConf.namelist('system').removeParam('c')
        qeConf.namelist('system').removeParam('cosab')
        qeConf.namelist('system').removeParam('cosbc')
        qeConf.namelist('system').removeParam('cosac')        
        qeConf.namelist('system').removeParam('celldm(1)')
        qeConf.namelist('system').removeParam('celldm(2)')
        qeConf.namelist('system').removeParam('celldm(3)')
        qeConf.namelist('system').removeParam('celldm(4)')
        qeConf.namelist('system').removeParam('celldm(5)')
        qeConf.namelist('system').removeParam('celldm(6)')        
        if 'cell_parameters' in qeConf.cards:
            qeConf.removeCard('cell_parameters')
        if self.__type == 'celldm':
            qeConf.namelist('system').addParam('ibrav', self.__ibrav)
            qeConf.namelist('system').addParam('celldm(1)', self.__a)
            qeConf.namelist('system').addParam('celldm(2)', self.__b/self.__a)
            qeConf.namelist('system').addParam('celldm(3)', self.__c/self.__a)
            if self.__ibrav < 14:
                qeConf.namelist('system').addParam('celldm(4)', self.__cAB)
            else:
                qeConf.namelist('system').addParam('celldm(4)', self.__cBC)
                qeConf.namelist('system').addParam('celldm(5)', self.__cAC)
                qeConf.namelist('system').addParam('celldm(6)', self.__cAB)
        else:
            if self.__type == 'traditional':
                qeConf.namelist('system').addParam('ibrav', self.__ibrav)
                qeConf.namelist('system').addParam('A', self.__a)
                qeConf.namelist('system').addParam('B', self.__a)
                qeConf.namelist('system').addParam('C', self.__a)
                qeConf.namelist('system').addParam('cosAB', self.__cAB)
                qeConf.namelist('system').addParam('cosAC', self.__cAC)
                qeConf.namelist('system').addParam('cosBC', self.__cBC)
            else:
                if 'generic' in self.__type:
                    qeConf.namelist('system').addParam('celldm(1)', 1.0)
                    self.__ibrav = 0
                    qeConf.namelist('system').addParam('ibrav', self.__ibrav)
                    if self.__type == 'generic hexagonal':
                        cardArg = 'hexagonal'
                    if self.__type == 'generic cubic' or self.__type == None:
                        cardArg = 'cubic'
                    qeConf.createCard('cell_parameters')
                    qeConf.card('cell_parameters').setArgument(cardArg)
                    qeConf.card('cell_parameters').removeLines()
                    for i in range(3):
                        v = self.__primitiveLattice.base[i,:]
                        qeConf.card('cell_parameters').addLine(str(v)[1:-1])                                
        qeConf.save(filename)


    def __getQEBaseFromParCos( self, ibrav = 1, a = 1, b = 1, c = 1,
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

    def _get_a(self):
        return self.__a

    def _set_a(self, value):
        self.__a = value
        self.setLattice(ibrav = self.__ibrav, a = self.__a)

    a = property(_get_a, _set_a, doc ="lattice parameter a")


    def _get_b(self):
        return self.__b

    def _set_b(self, value):
        self.__b = value
        self.setLattice(ibrav = self.__ibrav, b = self.__b)

    b = property(_get_b, _set_b, doc ="""lattice parameter b""")


    def _get_c(self):
        return self.__c

    def _set_c(self, value):
        self.__c = value
        self.setLattice(ibrav = self.__ibrav, c = self.__c)

    c = property(_get_c, _set_c, doc ="""lattice parameter c""")


    def _get_cBC(self):
        return self.__cBC

    def _set_cBC(self, value):
        self.__cBC = value
        self.setLattice(ibrav = self.__ibrav, cBC = self.__cBC)

    cBC = property(_get_cBC, _set_cBC, doc ="""lattice parameter cBC""")


    def _get_cAC(self):
        return self.__cAC

    def _set_cAC(self, value):
        self.__cAC = value
        self.setLattice(ibrav = self.__ibrav, cAC = self.__cAC)

    cAC = property(_get_cAC, _set_cAC, doc ="""lattice parameter cAC""")


    def _get_cAB(self):
        return self.__cAB

    def _set_cAB(self, value):
        self.__cAB = value
        self.setLattice(ibrav = self.__ibrav, cAB = self.__cAB)

    cAB = property(_get_cAB, _set_cAB, doc ="""lattice parameter cAB""")


    def _get_ibrav(self):
        return self.__ibrav

    def _set_ibrav(self, value):
        if value < 0: value = 0
        self.__ibrav = value
        if value == 0:
            base = self.__base/self.__a
            self.__type = 'generic cubic'
            self.setLattice(ibrav = self.__ibrav, a = self.__a, base = base)
        else:
            if 'generic' in self.__type:
                self.__type = 'celldm'
            self.setLattice(ibrav = self.__ibrav)

    ibrav = property(_get_ibrav, _set_ibrav, doc ="""Lattice symmetry parameter
                    ibrav""")


    def _get_type(self):
        return self.__type

    def _set_type(self, value):
        if 'generic' in value:
            self.__type = value
            self.__ibrav = 0
            base = self.__base/self.__a
            self.setLattice(ibrav = self.__ibrav, a = self.__a, base = base)
        else:
            if self.__ibrav == 0:
                pass
            else:
                self.__type = value

    type = property(_get_type, _set_type, doc ="""QE lattice type: 'celldm',
    'traditional' or 'generic cubic', 'generic hexagonal'(implies ibrav = 0""")

if __name__ == '__main__':

    qeLattice = QELattice(fname = 'zro2.scf.in')
#    qeLattice = QELattice(fname = 'qwe.in')
 #   qeLattice.setLatticeFromPrimitiveVectors(12,qeLattice.lattice().base)
 #   print qeLattice.latticeParams()
#    qeLattice.latticeType = 'generic hexagonal'
#    qeLattice.bb = 12
    qeLattice.a = 13.0
    qeLattice.b = 24.0
    qeLattice.c = 3.
    qeLattice.ibrav = 4
    print qeLattice.b
    qeLattice.saveLatticeToPWSCF('./qwe.in')
#    qeLattice2 = QELattice()
#    qeLattice2.setLatticeFromPrimitiveVectors(qeLattice.ibrav, qeLattice.lattice().base )
    #print qeLattice.lattice().base
    #testLattice = Lattice(5,5,5,90,90,90)
    #print testLattice.base
