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


import unittest
import os
from qecalc.qetask.qeparser.pwinput import PWInput

# variables
tests_dir = os.path.dirname(os.path.abspath('testStructure.py'))
testdata_dir = os.path.join(tests_dir, 'data')

useStringConfig = True

stringConfig = """
&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    tstress = .true.,
    tprnfor = .true.,
    prefix = 'mgalb4',
    pseudo_dir = '/home/user/pslib',
    outdir = '/scratch/user',
/
&SYSTEM
    ibrav = 4,
    nbnd = 21,
    nspin = 1,
    occupations = 'smearing',
    degauss = 0.025,
    smearing = 'methfessel-paxton',
    ecutwfc = 64.0,
    ecutrho = 256.0,
    celldm(1) = 5.7889,
    celldm(3) = 2.22411,
    ntyp = 3,
    nat = 6,
/
&ELECTRONS
    conv_thr = 1.0d-10,
    mixing_beta = 0.4,
/
ATOMIC_SPECIES
 Al  26.9825 al_1.ncpp
 B   11.0000 B.pbe-n-van_ak.UPF
 Mg  24.3050 mg.ncpp
ATOMIC_POSITIONS (alat)
 Al      0.00000000  0.00000000  0.00000000
 B       0.50000000  0.28867510  0.52052140
 B      -0.00000000  0.57735030  0.52052140
 Mg      0.00000000  0.00000000  1.11205140
 B       0.50000000  0.28867510  1.70358160
 B      -0.00000000  0.57735030  1.70358160
K_POINTS (automatic)
 32 32 16 0 0 0
"""

class TestStructureMethods(unittest.TestCase):
    
    def setUp(self):
                
        
        if useStringConfig:
            self.input = PWInput( config = stringConfig )
        else:
            self.input = PWInput( filename = 'pw.in' )
            
        #self.input.autoUpdate = False
            
        self.input.parse()        


    def test_read_pwoutput(self):
        
        filename = os.path.join(testdata_dir, 'mgalb4_pw.out')        
        self.input.structure.read(filename, 'pwoutput')
        
        #answer = """"Simple Hexagonal or Trigonal(P)" cell    :\n 5.70000000  0.00000000  0.00000000\n-2.85000000  4.93634480  0.00000000\n 0.00000000  0.00000000  12.54000000\n\nAtomic positions in units of lattice parametr "a":\nAl      0.00000000  0.00000000  0.00000000  \nB       0.50000000  0.28867510  0.52052140  \nB       0.00000000  0.57735030  0.52052140  \nMg      0.00000000  0.00000000  1.11205140  \nB       0.50000000  0.28867510  1.70358160  \nB       0.00000000  0.57735030  1.70358160  \n\nAl  26.9825 al.ncpp\nB   11.0000 b.ncpp\nMg  24.3050 mg.ncpp\n"""
        answer = """"Simple Hexagonal or Trigonal(P)" cell:
 5.70000000  0.00000000  0.00000000
-2.85000000  4.93634480  0.00000000
 0.00000000  0.00000000  12.54000000

Atomic positions in units of lattice parametr "a":
Al      0.00000000  0.00000000  0.00000000  
B       0.50000000  0.28867510  0.52052140  
B       0.00000000  0.57735030  0.52052140  
Mg      0.00000000  0.00000000  1.11205140  
B       0.50000000  0.28867510  1.70358160  
B       0.00000000  0.57735030  1.70358160  

Al  26.9825 al.ncpp
B   11.0000 b.ncpp
Mg  24.3050 mg.ncpp
"""
    
        #print '*', str(self.input.structure), '*'
        #print '*', str(self.input.structure), '*'
        self.assertEqual(str(self.input.structure), answer)
        filename = os.path.join(testdata_dir, 'fev3_pwgeom.out')  
        self.input.structure.read(filename, 'pwoutput')

        # after geometry optimization:
        answer = """"Face Centered Cubic" cell:
-5.50788176  0.00000000  5.50788176
 0.00000000  5.50788176  5.50788176
-5.50788176  5.50788176  0.00000000

Atomic positions in units of lattice parametr "a":
V       0.00000000  0.00000000  0.00000000  
V       0.50000000  0.00000000  0.00000000  
V       0.25000000  0.25000000  0.25000000  
Fe      0.75000000  0.25000000  0.25000000  

V   50.9415 V.pbe-sp-van.UPF
Fe  55.8470 Fe.pbe-sp-van_ak.UPF
"""

        answer_toString="""&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    tstress = .true.,
    tprnfor = .true.,
    prefix = 'mgalb4',
    pseudo_dir = '/home/user/pslib',
    outdir = '/scratch/user',
/
&SYSTEM
    ibrav = 2,
    nbnd = 21,
    nspin = 1,
    occupations = 'smearing',
    degauss = 0.025,
    smearing = 'methfessel-paxton',
    ecutwfc = 64.0,
    ecutrho = 256.0,
    celldm(1) = 11.0157635161,
    celldm(2) = 1.0,
    celldm(3) = 1.0,
    celldm(4) = 0.0,
    ntyp = 2,
    nat = 4,
/
&ELECTRONS
    conv_thr = 1.0d-10,
    mixing_beta = 0.4,
/
ATOMIC_SPECIES
 V   50.9415 V.pbe-sp-van.UPF
 Fe  55.8470 Fe.pbe-sp-van_ak.UPF
ATOMIC_POSITIONS (alat)
 V       0.00000000  0.00000000  0.00000000
 V       0.50000000  0.00000000  0.00000000
 V       0.25000000  0.25000000  0.25000000
 Fe      0.75000000  0.25000000  0.25000000
K_POINTS (automatic)
 32 32 16 0 0 0
"""
        self.assertEqual(str(self.input.structure), answer)
        self.assertEqual(self.input.toString(), answer_toString)


    def test_read_pwinput(self):
        filename = os.path.join(testdata_dir, 'al_pw.in')
        self.input.structure.read(filename, 'pwinput')
        
        answer = """"Face Centered Cubic" cell:
-3.85000000  0.00000000  3.85000000
 0.00000000  3.85000000  3.85000000
-3.85000000  3.85000000  0.00000000

Atomic positions in units of lattice parametr "a":
Al      0.00000000  0.00000000  0.00000000  

Al  26.9800 Al.pz-vbc.UPF
"""

        answer_toString = """&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    tstress = .true.,
    tprnfor = .true.,
    prefix = 'mgalb4',
    pseudo_dir = '/home/user/pslib',
    outdir = '/scratch/user',
/
&SYSTEM
    ibrav = 2,
    nbnd = 21,
    nspin = 1,
    occupations = 'smearing',
    degauss = 0.025,
    smearing = 'methfessel-paxton',
    ecutwfc = 64.0,
    ecutrho = 256.0,
    celldm(1) = 7.7,
    celldm(2) = 1.0,
    celldm(3) = 1.0,
    celldm(4) = 0.0,
    ntyp = 1,
    nat = 1,
/
&ELECTRONS
    conv_thr = 1.0d-10,
    mixing_beta = 0.4,
/
ATOMIC_SPECIES
 Al  26.9800 Al.pz-vbc.UPF
ATOMIC_POSITIONS (alat)
 Al      0.00000000  0.00000000  0.00000000
K_POINTS (automatic)
 32 32 16 0 0 0
"""
        
        self.assertEqual(str(self.input.structure), answer)
        #print self.input.structure.toString()
        self.assertEqual(self.input.structure.toString(), answer_toString)
        
        
    def test_read_cif(self):
        filename = os.path.join(testdata_dir, 'PbTe.cif')
        self.input.structure.read(filename, 'cif')
        
        answer = """"generic" cell:
 12.20951961  0.00000000  0.00000000
 0.00000000  12.20951961  0.00000000
 0.00000000  0.00000000  12.20951961

Atomic positions in crystal coordinates:
Pb2+     0.50000000  0.50000000  0.50000000  
Pb2+     0.50000000  0.00000000  0.00000000  
Pb2+     0.00000000  0.50000000  0.00000000  
Pb2+     0.00000000  0.00000000  0.50000000  
Te      0.00000000  0.00000000  0.00000000  
Te      0.00000000  0.50000000  0.50000000  
Te      0.50000000  0.00000000  0.50000000  
Te      0.50000000  0.50000000  0.00000000  

Pb2+ 0.0000 
Te  0.0000 
"""
        self.assertEqual(str(self.input.structure), answer)

        filename = os.path.join(testdata_dir, 'LiNbO3.cif')
        self.input.structure.read(filename, 'cif')
        #self.input.structure.reduce(ibrav = 3)
        print self.input.toString()                   
   
        
    def test_load_diffpy(self): 
        from diffpy.Structure.structure import Structure
        from diffpy.Structure.atom import Atom
        from diffpy.Structure.lattice import Lattice        
        
        at1 = Atom('V', [0., 0., 0.])
        at2 = Atom('V', [0.5, 0., 0.])
        at3 = Atom('V', [0., 0.5, 0.])
        at4 = Atom('V', [0., 0., 0.5])
        at5 = Atom('V', [0.5, 0.5, 0.])
        at6 = Atom('V', [0., 0.5, 0.5])
        at7 = Atom('V', [0.5, 0., 0.5])
        at8 = Atom('V', [0.5, 0.5, 0.5])
    
        at9 = Atom('V', [0.25, 0.25, 0.25])
        at10 = Atom('Fe', [0.75, 0.25, 0.25])
        at11 = Atom('V', [0.75, 0.75, 0.25])
        at12 = Atom('Fe', [0.25, 0.75, 0.25])
    
        at13 = Atom('Fe', [0.25, 0.25, 0.75])
        at14 = Atom('V', [0.75, 0.25, 0.75])
        at15 = Atom('Fe', [0.75, 0.75, 0.75])
        at16 = Atom('V', [0.25, 0.75, 0.75])          
# set a in angstrom
        a =  2.*5.663/1.889725989
        struct = Structure( [ at1, at2, at3, at4, at5, at6, at7, at8, at9, \
                             at10, at11, at12, at13, at14, at15, at16], \
                             lattice = Lattice(a, a, a, 90, 90, 90) )
        #print struct
        massList = [50.9415, 55.847]
        psList  = ['V.pbe-n-van.UPF', 'Fe.pbe-nd-rrkjus.UPF']

        self.input.structure.load(source = 'diffpy', structure = struct, \
                                ibrav = 2, massList = massList, psList = psList)                
        
        answer1 = """"Face Centered Cubic" cell:
-5.66300000  0.00000000  5.66300000
 0.00000000  5.66300000  5.66300000
-5.66300000  5.66300000  0.00000000

Atomic positions in units of lattice parametr "a":
V       0.00000000  0.00000000  0.00000000  
V       0.50000000  0.00000000  0.00000000  
V       0.25000000  0.25000000  0.25000000  
Fe      0.75000000  0.25000000  0.25000000  

V   50.9415 V.pbe-n-van.UPF
Fe  55.8470 Fe.pbe-nd-rrkjus.UPF
"""
        
        #print self.input.toString()
        self.assertEqual(str(self.input.structure), answer1)

        self.input.structure.load(source = 'diffpy', structure = struct, \
                                massList = massList, psList = psList)
        
        answer2 = """"generic" cell:
 11.32600000  0.00000000  0.00000000
 0.00000000  11.32600000  0.00000000
 0.00000000  0.00000000  11.32600000

Atomic positions in units of lattice parametr "a":
V       0.00000000  0.00000000  0.00000000  
V       2.99673076  0.00000000  0.00000000  
V       0.00000000  2.99673076  0.00000000  
V       0.00000000  0.00000000  2.99673076  
V       2.99673076  2.99673076  0.00000000  
V       0.00000000  2.99673076  2.99673076  
V       2.99673076  0.00000000  2.99673076  
V       2.99673076  2.99673076  2.99673076  
V       1.49836538  1.49836538  1.49836538  
Fe      4.49509614  1.49836538  1.49836538  
V       4.49509614  4.49509614  1.49836538  
Fe      1.49836538  4.49509614  1.49836538  
Fe      1.49836538  1.49836538  4.49509614  
V       4.49509614  1.49836538  4.49509614  
Fe      4.49509614  4.49509614  4.49509614  
V       1.49836538  4.49509614  4.49509614  

V   50.9415 V.pbe-n-van.UPF
Fe  55.8470 Fe.pbe-nd-rrkjus.UPF
"""
        self.assertEqual(str(self.input.structure), answer2)


    def test_load_matter(self): 
        try:
            from matter import Structure, Atom, Lattice
        except ImportError:
            return       
        
        at1 = Atom('V', [0., 0., 0.])
        at2 = Atom('V', [0.5, 0., 0.])
        at3 = Atom('V', [0., 0.5, 0.])
        at4 = Atom('V', [0., 0., 0.5])
        at5 = Atom('V', [0.5, 0.5, 0.])
        at6 = Atom('V', [0., 0.5, 0.5])
        at7 = Atom('V', [0.5, 0., 0.5])
        at8 = Atom('V', [0.5, 0.5, 0.5])
    
        at9 = Atom('V', [0.25, 0.25, 0.25])
        at10 = Atom('Fe', [0.75, 0.25, 0.25])
        at11 = Atom('V', [0.75, 0.75, 0.25])
        at12 = Atom('Fe', [0.25, 0.75, 0.25])
    
        at13 = Atom('Fe', [0.25, 0.25, 0.75])
        at14 = Atom('V', [0.75, 0.25, 0.75])
        at15 = Atom('Fe', [0.75, 0.75, 0.75])
        at16 = Atom('V', [0.25, 0.75, 0.75])          
# set a in angstrom
        a =  2.*5.663/1.889725989
        struct = Structure( [ at1, at2, at3, at4, at5, at6, at7, at8, at9, \
                             at10, at11, at12, at13, at14, at15, at16], \
                             lattice = Lattice(a, a, a, 90, 90, 90) )
        #print struct
        massList = [50.9415, 55.847]
        psList  = ['V.pbe-n-van.UPF', 'Fe.pbe-nd-rrkjus.UPF']

        self.input.structure.load(source = 'matter', structure = struct, \
                                ibrav = 2, massList = massList, psList = psList)                
        
        answer1 = """"Face Centered Cubic" cell:
-5.66300000  0.00000000  5.66300000
 0.00000000  5.66300000  5.66300000
-5.66300000  5.66300000  0.00000000

Atomic positions in units of lattice parametr "a":
V       0.00000000  0.00000000  0.00000000  
V       0.50000000  0.00000000  0.00000000  
V       0.25000000  0.25000000  0.25000000  
Fe      0.75000000  0.25000000  0.25000000  

V   50.9415 V.pbe-n-van.UPF
Fe  55.8470 Fe.pbe-nd-rrkjus.UPF
"""
        
        self.assertEqual(str(self.input.structure), answer1)

        self.input.structure.load(source = 'matter', structure = struct, \
                                massList = massList, psList = psList)
        
        answer2 = """"generic" cell:
 11.32600000  0.00000000  0.00000000
 0.00000000  11.32600000  0.00000000
 0.00000000  0.00000000  11.32600000

Atomic positions in units of lattice parametr "a":
V       0.00000000  0.00000000  0.00000000  
V       2.99673076  0.00000000  0.00000000  
V       0.00000000  2.99673076  0.00000000  
V       0.00000000  0.00000000  2.99673076  
V       2.99673076  2.99673076  0.00000000  
V       0.00000000  2.99673076  2.99673076  
V       2.99673076  0.00000000  2.99673076  
V       2.99673076  2.99673076  2.99673076  
V       1.49836538  1.49836538  1.49836538  
Fe      4.49509614  1.49836538  1.49836538  
V       4.49509614  4.49509614  1.49836538  
Fe      1.49836538  4.49509614  1.49836538  
Fe      1.49836538  1.49836538  4.49509614  
V       4.49509614  1.49836538  4.49509614  
Fe      4.49509614  4.49509614  4.49509614  
V       1.49836538  4.49509614  4.49509614  

V   50.9415 V.pbe-n-van.UPF
Fe  55.8470 Fe.pbe-nd-rrkjus.UPF
"""
        self.assertEqual(str(self.input.structure), answer2)


    def test_list(self):
        
        answer = """Al   0.000000 0.000000 0.000000 26.982500 al_1.ncpp
B    0.666667 0.333333 0.234036 11.000000 B.pbe-n-van_ak.UPF
B    0.333333 0.666667 0.234036 11.000000 B.pbe-n-van_ak.UPF
Mg   0.000000 0.000000 0.499998 24.305000 mg.ncpp
B    0.666667 0.333333 0.765961 11.000000 B.pbe-n-van_ak.UPF
B    0.333333 0.666667 0.765961 11.000000 B.pbe-n-van_ak.UPF
"""
        s = ''
        for a in self.input.structure:
            s = s + str(a) + '\n'
        self.assertEqual(s, answer)     


    def test_diffpy(self):
        answer = """lattice=Lattice(base=array([[  5.7889    ,   0.        ,   0.        ],
       [ -2.89445   ,   5.01333446,   0.        ],
       [  0.        ,   0.        ,  12.87515038]]))
Al   0.000000 0.000000 0.000000 1.0000
B    0.666667 0.333333 0.234036 1.0000
B    0.333333 0.666667 0.234036 1.0000
Mg   0.000000 0.000000 0.499998 1.0000
B    0.666667 0.333333 0.765961 1.0000
B    0.333333 0.666667 0.765961 1.0000"""
        self.assertEqual(str(self.input.structure.diffpy()), answer)
        
        
    def test_readStr(self):
        al_pw_str = """&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    pseudo_dir = '/home/user/pslib',
    outdir = 'temp/',
    prefix = 'al',
    tprnfor = .true.,
    tstress = .true.,
/
&SYSTEM
    ibrav = 2,
    ecutwfc = 17.5,
    occupations = 'smearing',
    smearing = 'marzari-vanderbilt',
    degauss = 0.05,
    celldm(1) = 7.7,
    celldm(2) = 1.0,
    celldm(3) = 1.0,
    celldm(4) = 0.0,
    ntyp = 1,
    nat = 1,
/
&ELECTRONS
    diagonalization = 'david',
    mixing_beta = 0.7,
/
ATOMIC_SPECIES
 Al  26.9800 Al.pz-vbc.UPF
ATOMIC_POSITIONS (alat)
 Al      0.00000000  0.00000000  0.00000000
K_POINTS
 3
 0.0625000  0.0625000  0.0625000   1.00
 0.0625000  0.0625000  0.1875000   3.00
 0.0625000  0.0625000  0.3125000   3.00"""
        answer = """"Face Centered Cubic" cell:
-3.85000000  0.00000000  3.85000000
 0.00000000  3.85000000  3.85000000
-3.85000000  3.85000000  0.00000000

Atomic positions in units of lattice parametr "a":
Al      0.00000000  0.00000000  0.00000000  

Al  26.9800 Al.pz-vbc.UPF
"""
        self.input.structure.readStr(al_pw_str)
    
        self.assertEqual(str(self.input.structure), answer) 


    def test_fileInit(self):
        from qecalc.qetask.qeparser.qestructure import QEStructure
        filename = os.path.join(testdata_dir, 'al_pw.in')
        stru = QEStructure(filename = filename)
        answer = """"Face Centered Cubic" cell:
-3.85000000  0.00000000  3.85000000
 0.00000000  3.85000000  3.85000000
-3.85000000  3.85000000  0.00000000

Atomic positions in units of lattice parametr "a":
Al      0.00000000  0.00000000  0.00000000  

Al  26.9800 Al.pz-vbc.UPF
"""
        
        self.assertEqual(str(stru), answer)
        
        stru = QEStructure()
        filename = os.path.join(testdata_dir, 'PbTe.cif')
        stru.read(filename = filename, format = 'cif')
        
        answer ="""&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    tstress = .true.,
    tprnfor = .true.,
    prefix = 'mgalb4',
    pseudo_dir = '/home/user/pslib',
    outdir = '/scratch/user',
/
&SYSTEM
    ibrav = 0,
    nbnd = 21,
    nspin = 1,
    occupations = 'smearing',
    degauss = 0.025,
    smearing = 'methfessel-paxton',
    ecutwfc = 64.0,
    ecutrho = 256.0,
    celldm(1) = 1.889725989,
    ntyp = 2,
    nat = 8,
/
&ELECTRONS
    conv_thr = 1.0d-10,
    mixing_beta = 0.4,
/
ATOMIC_SPECIES
 Pb2+ 0.0000
 Te  0.0000
ATOMIC_POSITIONS (crystal)
 Pb2+     0.50000000  0.50000000  0.50000000
 Pb2+     0.50000000  0.00000000  0.00000000
 Pb2+     0.00000000  0.50000000  0.00000000
 Pb2+     0.00000000  0.00000000  0.50000000
 Te      0.00000000  0.00000000  0.00000000
 Te      0.00000000  0.50000000  0.50000000
 Te      0.50000000  0.00000000  0.50000000
 Te      0.50000000  0.50000000  0.00000000
K_POINTS (automatic)
 32 32 16 0 0 0
CELL_PARAMETERS (cubic)
 6.46100000  0.00000000  0.00000000
 0.00000000  6.46100000  0.00000000
 0.00000000  0.00000000  6.46100000
"""
        #print stru.toString( stringConfig )
        self.assertEqual(stru.toString( stringConfig ), answer )

        
    def test_writeStr(self):
        answer = """CRYST1    5.789    5.789   12.875  90.00  90.00 120.00                          
ATOM      1 Al           1       0.000   0.000   0.000  1.00  0.00          Al  
ATOM      2 B            1       2.894   1.671   3.013  1.00  0.00           B  
ATOM      3 B            1      -0.000   3.342   3.013  1.00  0.00           B  
ATOM      4 Mg           1       0.000   0.000   6.438  1.00  0.00          Mg  
ATOM      5 B            1       2.894   1.671   9.862  1.00  0.00           B  
ATOM      6 B            1      -0.000   3.342   9.862  1.00  0.00           B  
TER       7              1                                                      
END                                                                             
"""
        self.assertEqual(self.input.structure.writeStr(format = 'pdb'), answer)
        
    def test_constructor(self):
        from qecalc.qetask.qeparser.qelattice import QELattice
        from qecalc.qetask.qeparser.qestructure import QEStructure
        from qecalc.qetask.qeparser.qeatom import QEAtom
        
        filename = os.path.join(testdata_dir, 'fev3_pwgeom.out')  
        self.input.structure.read(filename, 'pwoutput')        
        
        vmass = 50.94150
        vpot = 'V_potential'
        femass = 55.84700
        fepot = 'Fe_potential'
        at1 = QEAtom('V', [0., 0., 0.], vmass, vpot)
        at2 = QEAtom('V', [0.5, 0., 0.], vmass, vpot)
        at3 = QEAtom('V', [0., 0.5, 0.], vmass, vpot)
        at4 = QEAtom('V', [0., 0., 0.5], vmass, vpot)
        at5 = QEAtom('V', [0.5, 0.5, 0.], vmass, vpot)
        at6 = QEAtom('V', [0., 0.5, 0.5], vmass, vpot)
        at7 = QEAtom('V', [0.5, 0., 0.5], vmass, vpot)
        at8 = QEAtom('V', [0.5, 0.5, 0.5], vmass, vpot)
    
        at9 = QEAtom('V', [0.25, 0.25, 0.25], vmass, vpot)
        at10 = QEAtom('Fe', [0.75, 0.25, 0.25], femass, fepot)
        at11 = QEAtom('V', [0.75, 0.75, 0.25], vmass, vpot)
        at12 = QEAtom('Fe', [0.25, 0.75, 0.25], femass, fepot)
    
        at13 = QEAtom('Fe', [0.25, 0.25, 0.75], femass, fepot)
        at14 = QEAtom('V', [0.75, 0.25, 0.75], vmass, vpot)
        at15 = QEAtom('Fe', [0.75, 0.75, 0.75], femass, fepot)
        at16 = QEAtom('V', [0.25, 0.75, 0.75], vmass, vpot)
        
        lattice = QELattice(lattice = self.input.structure.lattice)
        #print lattice
        new_struct = QEStructure( [ at1, at2, at3, at4, at5, at6, at7, at8, at9, \
                             at10, at11, at12, at13, at14, at15, at16], \
                             lattice = lattice )
        answer = """"Face Centered Cubic" cell:
-5.50788176  0.00000000  5.50788176
 0.00000000  5.50788176  5.50788176
-5.50788176  5.50788176  0.00000000

Atomic positions in crystal coordinates:
V       0.00000000  0.00000000  0.00000000  
V       0.50000000  0.00000000  0.00000000  
V       0.00000000  0.50000000  0.00000000  
V       0.00000000  0.00000000  0.50000000  
V       0.50000000  0.50000000  0.00000000  
V       0.00000000  0.50000000  0.50000000  
V       0.50000000  0.00000000  0.50000000  
V       0.50000000  0.50000000  0.50000000  
V       0.25000000  0.25000000  0.25000000  
Fe      0.75000000  0.25000000  0.25000000  
V       0.75000000  0.75000000  0.25000000  
Fe      0.25000000  0.75000000  0.25000000  
Fe      0.25000000  0.25000000  0.75000000  
V       0.75000000  0.25000000  0.75000000  
Fe      0.75000000  0.75000000  0.75000000  
V       0.25000000  0.75000000  0.75000000  

V   50.9415 V_potential
Fe  55.8470 Fe_potential
"""
        self.assertEqual(str(new_struct), answer)
        
        new_struct = QEStructure(self.input.structure)
        
        # Check if the constructor messes _qeInput        
        self.input.structure.lattice.a = 10              
        new_struct.lattice.a = 12
        s1 = self.input.structure.toString()
        s2 =  new_struct.toString()
        self.assertNotEqual(s1, s2)


    def test_reduce(self):
        filename = os.path.join(testdata_dir, 'PbTe.cif')
        self.input.structure.read(filename, 'cif')
        
        self.input.structure.reduce(ibrav = 2)
        
        answer = """"Face Centered Cubic" cell:
-6.10475981  0.00000000  6.10475981
 0.00000000  6.10475981  6.10475981
-6.10475981  6.10475981  0.00000000

Atomic positions in crystal coordinates:
Pb2+    -0.50000000  1.50000000 -0.50000000  
Te      0.00000000  0.00000000  0.00000000  

Pb2+ 0.0000 
Te  0.0000 
"""
        self.assertEqual(str(self.input.structure), answer)        
        
if __name__ == '__main__':
    #unittest.main()
    
    pwin_suite = unittest.TestSuite()
    pwin_suite.addTest(TestStructureMethods( "test_fileInit" ))
    
    pwinput_suite = unittest.TestSuite()
    pwinput_suite.addTest(TestStructureMethods( "test_read_pwinput" ))
    
    pwoutput_suite = unittest.TestSuite()
    pwoutput_suite.addTest(TestStructureMethods( "test_read_pwoutput" ))    
        
    unittest.main()
        
