import unittest
from qecalc.qetask.qeparser.pwinput import PWInput

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
            
        self.input.parse()
    
    
    def test_load_pwoutput(self):
        
        self.input.structure.load(source = 'pwoutput', \
                                  pwoutput = 'data/mgalb4_pw.out')
        
        answer = """"Simple Hexagonal or Trigonal(P)" cell:\n 5.70000000  0.00000000  0.00000000\n-2.85000000  4.93634480  0.00000000\n 0.00000000  0.00000000  12.54000000\n\nAtomic positions in units of lattice parametr "a":\nAl      0.00000000  0.00000000  0.00000000  \nB       0.50000000  0.28867510  0.52052140  \nB       0.00000000  0.57735030  0.52052140  \nMg      0.00000000  0.00000000  1.11205140  \nB       0.50000000  0.28867510  1.70358160  \nB       0.00000000  0.57735030  1.70358160  \n\nAl  26.9825 al.ncpp\nB   11.0000 b.ncpp\nMg  24.3050 mg.ncpp\n"""
                        
        self.assertEqual(str(self.input.structure), answer)
        
        self.input.structure.load(source = 'pwoutput', \
                                  pwoutput = 'data/fev3_pwgeom.out')
        
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
        self.assertEqual(str(self.input.structure), answer)
        
        
    def test_load_diffpy(self): pass
        

 
if __name__ == '__main__':
    unittest.main()
        