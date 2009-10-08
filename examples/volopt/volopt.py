#from task.writetopwscf import varnameValue, atomic_positions
#from task.task import Task
from qecalc import QECalc
import numpy as np
from scipy.optimize import brent
import scipy

def getHexEnergy(c, *args ):
    """ Total energy launcher for scipy "brent" routine
        'args' is a tuple with volume and task name 
        Volume value should be directly related to lattice constants: 
        E.g.: Vhex = a^2*c ommiting all the constant factors """
    volume = args[0]
    qe = args[1]
    qe.structure.lattice.a = np.sqrt(volume/c)
    qe.structure.lattice.c = c
    qe.structure.saveStructureToPWSCF()
    qe.structure.lattice.printBase()
    qe.pwscfLauncher()
    return qe.getTotalEnergy()[0]

# will find optimal a and c of hexagonal lattice of fixed volume:
def hexVolOpt(a0, c0_a0, volumeExpansion):
    '''provide initial(equilibrium) lattice parameters a0, c0/a0, and \
    desired volume expansion in percents at which one should run \
    the optimization '''
    qe = QECalc('config.ini')
    if qe.structure.lattice.ibrav != 4:
        raise Exception("The lattice must be hexagonal")
#   Initial(equilibrium) volume:    
    c0 = a0*c0_a0    
    volume = a0*a0*c0
            
    volume = volume + volume*volumeExpansion/100.
#   initial assumption: all latice parameters expand equally in percents    
    cExpansion = (1.+volumeExpansion/100.)**(1./3.)
    c = c0*cExpansion
    
    prcntC = 0.2 # percent of c we want to bracket around it for minima search(does not have to warantee the minima is inside)s

    brentOut = brent(getHexEnergy, (volume, qe), (c-c*prcntC/100, c+c*prcntC/100), tol = 1.e-7, full_output = 1)
    print brentOut
    c = brentOut[0]
    energy = brentOut[1]
    a = np.sqrt(volume/c)
    return a, c/a, energy

if __name__ == '__main__':

#    volPercRange = scipy.linspace(0.1, 3.0, 29)
    volPercRange = scipy.linspace(0.2, 2.4 , 12)
    # !!!!!!  Make sure you have correct starting scf.in at equilibrium
    qe = QECalc('config.ini')
    a0 = qe.structure.lattice.a
    c0 = qe.structure.lattice.c
    apar = [a0]
    c_apar = [c0/a0]
    # obtain total energy at equilibrium:
    qe.pwscfLauncher()
    energy = qe.getTotalEnergy()
    print volPercRange
    for volPrcnt in volPercRange:
         print "Optimizing volume at " + str(volPrcnt) + "% expansion"
         a, c_a, e = hexVolOpt(a0, c0/a0, volPrcnt)
         print e, a, c_a
         apar.append(a)
         c_apar.append(c_a)
         energy.append(e)
    print "a:"
    print apar
    print "c/a:"
    print c_apar
    print "Energy:"
    print energy
