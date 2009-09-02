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
#    a = np.sqrt(volume/c)
#    geometry = ['Al       0.000000000   0.0000000000000000   0.000000000',
#                'B        0.500000000   0.2886751345948129   '+str(c/a/2.),
#                'B        0.000000000   0.5773502691896257   '+str(c/a/2.)]
#    varnameValue(task.pwscfInput, 'celldm(1)', a)
#    varnameValue(task.pwscfInput, 'celldm(2)', a)
#    varnameValue(task.pwscfInput, 'celldm(3)', c/a)
#    atomic_positions(task.pwscfInput, geometry)
    
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

    brentOut = brent(getHexEnergy, (volume, qe), (c-c*prcntC/100, c+c*prcntC/100), tol = 1.e-5, full_output = 1)
    print brentOut
    c = brentOut[0]
    energy = brentOut[1]
    a = np.sqrt(volume/c)
    return a, c/a, energy

#    stepPrcntC = 0.1 # percent
#    cRange = np.r_[c-c*prcntC/100:c+c*prcntC/100:c*stepPrcntC/100]
#    aRange = np.zeros(len(cRange))
#    aRange = aRange + volume 
#    aRange = np.sqrt(aRange/cRange)
#    print cRange
#    print aRange
#    energies = []    
#    for a, c in zip(aRange, cRange):
#        varnameValue(task.pwscfInput, 'celldm(1)', a)
#        varnameValue(task.pwscfInput, 'celldm(2)', a)
#        varnameValue(task.pwscfInput, 'celldm(3)', c/a)
#        task.getLauncher()
#        energies.append(task.getEstimator())
#    file = open('energies'+strvolumeExpansion)+'.txt','w')
#    for a,c, e in zip(aRange, cRange, energies):
#        print a, c, e
        
if __name__ == '__main__':

#    volPercRange = scipy.linspace(0.1, 3.0, 29)
    volPercRange = scipy.linspace(0.2, 1.2 , 6)

    qe = QECalc('config.ini')
    a0 = qe.structure.lattice.a
    c0 = qe.structure.lattice.c
    apar = []
    c_apar = []
    print volPercRange
    for volPrcnt in volPercRange:
         print "Optimizing volume at " + str(volPrcnt) + "% expansion"
         a, c_a, vol = hexVolOpt(a0, c0/a0, volPrcnt)
         print vol, a, c_a
         apar.append(a)
         c_apar.append(c_a)
    print "a:"
    print apar
    print "c/a:"
    print c_apar
