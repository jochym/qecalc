from qecalc.pwcalc import PWCalc
import numpy as np
from scipy.optimize import brent
import scipy
import os

def getHexEnergy(c, *args ):
    """ Total energy launcher for scipy "brent" routine
        'args' is a tuple with volume and task name
        Volume value should be directly related to lattice constants:
        E.g.: Vhex = a^2*c ommiting all the constant factors """
    volume = args[0]
    qe = args[1]
    qe.pw.input.structure.lattice.a = np.sqrt(volume/c)
    qe.pw.input.structure.lattice.c = c
    qe.pw.input.structure.save()
    qe.launch()
    qe.pw.input.structure.parseOutput(qe.pw.setting.pwscfOutput)
    qe.pw.input.structure.save()
    return qe.pw.output.property('total energy')[0]

# will find optimal a and c of hexagonal lattice of fixed volume:
def hexVolOpt(a0, c0_a0, volumeExpansion):
    '''provide initial(equilibrium) lattice parameters a0, c0/a0, and \
    desired volume expansion in percents at which one should run \
    the optimization '''
    qe = PWCalc('config.ini')
    qe.pw.input.parse()
    if qe.pw.input.structure.lattice.ibrav != 4:
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
    os.system('cp ' + qe.pw.setting.pwscfOutput + ' ' +  str(c) + qe.pw.setting.pwscfOutput)
    os.system('cp ' + qe.pw.setting.pwscfInput + ' ' +  str(c) + qe.pw.setting.pwscfInput)
    return a, c/a, energy

if __name__ == '__main__':

#    volPercRange = scipy.linspace(0.1, 3.0, 29)
#    volPercRange = scipy.linspace(-0.2, -1.0 , 5)
#    volPercRange = [2.6, 2.8, 3.0,3.2, 3.4, 3.6, 3.8, 4.0 ]
    volPercRange = [ 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0]
    # !!!!!!  Make sure you have correct starting scf.in at equilibrium
    pwcalc = PWCalc('config.ini')
    # name of file at eqilibrium is pwscfInput.eqv:
    eqvFileName = pwcalc.pw.setting.pwscfInput + '.eqv'
    if not os.path.exists(eqvFileName):
        raise Exception('Should provide PW input file at equilibrium')
    os.system('cp ' + eqvFileName + ' ' + pwcalc.pw.setting.pwscfInput )
    pwcalc.pw.input.parse()
    if pwcalc.pw.input.namelist('control').param('calculation') != "'relax'":
        print pwcalc.pw.input.namelist('control').param('calculation')
        raise Exception("""Should use "calculation = 'relax'" in "control" namelist""")
    a0 = pwcalc.pw.input.structure.lattice.a
    c0 = pwcalc.pw.input.structure.lattice.c
    apar = [a0]
    c_apar = [c0/a0]
    print apar
    print c_apar
    # obtain total energy at equilibrium:
    pwcalc.launch()
    energy = pwcalc.pw.output.property('total energy')
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
