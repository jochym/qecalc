from qecalc.qetask.pwtask import PWTask
from qecalc.qetask.qeparser.pwinput import PWInput
import numpy
from scipy.optimize import brent
import scipy
import os

class VolumeOptimizer:
    """
    Optimizes lattice parameters and performs atomic relaxation 
    at series of volume expansions
    
    pwInput.eqv (where pwInput is, e.g., scf.in) - template input file for 
    the task. Should have  "calculation = 'relax'" in the "control" namelistIt.
    This template is used at each combination of lattice parameters
    
    Optimized configurations, energies, and volume expansions are saved to a
    pickle file in a form of a dictionary:
        output = {
                  "volume": [],
                  "energy": [],
                  "config": []
                  }        
    """
    
    def __init__(self, pwTask):
        
        self.pickleName = 'volume_expansion.pkl'
        self.pw = pwTask       

    def launch(self, volPercRange):
        """
        Performs volume expansion
        volPercRange - list of percentage volume expansions relative to 
                       equilibrium volume. Can be negative.
        Equilibrium lattice parameters are provided through pwInput.eqv file \
        (e.g. If pwInput = scf.in, then scf.in.eqv)
        """
        # output container
        output = {
                  "volume": [],
                  "energy": [],
                  "config": []
                  }
        
        # name of file at eqilibrium is pwInput.eqv:
        eqvFileName = self.pw.setting.get('pwInput') + '.eqv'
        if not os.path.exists(eqvFileName):
            raise Exception('Should provide PW input file at equilibrium')
        os.system('cp ' + eqvFileName + ' ' + self.pw.setting.get('pwInput') )        
        self.pw.input.parse()
        if self.pw.input.namelist('control').param('calculation') != "'relax'":
            print self.pw.input.namelist('control').param('calculation')
            raise Exception("""Should use "calculation = 'relax'" in "control" \
                                                                    namelist""")
        # initialize equilibrium lattice
        pwInput0 = PWInput( filename = eqvFileName )
        pwInput0.parse()
        self.lattice0 = pwInput0.structure.lattice        
        
        output["config"].append( self.pw.input.toString() )      
        # obtain total energy at equilibrium:
        self.pw.launch()
        output["energy"] = self.pw.output.property('total energy')
        output["volume"] = [0.0]
        print volPercRange
        for volPrcnt in volPercRange:
             print "Optimizing volume at " + str(volPrcnt) + "% expansion"
             cfg, e = self._volOpt(volPrcnt)
             print cfg
             print e
             output["config"].append( cfg )
             output["energy"].append( e )
             output["volume"].append( volPrcnt )
        
        print output                
        import pickle       
        pickle.dump( output, open(self.pickleName, 'wb') )


    def _volOpt(self, volumeExpansion):
        '''
           volumeExpansion -  percent volume expansion relative to equilibrium
           at which the optimization is performed
        '''        
        ibrav = self.lattice0.ibrav
        if ibrav < 1 or ibrav > 4:
            raise Exception("This lattice type is not implemented")
        
        a0 = self.lattice0.a
        c0 = self.lattice0.c
        
        if ibrav == 4:
        # will find optimal a and c of hexagonal lattice of fixed volume:
            volume = a0*a0*c0
            volume = volume + volume*volumeExpansion/100.
            # initial assumption: all latice parameters expand equally in percents
            c = c0*(1.+volumeExpansion/100.)**(1./3.)
    
            # percent of c we want to bracket around it for minima 
            # search(does not have to guarantee the minima is inside
            prcntC = 0.2 
    
            brentOut = brent(self._getHexEnergy, (volume,), (c-c*prcntC/100, \
                                  c+c*prcntC/100), tol = 1.e-7, full_output = 1)
            print brentOut
            c = brentOut[0]
            energy = brentOut[1]            
            # relax structure at optimized parameters to get optimized atomic positions
            relax_energy = self._relax( a = numpy.sqrt(volume/c), c = c )
            
        if ibrav > 0 and ibrav < 4:
            aExpansion = (1.+volumeExpansion/100.)**(1./3.)
            a = a0 * (1.+volumeExpansion/100.)**(1./3.)
            volume = a0*a0*c0
            volume = volume + volume*volumeExpansion/100.
            prcntA = 0.2 
            brentOut = brent(self._getCubicEnergy, (volume,), (a-a*prcntA/100, \
                                  a+a*prcntA/100), tol = 1.e-7, full_output = 1)
            energy = brentOut[1]
            
        print "Double check: Brent energy = %f, Relax energy = %f"%(energy, relax_energy)       
        
        os.system('cp ' + self.pw.setting.get('pwOutput') + ' ' +  \
                             str(volumeExpansion) + self.pw.setting.get('pwOutput'))
        os.system('cp ' + self.pw.setting.get('pwInput') + ' ' +  \
                              str(volumeExpansion) + self.pw.setting.get('pwInput'))            
            
        return self.pw.input.toString(), energy


    def _getHexEnergy(self, c, *args ):
        """ Total energy launcher for scipy "brent" routine
            'args' is a tuple with volume and task name
            Volume value should be directly related to lattice constants:
            E.g.: Vhex = a^2*c omitting all the constant factors """
        volume = args[0]
        return self._relax(a = numpy.sqrt(volume/c), c = c)
    
    
    def _getCubicEnergy(self, a, *args):
        """ 
        Total energy launcher for scipy "brent" routine for cubic lattice
        """
        return self._relax(a = a)
    
        
    def _relax(self, ibrav = None, a = None, b = None, c = None, \
                                     cosBC = None ,cosAC = None, cosAB = None):
        """
        Runs geometry relaxation with given lattice parameters. Modifies structure
        and returns total energy
        """
        if ibrav == None:
            ibrav = self.pw.input.structure.lattice.ibrav
        self.pw.input.structure.lattice.setLattice(ibrav = ibrav, a = a, b = b,\
                                  c = c, cBC = cosBC, cAC = cosAC, cAB = cosAB)
        self.pw.input.structure.save()
        self.pw.launch()
        self.pw.input.structure.read(filename = self.pw.setting.get('pwOutput'), format='pwoutput')
        self.pw.input.save()
        return self.pw.output.property('total energy')[0] 


if __name__ == '__main__':

    volPercRange = [ 2.0, 2.4, 2.8, 3.2, 3.6, 4.0]   
    VolumeOptimizer(pwTask = PWTask('config.ini')).launch( volPercRange = volPercRange )
