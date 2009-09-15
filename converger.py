from qecalc import QECalc
import numpy

class Converger(QECalc):
# value to converge with respect to k-points or different parameters in 'system'
# namelist
# currently can be 'total energy', 'single phonon', or 'geometry':
    def __init__(self, fname=None, taskName = 'total energy', tolerance = 1, nMaxSteps = 10):
        """taskName - currently can be 'total energy', 'single phonon',
           or 'geometry'
           tolerance -  task convergence criteria in percents
           nMaxSteps =  maximum number of optimization steps for
           the optimization routines"""

        # Default values, see explanations below:
#        convergerDic = {
#        'taskName': 'total energy',
#        'tolerance': '1',
#        'nMaxSteps': '10'
#        }
        QECalc.__init__(self,fname)
        # value to converge with respect to k-points or energy cutoffs
        # currently can be 'total energy', 'single phonon', or 'geometry':
#        self.taskName = self.config.get('Converger', 'taskName')

        # convergence criteria in percents:
#        self.tolerance = self.config.getfloat('Converger','tolerance')

        # maximum number of optimization steps:
#        self.nMaxSteps = self.config.getint('Converger','nMaxSteps')

        self.taskName = taskName
        self.tolerance = tolerance
        self.nMaxSteps = nMaxSteps
        
        self.lookupTable = {
        'total energy' : (self.pwscfLauncher, self.getTotalEnergy),
        'single phonon': (self.singlePhononLauncher, self.getSinglePhonon),
        'geometry'     : (self.pwscfLauncher, self.getLatticeParameters)
#        'multiple phonon': (self.multiPhononLauncher, self.getMultiPhonon)
        }
        assert self.lookupTable.has_key(self.taskName), "Convergence \
        estimator's name is not known"


    def converge(self, what, startValue, step, *params):
        """what - variable name from pwscf input, in case of k-points,
           what = 'kpoints'
           params - assumes starting value, increment, extraparameters(if
           required for given property)"""
        whatPossible = {'nbnd'         : 'system',
                        'degauss'      : 'system',
                        'ecutwfc'      : 'system',
                        'ecutrho'      : 'system',
                        'conv_thr'     : 'electrons',
                        'etot_conv_thr': 'control',
                        'forc_conv_thr': 'control',
                        'path_thr'     : 'ions',
                        'kpoints'      : 'k_points'
                        }
        if what not in whatPossible:
            raise Exception('Do not know how to converge that value!')
        value = numpy.array(startValue)
        step = numpy.array(step)
        runHistory = []
        for iStep in range(self.nMaxSteps):
            if what == 'kpoints':
                self.setkPointsAutomatic(value)
            else:
                self.qeConfig.namelist(whatPossible[what]).addParam(what, value)
            self.qeConfig.save()
            self.getLauncher()
            runHistory.append( self.getEstimator() )
            if iStep >= 2:
                if self.isConverged(runHistory): break
            value = value + step

        print 'optimized ' + what + ' value : ', value, '\n'
        print "Printing run history:\n", runHistory, '\n'
        print "End of convergence test\n"
        self.convergedValue = value
        return value

        

    def isConverged(self,runHistory):
        import math
        """Check for convergence:  two last runs should be less than
           the tolerance value if there is a list of values, the code will
           choose one with maximum error"""
        tol1 = []
        tol2 = []
        valTol = 1e-7
        for i in range( len(runHistory[-1]) ):
            # check if the denominator is not zerro:
            if math.fabs( runHistory[-2][i] ) > valTol and \
               math.fabs( runHistory[-2][i] ) > valTol :
                tol1.append( math.fabs( runHistory[-1][i]/runHistory[-2][i] - 1.0) )
                tol2.append( math.fabs( runHistory[-2][i]/runHistory[-3][i] - 1.0) )
        if max(tol1) < self.tolerance/100. and max(tol2) < self.tolerance/100.:
            print "\nSuccess! ",self.taskName,\
            " estimator value in two consecutive runs \ndiffers less than ", \
            self.tolerance, ' percent: ', max(tol2)*100,'%, ', \
            max(tol1)*100, '%\n'
            return True
        else:
            print runHistory[-1]
            return False

    def getLauncher(self):
        return self.lookupTable[self.taskName][0]()

    def getEstimator(self):
        return self.lookupTable[self.taskName][1]()
    

class EcutConverger(Converger):
    def __init__(self, fname=None):
        # Default values, see explanations below:
        configDic = {
        'isNormConserving': 'True',
        'ecutInit': '32',
        'ecutStep': '4'
        }
        Converger.__init__(self, fname)

        self.isNormConserving = self.config.getboolean('EcutConverger','isNormConserving')
        self.ecutInit = self.config.getfloat('EcutConverger','ecutInit')
        self.ecutStep = self.config.getfloat('EcutConverger','ecutStep')
        self.ecutConverger()

    def ecutConverger(self):
        #from parser.writetopwscf import varnameValue

        if self.isNormConserving:
            ecutrhoMult = 4.
        else:
            ecutrhoMult = 8.
        ecutwfc = self.ecutInit
        runHistory = []
        for iStep in range(self.nMaxSteps):
            ecutrho = ecutrhoMult*ecutwfc
            self.qeConfig.namelist('system').addParam('ecutwfc', ecutwfc)
            self.qeConfig.namelist('system').addParam('ecutrho', ecutrho)
            self.qeConfig.save()
#            varnameValue(self.pwscfInput,"ecutwfc", ecutwfc)
#            varnameValue(self.pwscfInput,"ecutrho", ecutrho)
            self.getLauncher()
            runHistory.append( self.getEstimator() )
            if iStep >= 2:
                if self.isConverged(runHistory): break
            ecutwfc = ecutwfc + self.ecutStep

        print "optimized ecut value : ", ecutwfc
        print runHistory


class KConverger(Converger):
    def __init__(self, fname=None):
        from string import split
        Converger.__init__(self, fname)

        self.isMetallic =  self.config.getboolean('KConverger','isMetallic')
        self.kInit = [ int(k) for k in split(self.config.get('KConverger','kInit')) ]
        self.kStep = [ int(k) for k in split(self.config.get('KConverger','kStep')) ]

        self.kConverger()

    def kConverger(self):
#        import parser.writetopwscf

        k_points = self.kInit
        runHistory = []
        for iK in range(self.nMaxSteps):
            self.setkPointsAutomatic(k_points)
#            writetopwscf.k_points(self.pwscfInput,k_points)
            self.getLauncher()
            runHistory.append( self.getEstimator() )
            if iK >= 2:
                if self.isConverged(runHistory): break
            for i in range(len(self.kStep)):
                k_points[i] = k_points[i] + self.kStep[i]
        print "optimized kpoints : ", k_points
        print runHistory
        return