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
from qecalc.pwcalc import PWCalc
from qecalc.singlephononcalc import SinglePhononCalc
import numpy
from qecalc.qetask.setting import Setting

class Converger(Setting):
# value to converge with respect to k-points or different parameters in 'system'
# namelist
# currently can be 'total energy', 'single phonon', or 'geometry':
    def __init__(self, filename = None, configString = None, sectionName = 'Converger',\
                      taskName = None, tolerance = 1, nMaxSteps = 10, cleanOutDir = True):
        """taskName - currently can be 'total energy', 'single phonon',
           tolerance -  task convergence criteria in percents
           nMaxSteps =  maximum number of optimization steps for
           the optimization routines"""

        Setting.__init__( self, filename, configString )
        self.history = []

        self._filename = filename
        self._configString = configString

        if taskName == None:
            import ConfigParser
            configDic = {
                            'taskName'  : 'total energy',
                            'tolerance ': '0.1',
                            'nMaxSteps' : '10',
                            'what'      : 'ecutwfc',
                            'startValue': None,
                            'step'      : None,
                            'multiply'  : None,
                        }
        
            self.section(sectionName, configDic)
            self.tolerance = float(self.tolerance)
            self.nMaxSteps = int(self.nMaxSteps)
            if self.what == 'kpoints':
                self.startValue = map(int, self.startValue.split())
                if self.step != None:
                    self.step = map(int, self.step.split())
            else:
                self.startValue = float(self.startValue)
                if self.step != None:
                    self.step = float(self.step)
            if self.multiply != None:
                self.multiply = float(self.multiply)

        else:
            self.taskName = taskName
            self.tolerance = tolerance
            self.nMaxSteps = nMaxSteps

        # defines lookupTable: specifies capable Calc and property of interest
        # which can be looked up in output parsers
        self.lookupTable = {
        'total energy' : (PWCalc, 'total energy'),
        'fermi energy' : (PWCalc, 'fermi energy'),
        'single phonon': (SinglePhononCalc, 'single phonon'),
        'geometry'     : (PWCalc, 'lattice parameters')
        }
        assert self.lookupTable.has_key(self.taskName), "Convergence \
        estimator's name is not known"


    def converge(self, what = None, startValue = None, step = None, multiply = None):
        """what - variable name from pw input, in case of k-points,
           what = 'kpoints'
           params -  extraparameters(if required for given property)"""
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
        self.runHistory = []
        if what == None:
            # look up in config file
            what = self.what
            startValue = self.startValue
            step = self.step
            multiply = self.multiply
        else:
            if startValue == None and step == None or\
               startValue == None and multiply == None:
                raise('Converger.converge: starting value and increment\
                                                         were not properly set')
                
        self._initCalc()
        calc = self.calc

        if what not in whatPossible:
            raise Exception('Do not know how to converge that value!')
        if step == None and multiply == None:
            raise Exception("Should set either 'step' or 'multiply'")

        value = numpy.array(startValue)
        step = numpy.array(step)

        for iStep in range(self.nMaxSteps):
            if what == 'kpoints':
                calc.pw.input.kpoints.setAutomatic(value)
            else:
                calc.pw.input.namelist(whatPossible[what]).set(what, value)
            print '\n\nStep: ', iStep
            print what, ': ', value
            calc.pw.input.save()
            calc.launch()
            propertyName = self.lookupTable[self.taskName][1]
            print '\n', propertyName + ': ', calc.lookupProperty(propertyName)
            self.runHistory.append( calc.lookupProperty(propertyName) )
            if iStep >= 2:
                if self.isConverged(self.runHistory): break
            if multiply !=None:
                value = value*numpy.array(multiply)
            else:
                value = value + step
            # if the run includes geometry optimization - import optimized
            # structure othervise it will reimport existing structure:
            calc.pw.input.structure.read(filename = \
                           calc.pw.setting.get('pwOutput'), format = 'pwoutput')
            calc.pw.input.structure.save()

        print 'optimized ' + what + ' value : ', value, '\n'
        print "Printing run history:\n", self.runHistory, '\n'
        print "End of convergence test\n"
        self.convergedValue = value
        self.history.append(self.runHistory)
        return value


    def _initCalc(self):
        
        self.calc = self.lookupTable[self.taskName][0]( \
                   filename = self._filename, configString = self._configString)        
        
        self.calc.pw.cleanOutDir = True
        
    

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
            #print runHistory[-1]
            return False
