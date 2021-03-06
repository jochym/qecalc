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
from pwcalc import PWCalc
from singlephononcalc import SinglePhononCalc
import numpy

class Converger():
# value to converge with respect to k-points or different parameters in 'system'
# namelist
# currently can be 'total energy', 'single phonon', or 'geometry':
    def __init__(self, fname=None, taskName = 'total energy', tolerance = 1, nMaxSteps = 10):
        """taskName - currently can be 'total energy', 'single phonon',
           or 'geometry'
           tolerance -  task convergence criteria in percents
           nMaxSteps =  maximum number of optimization steps for
           the optimization routines"""

        # define calcs:
        self.pwCalc = PWCalc(fname)
        self.singlePhononCalc = SinglePhononCalc(fname)

        self.taskName = taskName
        self.tolerance = tolerance
        self.nMaxSteps = nMaxSteps

        # defines lookupTable: specifies capable Calc and property of interest
        # which can be looked up in output parsers
        self.lookupTable = {
        'total energy' : (self.pwCalc, 'total energy'),
        'fermi energy' : (self.pwCalc, 'fermi energy'),
        'single phonon': (self.singlePhononCalc, 'single phonon'),
        'geometry'     : (self.pwCalc, 'lattice parameters')
        }
        assert self.lookupTable.has_key(self.taskName), "Convergence \
        estimator's name is not known"


    def converge(self, what, startValue, step = None, multiply = None):
        """what - variable name from pwscf input, in case of k-points,
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
        calc = self.lookupTable[self.taskName][0]
        # this implies all available calcs have pw task in them:
        calc.pw.input.parse()

        if what not in whatPossible:
            raise Exception('Do not know how to converge that value!')
        if step == None and multiply == None:
            raise Exception("Should set either 'step' or 'multiply'")

        value = numpy.array(startValue)
        step = numpy.array(step)
        
        runHistory = []
        for iStep in range(self.nMaxSteps):
            if what == 'kpoints':
                calc.pw.input.kpoints.setAutomatic(value)
            else:
                calc.pw.input.namelist(whatPossible[what]).addParam(what, value)
            calc.pw.input.save()
            calc.launch()
            propertyName = self.lookupTable[self.taskName][1]
            print '\nStep: ', iStep
            print what, ': ', value
            print propertyName + ': ', calc.lookupProperty(propertyName)
            runHistory.append( calc.lookupProperty(propertyName) )
            if iStep >= 2:
                if self.isConverged(runHistory): break
            if multiply !=None:
                value = value*numpy.array(multiply)
            else:
                value = value + step
            # if the run includes geometry optimization - import optimized
            # structure othervise it will reimport existing structure:
            calc.pw.input.structure.parseOutput(calc.pw.setting.pwscfOutput)
            calc.pw.input.structure.save()

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
            #print runHistory[-1]
            return False
    