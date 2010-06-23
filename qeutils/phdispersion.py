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

from qedispersion import QEDispersion
import numpy

class PHDispersion(QEDispersion):
    def __init__(self, lattice , matdynTask):
        QEDispersion.__init__(self,lattice)
        self.matdynTask = matdynTask

    def launch(self, *pathNPoints):
        #self.path, self.axis =
        self.setPath(*pathNPoints)

        self.matdynTask.input.qpoints.set(self.path, self.axis)
        self.matdynTask.input.save()

        self.matdynTask.launch()
        self.pol, self.dispersion, qpoints =  self.matdynTask.output.property('multi phonon')
        #self.pol = pol.reshape(pol.shape[0], pol.shape[1], pol.shape[2]*pol.shape[3])
        #self.dispersion = disp

    def setPhononPath(self, *pathNPoints):
        """
        Obsolete
        """
        self.launch(*pathNPoints)

    def _launchIndexFromLine(self, qstart, qend, nPoints = 40, thrsh_scl = 0.8):
        """
        Draws strait line between qstart and qend, runs matdyn job and solves for crossings
        along the line
        """
        print qstart, qend, nPoints
        x = numpy.linspace( qstart[0], qend[0], nPoints)
        y = numpy.linspace( qstart[1], qend[1], nPoints)
        z = numpy.linspace( qstart[2], qend[2], nPoints)

        path = numpy.column_stack((x,y,z))

        self.matdynTask.input.parse()
        self.matdynTask.input.qpoints.set(path)  
        self.matdynTask.input.save()

        self.matdynTask.launch()
        pol, dispersion, qpoints =  self.matdynTask.output.property('multi phonon')
        #print qpoints
        return  self._solveIndexFromLine(dispersion, pol, thrsh_scl)

    def _solveIndexFromLine(self, dispersion, pol, thrsh_scl = 0.8):
        """
        Provided dispersions and polarization vectors along a strait line,
        solves for indeces
        Returns index array of how to reindex the dispersions in order to take
        into account the crossings
        """
        idxs = [range(dispersion.shape[1])]
        print idxs
        for k in range(1, dispersion.shape[0]):
            idx = []
            scl = []
            for i, o1 in zip( idxs[k-1], dispersion[k-1]):
                for j, o2 in enumerate(dispersion[k]):
                    if ( abs((pol[k,j,:].conj()*pol[k-1,i,:]).sum()) >  thrsh_scl):
                   # if ( ( abs((self.pol[k,j,:].conj()*self.pol[k-1,i,:]).sum()) >  thrsh_scl) or \
                   # ( abs((self.pol[k,j,:]).real*(self.pol[k-1,i,:]).real).sum()  >  thrsh_scl) ):
                        #omg.append(o2)
                        scl.append( (pol[k,j,:].conj()*pol[k-1,i,:]).sum() )
                        idx.append(j)
                        break
            #if len(idx) != 9:
            if len(idx) != dispersion.shape[1]:
                print "OMG!!! ", k, self.dispersion[k-1], self.dispersion[k], idx, '\n'
                print 'scl = ', scl
                raise
            #dispersion.append(omg)
            idxs.append(idx)
        return idxs


    def solveAllCrossings(self, thrsh_scl = 0.8):
        self.pol = self.pol.reshape(self.pol.shape[0], self.pol.shape[1], self.pol.shape[2]*self.pol.shape[3])
        idxs = [range(self.dispersion.shape[1])]
        dk0 = self.path[1] - self.path[0]
        eps = 1.0e-10;
        for k in range(1, self.dispersion.shape[0]):
            dk = self.path[k] - self.path[k-1]
            #if ( abs((dk - dk0).sum()) > eps ):
            if 1 == 2:
                print "dk0!", dk0, dk
                dk0 = dk
                iii =  self._launchIndexFromLine( qstart = self.path[0], \
                       qend = self.path[k], nPoints = 300,  \
                       thrsh_scl = thrsh_scl)[-1]
                print k, self.path[k]
                print iii
                idxs.append( iii )
                continue
            idx = []
            scl = []
            for i, o1 in zip( idxs[k-1], self.dispersion[k-1]):
                for j, o2 in enumerate(self.dispersion[k]):
                    scl.append( (self.pol[k,j,:].conj()*self.pol[k-1,i,:]).sum() )
                    if ( abs((self.pol[k,j,:].conj()*self.pol[k-1,i,:]).sum()) >  thrsh_scl):
                        #scl.append( (self.pol[k,j,:].conj()*self.pol[k-1,i,:]).sum() )
                        idx.append(j)
                        break
            if len(idx) != 9:
                print "solveAllCrossings OMG!!! ", k, self.dispersion[k-1], \
                                                   self.dispersion[k], idx, '\n'
                print self.pol[k,idx,:]
                print self.pol[k-1,idxs[k-1],:]
                print 'scl = ', scl
                raise
            #dispersion.append(omg)
            idxs.append(idx)
        dispersion = []
        for k, idx in enumerate(idxs):
            dispersion.append( self.dispersion[k, idx])
    
        self.dispersion = numpy.array(dispersion)

    def solveAllCrossings2(self, thrsh_scl = 0.8):
        from numpy import real
        idxs = [range(self.dispersion.shape[1])]
        dk0 = self.path[1] - self.path[0]
        eps = 1.0e-10;
        for k in range(1, self.dispersion.shape[0]):
            dk = self.path[k] - self.path[k-1]
            #if ( abs((dk - dk0).sum()) > eps ):
            if 1 == 2:
                print "dk0!", dk0, dk
                dk0 = dk
                iii =  self._launchIndexFromLine( qstart = self.path[0], qend = self.path[k], nPoints = 300,  thrsh_scl = thrsh_scl)[-1]
                print k, self.path[k]
                print iii
                idxs.append( iii )
                continue
            idx = []
            #for ii in idxs[k-1]:
            #    print  (self.pol[k-1,4,:].conj()*self.pol[k-1,ii,:]).sum()
            for i, o1 in zip( idxs[k-1], self.dispersion[k-1]):
                mscl = 0.0
                mwdiff = 1.0
               # mscl = numpy.zeros(self.pol.shape[2])
                jj = -1
                for j, o2 in enumerate(self.dispersion[k]):
                    scl = []                    
                    weightk = []
                    weightkm1 = []
                    sc = 0.0
                    for iAtom in range(0, self.pol.shape[2]):
                        weightk.append( real(self.pol[k,j,iAtom,:].conj()*self.pol[k,j,iAtom,:]).sum() )
                        weightkm1.append( real(self.pol[k-1,i,iAtom,:].conj()*self.pol[k-1,i,iAtom,:]).sum() )
                        sc = sc + abs((self.pol[k,j,iAtom,:].conj()*self.pol[k-1,i,iAtom,:]).sum())

                    #if ( sc >  thrsh_scl ) :
                    #    scl.append( sc )
                    #else:
                        #scl = []
                    #    break


                    #wdiff = (abs( (numpy.array(weightk) - numpy.array(weightkm1))/numpy.array(weightkm1) )).sum()
                    wdiff = (abs( (numpy.array(weightk) - numpy.array(weightkm1)) )).sum()
                    #print wdiff
                    
                    #if wdiff < 0.1 and wdiff < mwdiff:
                    #    print 'weights  k miuns k-1:', wdiff
                    #    mwdiff = wdiff
                    #    jj = j
                    if  sc >  thrsh_scl:
                        jj = j
                        break
                    #if len(scl) == self.dispersion.shape[1]:
                    #    print len(scl), scl
                    #    if ((numpy.array(scl) - mscl).sum() > 0):
                    #    if ( sc > mscl ):
                    #        mscl = sc #numpy.array(scl)
                            #print mscl
                    #        jj = j

                    #sc = (((self.pol[k,j,:]*self.pol[k-1,i,:].conj()).sum()))
                    # implement weights criteria!
                    #if ( sc >  thrsh_scl ) : #and abs( (self.pol[k,j,:] - self.pol[k-1,i,:]).sum)):
                        #scl.append( (self.pol[k,j,:].conj()*self.pol[k-1,i,:]).sum() )
                        #if msc <= sc:
                        #    msc = sc
                        #    jj = j
                if jj >=0:
                    idx.append(jj)
                else:
                    break            
            if jj < 0:
                print "solveAllCrossings OMG!!! ", k, self.dispersion[k-1], self.dispersion[k], idx, '\n'
                #print self.pol[k,idx,:]
                #print self.pol[k-1,idxs[k-1],:]
                print 'mscl = ', mscl
                raise
            #dispersion.append(omg)
            print idx
            idxs.append(idx)
        dispersion = []
        for k, idx in enumerate(idxs):
            dispersion.append( self.dispersion[k, idx])
    
        self.dispersion = numpy.array(dispersion)


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 1:18:20 PM$"
