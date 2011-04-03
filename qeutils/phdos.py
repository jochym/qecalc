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
import numpy
from qedos import QEDOS
from qeutils import kmesh

class PhononDOS(QEDOS):

    def __init__(self, matdynTask, structure = None):
        QEDOS.__init__(self)

        self.structure = structure
        self._freqs = None
        self._modes = None
        self._qpts = None
        self.axis = []
        self.dos = []
        self.partdos = {} # OrderedDict on init in PartDOS {label: pdos}
        self.partdosElem = {} # OrderedDict on init in PartDOS {symbol: pdos}

        self.matdynTask = matdynTask


    def launch(self, nqpoints, partialDOS = False):
        """
        launches matdyn task with a grid provided through list 'nqpoints'
        """
        if partialDOS == False:
            self.matdynTask.qpoints.setAutomatic(nqpoints)
            self.matdynTask.input.save()

            self.matdynTask.launch()
            self.loadPhonons()
            self.axis, self.dos = self.matdynTask.output.property('phonon dos')
        else:
            if self.structure == None:
                raise('PartialDOS: structure was not set')
            qpoints = kmesh.kMeshCart(nqpoints, \
                                         self.structure.lattice.reciprocalBase())

            #update qpoints and launch matdyn
            self.matdynTask.syncSetting()
            self.matdynTask.input.qpoints.set(qpoints)
            self.matdynTask.input.save()
            self.matdynTask.launch()
            self.loadPhonons()
            #self.axis, self.dos = self.DOS()
            #self.axis, self.partdos = self.partDOS()

    def loadPhonons(self):
        self.matdynTask.syncSetting()
        self.matdynTask.output.parse()
        self._modes, self._freqs, self._qpts =  \
                                self.matdynTask.output.property('multi phonon')

    def setPhonons(self, modes, freqs, qpts):
        self._modes, self._freqs, self._qpts =   modes, freqs, qpts


    def setRange(self, minOmega, maxOmega, deltaOmega):
        if minOmega == None:
            minOmega = numpy.min(self._freqs)
        if maxOmega == None:
            maxOmega = numpy.max(self._freqs)
        if minOmega > maxOmega: minOmega, maxOmega = maxOmega, minOmega
        if deltaOmega == None:
            deltaOmega = (maxOmega - minOmega)/200.0
        return minOmega, maxOmega, deltaOmega


    def DOS(self, minOmega = None, maxOmega = None, deltaOmega = None):
        minOmega, maxOmega, deltaOmega =  \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histOmega = numpy.zeros(nPoints)
        norm = 0.0
        for cell_freqs in self._freqs:
            for omega in cell_freqs:
                idx = int( (abs(omega) - minOmega)/deltaOmega )
                if idx < len(histOmega):
                    histOmega[idx] = histOmega[idx] + 1.0
                    norm = norm + 1.0

        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        self.axis, self.dos = axis, histOmega/norm
        return  self.axis, self.dos

    def partDOSType(self, atomSymbol, minOmega = None, maxOmega = None, \
                                                            deltaOmega = None):
        from numpy import real
        raise NonImplementedError
        minOmega, maxOmega, deltaOmega   =      \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histPartOmega = numpy.zeros(nPoints)
        norm = 0.0
        for iAtom, atom in enumerate(self.structure.matter()):
            if atomSymbol == atom.element:
                for cell_freqs, vectors in zip(self._freqs, self._modes):
                    for omega, vector in zip(cell_freqs, vectors[:,iAtom,:]):
                        idx = int( (abs(omega) - minOmega)/deltaOmega )
                        if idx < len(histPartOmega):
                            weight = (real(vector*vector.conjugate())).sum()
                            histPartOmega[idx] = histPartOmega[idx] + weight
                            norm = norm + weight
        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        self.axis, self.partdos = axis, histOmega/norm
        return axis, histPartOmega/norm

    def partDOSAtom(self, iAtom, minOmega = None, maxOmega = None, \
                                                            deltaOmega = None):
        from numpy import real
        from qecalc.qetask.qeparser.orderedDict import OrderedDict
        minOmega, maxOmega, deltaOmega   =      \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histPartOmega = numpy.zeros(nPoints)
        labels = self.structure.matter().getLabels()
        norm = 0.0

        if self.partdos == {}:
            self.partdos = OrderedDict()
            for i in range( len(self.structure.matter()) ):
                self.partdos[labels[i]] = []
        atom = self.structure.matter()[iAtom]
        for cell_freqs, vectors in zip(self._freqs, self._modes):
            for omega, vector in zip(cell_freqs, vectors[:,iAtom,:]):
                idx = int( (abs(omega) - minOmega)/deltaOmega )
                if idx < len(histPartOmega):
                    weight = (real(vector*vector.conjugate())).sum()
                    histPartOmega[idx] = histPartOmega[idx] + weight
                    norm = norm + weight
        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        self.axis, self.partdos[labels[iAtom]] = axis, histPartOmega/norm
        return self.axis, self.partdos[labels[iAtom]]

    def partDOS(self, iAtom = None, atomSymbol = None, minOmega = None, \
                                            maxOmega = None, deltaOmega = None):
        """
        calculates partial phonon DOS.
        iAtom - number of atom in structure (partDOS will calculate
                                             its contribution only)
        atomSymbol - atomic symbol (partDOS will calculate
                     contribution of all atoms with this symbol)
        iAtom and atomSymbol are not provided, ordered dictionary
        of histograms with unic labels is created: { 'label1': hist
                                                         ...
                                                    }
        """
        from qecalc.qetask.qeparser.orderedDict import OrderedDict
        minOmega, maxOmega, deltaOmega   =      \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histPartOmega = numpy.zeros(nPoints)
        if iAtom != None:
            return self.partDOSAtom(iAtom, minOmega, maxOmega, deltaOmega)

        if atomSymbol != None:
            for iAtom, atom in enumerate(self.structure.matter()):
                if atomSymbol == atom.element:
                    axis, hist = self.partDOSAtom(iAtom, minOmega, maxOmega, \
                                                                   deltaOmega)
                    histPartOmega = histPartOmega + hist
            return axis, histPartOmega
        # will return OrderedDict of histograms:
        self.dos = numpy.zeros(nPoints)
        self.partdosElem = OrderedDict()
        for atom in self.structure.matter():
            self.partdosElem[atom.element] = numpy.zeros(nPoints)
        for iAtom, atom in enumerate(self.structure.matter()):
            axis, hist = self.partDOSAtom(iAtom, minOmega, maxOmega, deltaOmega)
            self.dos = self.dos + hist
            self.partdosElem[atom.element] = self.partdosElem[atom.element] + \
                                                                            hist
        return axis, self.partdos

    def save(self):
        """
        saves all DOSes
        """
        if self.dos != []:
            numpy.savetxt('dos_total.dat', numpy.column_stack( (self.axis, \
                                                                    self.dos)) )
        if self.partdos != {}:
            for label in self.partdos:
                numpy.savetxt('dos_' + label + '.dat', \
                        numpy.column_stack( (self.axis, self.partdos[label]) ) )

        if self.partdosElem != {}:
            for elem in self.partdosElem:
                numpy.savetxt('dos_' + elem + '.dat', \
                     numpy.column_stack( (self.axis, self.partdosElem[elem]) ) )

    def get(self):
        return self._modes, self._freqs, self._qpts

    def qpoints(self):
        return self._qpts

    def freqs(self):
        return self._freqs

    def modes(self):
        return self._modes

    def plot(self): pass

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 2, 2009 12:14:12 PM$"
