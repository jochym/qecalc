from qecalc import QECalc
import numpy

class QEPhon(QECalc):
    def __init(self, fname)__:
        QECalc.__init__(fname)
        self.__freqs = None
        self.__modes = None
        self.__qpts = None

    def launch(self):
        self.multiPhononLauncher()
        self.load()
        
    def load(self, fname = None):
        self.__modes, self.__freqs, self.__qpts = self.getMultiPhonon(fname)

    def getall(self):
        return self.__modes, self.__freqs, self.__qpts

    def qpoints(self):
        return self.__qpts

    def freqs(self):
        return self.__freqs

    def modes(self):
        return self.__modes

    def setRange(self, minOmega, maxOmega, deltaOmega):
        if minOmega == None:
            minOmega = numpy.min(self.__freqs)
        if maxOmega == None:
            maxOmega = numpy.max(self.__freqs)
        if minOmega > maxOmega: minOmega, maxOmega = maxOmega, minOmega
        if deltaOmega == None:
            deltaOmega = (maxOmega - minOmega)/200.0
        return minOmega, maxOmega, deltaOmega


    def DOS(self, minOmega = None, maxOmega = None, deltaOmega = None):
        minOmega, maxOmega, deltaOmega =  \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        histOmega = numpy.zeros(int((maxOmega-minOmega)/deltaOmega))
        norm = 0.0
        for cell_freqs in self.__freqs:
            for omega in cell_freqs:
                idx = int( (abs(omega) - minOmega)/deltaOmega )
                if idx < len(histOmega):
                    histOmega[idx] = histOmega[idx] + 1.0
                    norm = norm + 1.0
        return histOmega/norm

    def partDOS(iAtom, minOmega = None, maxOmega = None, deltaOmega = None):
        minOmega, maxOmega, deltaOmega =  \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        histPartOmega = numpy.zeros(int((maxOmega-minOmega)/deltaOmega))
        norm = 0.0
        for cell_freqs, vectors in zip(self.__freqs, self.__modes):
            for omega, vector in zip(cell_freqs, vectors[:,iAtom,:]):
                idx = int( (abs(omega) - minOmega)/deltaOmega )
                if idx < len(histPartOmega):
                    weight = (real(vector*vector.conjugate())).sum()
                    histPartOmega[idx] = histPartOmega[idx] + weight
                    norm = norm + weight
        return histPartOmega/norm


