from qecalc import QECalc
import numpy

class QEPhon(QECalc):
    def __init__(self, fname):
        QECalc.__init__(self, fname)
        self.__freqs = None
        self.__modes = None
        self.__qpts = None

#    def calculatePhonons(self):
#        self.multiPhononTaskLauncher()
#        self.loadPhonons()
        
    def loadPhonons(self, fname = None):
        self.__modes, self.__freqs, self.__qpts = self.getMultiPhonon(fname)

    def getPhonons(self):
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
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histOmega = numpy.zeros(nPoints)
        norm = 0.0
        for cell_freqs in self.__freqs:
            for omega in cell_freqs:
                idx = int( (abs(omega) - minOmega)/deltaOmega )
                if idx < len(histOmega):
                    histOmega[idx] = histOmega[idx] + 1.0
                    norm = norm + 1.0

        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        return histOmega/norm, axis

    def partDOS(self, atomSymbol, minOmega = None, maxOmega = None, deltaOmega = None):
        from numpy import real
        minOmega, maxOmega, deltaOmega   =      \
                                self.setRange(minOmega, maxOmega, deltaOmega)
        nPoints = int((maxOmega - minOmega)/deltaOmega)
        histPartOmega = numpy.zeros(nPoints)
        norm = 0.0
        for iAtom, atom in enumerate(self.structure.diffpy()):
            if atomSymbol == atom.element:
                for cell_freqs, vectors in zip(self.__freqs, self.__modes):
                    for omega, vector in zip(cell_freqs, vectors[:,iAtom,:]):
                        idx = int( (abs(omega) - minOmega)/deltaOmega )
                        if idx < len(histPartOmega):
                            weight = (real(vector*vector.conjugate())).sum()
                            histPartOmega[idx] = histPartOmega[idx] + weight
                            norm = norm + weight
        axis = numpy.linspace(minOmega, maxOmega, nPoints)
        return histPartOmega/norm, axis



class QEPhonQHA(QECalc):
    def __init__(self, fname):
        QECalc.__init__(fname)


    def dosLauncher(self):
        """What QHA (By Eyvaz Isaev) needs:
        1) set proper Vertecies
        2) Run tetra.x to generate k-points
        3) Prepare matdyn.in
        4) Run matdyn.x
        5) Run Partial_phonon DOS.x < phdos1.in Needs ttrinp, phdos1.in, matdyn.modes
           in current dir
           Will generate partial_DOS file
        6) Run phonon_dos.x < frequency. Needs phdos.in, frequency"""
        import os
        from parser.qe_io_dict import *
        self.__setVertecies()
        os.system('./tetra.x')

        matdynIn = read_file(self.matdynInput)
        keyStart = find_key_from_string(matdynIn, '/')
        newDic = {}
        for i in range(keyStart+1)[1:]:
            newDic[i] = matdynIn[i]
        save_dic(newDic, self.matdynInput)

        kptsFile = file('kpts_out','r')
        kptsString = kptsFile.read()

        file = open(self.matdynInput, 'a')
        
        file.write(kptsString)

        file.close()
        kptsFile.close()

        self.matdynLauncher()


        os.system('./Partial_phonon DOS.x < phdos1.in')
        os.system('./phonon_dos.x < matdyn.freq')

    def __load(self):

        pDOS = []
        for atom in structure.diffpy():
            #for atom in mp.structure.diffpy()
            fname = 'projected_DOS.'+ atom.element
            pdos = read_table(fname)
            for i in range(pdos.shape[0]):
                for j in range(pdos.shape[1]):
                    if isNaN( pdos[i,j] ): pdos[i,j] = pdos[i-1,j] + pdos[i+1,j]
            pDOS.append(pdos)
            axis = pDOS[0][:, 0]
            return numpy.array(pDOS), axis


    def __isNaN(self, x):
        return (x == x) == False

    def __setVertecies(self):
        import os
        if structure.lattice.ibrav != 4:
            raise Exception('This lattice type is not supported')

        if structure.lattice.ibrav == 4:
            c_a_2 = structure.lattice.c/structure.lattice.a/2.0
            os.system("sed 's/X/'" + str(c_a_2) + "'/g' ttrinp_hcp > ttrinp")

    def __read_table(self,fname):
        pdos = []
        file = open(fname,'r')
        line = file.readline()
        while line:
            if '#' not in line:
                # also convert to eV
                pdos.append([float(w)*0.1239 for w in line.split()])
            line = file.readline()
        return numpy.array(pdos)
