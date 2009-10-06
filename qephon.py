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



from parser.qe_io_dict import *
class QEPhonQHA(QECalc):
    """This is a wrapper class for Eyvaz Isaev QHA code"""
    def __init__(self, fname):
        QECalc.__init__(self,fname)


    def dosLauncher(self):
        """What QHA (By Eyvaz Isaev) needs:
        1) set proper Vertecies
        2) Run tetra.x to generate k-points
        3) Prepare matdyn.in
        4) Run matdyn.x
        5) Run Partial_phonon DOS.x < phdos1.in Needs ttrinp, phdos1.in, matdyn.modes
           in current dir
           Will generate partial_DOS file
        6) Run phonon_dos.x < frequency. Needs phdos.in, frequency
        7) To run total DOS with new frequencies just provide new frequency file
           Partial DOSes will be wrong of course"""
        import os

        self.prepareMatdyn()
        
        os.system('./Partial_phonon DOS.x < phdos1.in')
        os.system('./phonon_dos.x < ' + self.matdynFreqs)

    def prepareMatdyn(self):
        """What QHA (By Eyvaz Isaev) needs:
        1) set proper Vertecies
        2) Run tetra.x to generate k-points
        3) Prepare matdyn.in
        4) Run matdyn.x
        """
        import os
        self.__setVertecies()
        os.system('./tetra.x')

        matdynIn = read_file(self.matdynInput)
        keyStart = find_key_from_string(matdynIn, '/')
        newDic = {}
        for i in range(keyStart+1)[1:]:
            newDic[i] = matdynIn[i]
        save_dic(newDic, self.matdynInput)

        kptsFile = open('kpts_out','r')
        kptsString = kptsFile.read()

        # write new kpoints into file
        file = open(self.matdynInput, 'a')
        file.write(kptsString)

        file.close()
        kptsFile.close()
        self.matdynLauncher()


    def dosRelauncher(self):
        """Will relaunch DOS calculation with updated frequencies (obtained e.g.
           from loadPhonons). Only total DOS will be correct. Should update
           lattice parameters before calling this routine"""
        import os
        #update vertices with new lattice parameters:
        self.__setVertecies()
        os.system('./tetra.x') # generate new q points in kpts_out
        self.__setQptsFromKptsOut()
        self.saveMatdynFreq()
        cmd = './phonon_dos.x < ' + self.matdynFreqs
        print cmd
        os.system(cmd)


    def loadPhonons(self, fname = None):
        self.__modes, self.__freqs, self.__qpts = self.getMultiPhonon(fname)


    def saveMatdynFreq(self, fname=None):
        """Will save frequencies into text file in matdyn.freq format.
           Needed for phonon_dos.x program"""
        if fname == None:
            fname = self.matdynFreqs
        file = open(fname, 'w')
        str = ' &plot nbnd=   %d, nks=%d /\n'% (self.structure.nat*3, len(self.__qpts))
        file.write(str)
        for i in range(len(self.__qpts)):
            str = '           %f  %f  %f\n'% (self.__qpts[i,0],  \
                                             self.__qpts[i,1], self.__qpts[i,2])
            file.write(str)
            str = ''
            for j in range(self.__freqs.shape[1]):
                if j % 5 == 0 and j != 0:
                    str = str + '  %*.4f\n'% (8, self.__freqs[i,j])
                else:
                    str = str + '  %*.4f'%  (8, self.__freqs[i,j])
            if (self.structure.nat*3-1) % 5 != 0:
                str = str + '\n'
            file.write(str)
        file.close()
            
                   

    def __load(self):

        pDOS = []
        for atom in self.structure.diffpy():
            #for atom in mp.structure.diffpy()
            fname = 'projected_DOS.'+ atom.element
            pdos = self.__read_table(fname)
            for i in range(pdos.shape[0]):
                for j in range(pdos.shape[1]):
                    if self.__isNaN( pdos[i,j] ): pdos[i,j] = pdos[i-1,j] + pdos[i+1,j]
            pDOS.append(pdos)
            axis = pDOS[0][:, 0]
            return axis, numpy.array(pDOS)


    def DOS(self):
        axis, pDOS = self.__load()
        return axis, pDOS[0,:,1]



    def __isNaN(self, x):
        return (x == x) == False

    def __setVertecies(self):
        import os
        if self.structure.lattice.ibrav != 4:
            raise Exception('This lattice type is not supported')

        if self.structure.lattice.ibrav == 4:
            c_a_2 = self.structure.lattice.c/self.structure.lattice.a/2.0
            os.system("sed 's/X/'" + str(c_a_2) + "'/g' ttrinp_hcp > ttrinp")

    def __read_table(self,fname):
        pdos = []
        file = open(fname,'r')
        line = file.readline()
        while line:
            if '#' not in line:
                pdos.append([float(w) for w in line.split()])
            line = file.readline()
        return numpy.array(pdos)

    def __setQptsFromKptsOut(self):
        kpts = []
        file = open('kpts_out', 'r')
        line = file.readline()
        while line:
            words = line.split()
            if len(words) > 1:
                kpts.append([float(w) for w in words])
            line = file.readline()
        self.__qpts = numpy.array(kpts)

    def setFreqs(self, freqs):
        self.__freqs = freqs


