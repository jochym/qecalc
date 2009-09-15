from launcher import Launcher
from property import Property
from parser.configParser import *
from qestructure import QEStructure
import dispersion

import numpy
from diffpy.Structure.structure import Structure
from diffpy.Structure.lattice import Lattice


class QECalc(Launcher, Property):
    def __init__(self, fname):
        
        Launcher.__init__(self, fname)
        Property.__init__(self, fname)

        self.qeConfig = QEConfig(self.pwscfInput)
        self.qeConfig.parse()

        self.structure = QEStructure(self.pwscfInput)

        self.dispersion = dispersion.Dispersion(self)

#        self._kpts = self.getkPoints()
        self._ekincutoff = self.getEkincutoff()

#    def getNamelistParameter(self, namelist, parameter):
#        self.__qeConfig.parse()
#        return self.__qeConfig.namelist(namelist).param(parameter)

#    def setNamelistParameter(self, namelist, parameter, value):
#        self.__qeConfig.namelists[namelist].params[parameter] = str(value)
#        self.__qeConfig.save(self.pwscfInput)
#
#    def getCard(self, name):
#        self.__qeConfig.parse()
#        return self.__qeConfig.cards[name].getLines()

    def _kMesh(self, nkx, nky, nkz):
        """Will generate k-point mesh in crystal coordinates in [0,1] box
        Mesh is uniform. I.e. lattice symmetry is not taken into account."""
        j = numpy.complex(0,1)
        data = numpy.mgrid[0.0:1.0:nkx*j,0:1.0:nky*j,0:1.0:nkz*j]
        data = numpy.transpose(data.reshape(3,nkx*nky*nkz))
        return data

    def kMeshCart(self, nq1, nq2, nq3):
        """Will generate k-point mesh in cartesian coordinates
        Lattice symmetry is not taken into account."""
        kpts = self._kMesh(nq1, nq2, nq3)
        kpts_cart = numpy.zeros(kpts.shape)
        for k in range(kpts.shape[0]):
            # convert into cartesian coordinates in units of inv lattice vector a
            kpts_cart[k,:] = self.structure.lattice.recipCartesian(kpts[k,:])
#            self.structure.lattice.diffpy().cartesian(kpt)/ \
#                        self.structure.lattice.a

        return kpts_cart

    def getkPointsFromPWSCF(self):
        """ Returns array of points. Does not have to be AUTOMATIC """
        kpoints = []
        for line in self.qeConfig.cards['k_points'].getLines():
            kpoints.append([float(k) for k in line.split()])
        return numpy.array(kpoints)

    def setkPointsAutomatic(self, kpoints, shifts = None):
        if shifts == None:
            shifts = [int(k) for k in self.getkPointsFromPWSCF()[0,3:]]
        self.qeConfig.cards['k_points'].removeLines()
        self.qeConfig.cards['k_points'].setArgument('AUTOMATIC')
        string = ""
        for k in kpoints:
            string = string + str(k) + " "
        for s in shifts:
            string = string + str(s) + " "
        self.qeConfig.cards['k_points'].addLine(string)
        self.qeConfig.save(self.pwscfInput)

    def getEkincutoff(self):
        self.qeConfig.parse()
        return float(self.qeConfig.namelist('system').param('ecutwfc'))


if __name__ == '__main__':
    qe = QECalc('config.ini')
    qe.qeConfig.setNamelistParameter('system', 'ecutwfc', 44)
    qe.qeConfig.save(qe.qeConfig.filename)
    print qe.getkPoints().shape
    qe.setkPointsAutomatic(numpy.array([17, 33, 12, 0, 0, 0]))
    print 'ibrav' in qe.qeConfig.namelists['system'].params