from launcher import Launcher
from property import Property
from parser.configParser import *
from qestructure import QEStructure

import numpy
from diffpy.Structure.structure import Structure
from diffpy.Structure.lattice import Lattice


class QECalc(Launcher, Property):
    def __init__(self, fname):
        
        Launcher.__init__(self, fname)
        Property.__init__(self, fname)

        print self.pwscfInput
        self.qeConfig = QEConfig(self.pwscfInput)
        self.qeConfig.parse()

        self.structure = QEStructure(self.pwscfInput)

#        self._kpts = self.getkPoints()
        self._ekincutoff = self.getEkincutoff()
        print self._ekincutoff

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

    def getkPointsFromPWSCF(self):
        """ Returns array of points. Does not have to be AUTOMATIC """
        kpoints = []
        for line in self.qeConfig.getCardLines('k_points'):
            kpoints.append([float(k) for k in line.split()])
        return numpy.array(kpoints)

    def setkPointsAutomatic(self, kpoints):
        self.qeConfig.cards['k_points'].removeLines()
        self.qeConfig.cards['k_points'].setArgument('AUTOMATIC')
        string = ""
        for k in kpoints:
            string = string + str(k) + " "
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