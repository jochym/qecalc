# This programm will calculate generalized phonon DOS from QHA outputed
# DOS/partial DOS files
from qephon import QEPhon
import pylab
import periodictable
import numpy
import math

def Gaussian(x, mu, FWHM):
    sigma = FWHM/(2.0*math.sqrt(2.0*math.log(2.0)))
    return 1.0/(sigma*math.sqrt(math.pi*2.0))*math.exp(-(((x - mu)/sigma)**2)/2.0)
    
def convoluteGaussian(axis, data, FWHM):
    conv = numpy.zeros(data.shape[0])
    for i, omega0 in enumerate(axis):
        for j, omega in enumerate(axis):
            conv[i] = conv[i] + data[j]*Gaussian(omega0 - omega, 0.0 , FWHM)
    return conv
#            print conv[i]
def isNaN(x):
     return (x == x) == False

def read_table(fname):
    pdos = []
    file = open(fname,'r')
    line = file.readline()
    while line:
        if '#' not in line:
            # also convert to eV
#            s = [float(w) for w in line.split()]
#            ss = [w for w in line.split()]
#            print s
            pdos.append([float(w)*0.1239 for w in line.split()])
        line = file.readline()
    return numpy.array(pdos)

mp = QEPhon('config.ini')
pDOS = []
for atom in mp.structure.diffpy():
#for atom in mp.structure.diffpy()
    fname = 'projected_DOS.'+ atom.element
    pdos = read_table(fname)
    for i in range(pdos.shape[0]):
        for j in range(pdos.shape[1]):
            if isNaN( pdos[i,j] ): pdos[i,j] = pdos[i-1,j] + pdos[i+1,j]
    pDOS.append(pdos)
axis = pDOS[0][:, 0]
natom = len(pDOS)

#Calculate Generlized DOS
neutronCrossSections = [3.63, 5.56, 5.56]
gDOS = numpy.zeros(pDOS[0].shape[0])
print gDOS.shape
for atom, sigma, pdos in zip(mp.structure.diffpy(), neutronCrossSections, pDOS):
    gDOS = gDOS + pdos[:,2]*sigma/mp.structure.atomicSpecies[atom.element].mass

gDOSConv = convoluteGaussian(axis, gDOS, 4.0)

str1 = ''
str2 = ''
for i, omega in enumerate(axis):
    str1 = str1 + str(omega) + '    ' + str(gDOS[i]) + '\n'
    str2 = str2 + str(omega) + '    ' + str(gDOSConv[i]) + '\n'


gDOS_file = open('gdos.dat', 'w')
gDOSConv_file = open('gdos_conv.dat', 'w')

gDOS_file.write(str1)
gDOSConv_file.write(str2)


# check total DOS:
totalDos = (pDOS[0][:,2] + 2.0*pDOS[1][:,2])/3.0
area = numpy.sum(totalDos[:])
area2 = numpy.sum( pDOS[0][:,1])
#pylab.plot(axis, totalDos/area)
#pylab.plot(axis, pDOS[0][:,1]/area2)
#gplot = [Gaussian(x, 20, 4) for x in axis]
pylab.plot(axis, gDOSConv)
pylab.show()
    

