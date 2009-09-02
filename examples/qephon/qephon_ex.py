from qephon import QEPhon
import pylab
# example of MgB2 phonon DOS

mphon = QEPhon('config.ini')

mphon.structure.lattice.printBase()

mphon.loadPhonons() # will load phonons from matdyn.modes file

hist1, axis = mphon.DOS()
hist2, axis = mphon.partDOS('B')
pylab.plot(axis, hist1)
pylab.show()
pylab.plot(axis, hist2)
pylab.show()

