QECalc overview
-----------------

`Quantum Espresso <http://www.quantum-espresso.org>`_ (QE) is widely used GNU distributed open source ab initio package
for plane wave Density Functional Theory (DFT) and molecular dynamics calculations.
Often users need to go beyond basic capabilities of an ab initio program and
use its outputs for more advanced tasks. Some examples:

* Convergence studies of a property of interest with respect to ranging values of different input parameters
* Various optimization and minimization problems
* Plotting and data processing

QECalc is a set of Quantum Espresso launchers and input/ouput parsers
organized  under a single API.
Its primary goal is to use its classes to streamline user's work flow,
offer new functionality and provide the machinery  to build new  features using
numpy, scipy, and matplotlib. 

QECalc follows the modular design paradigm of Quantum
Espresso by wrapping its executables into 'QETask' objects. QECalc offers
a library of such QETasks for several QE modules (e.g. pw.x, ph.x etc) and defines
several Calcs based on lists of these Tasks. One convenient example
is the class Converger in qeutils directory  which can be  used to converge such
properties as 'total energy' and 'single phonon' with respect to
many iterable variables of pw.x config file.

Please, read INSTALL file for installation instructions and/or go to the project's
documentation at http://docs.danse.us/AbInitio/espresso/qecalc

You may also check for new QECalc version at

    http://pypi.python.org/pypi/qecalc