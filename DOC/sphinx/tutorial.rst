Tutorial
========

Introduction
------------
`Quantum Espresso <http://www.quantum-espresso.org>`_ (QE) is widely used GNU distributed open source ab initio package
for plane wave Density Functional Theory (DFT) and molecular dynamics calculations.
Often users need to go beyond basic capabilities of an ab initio program and
use its outputs for more advanced tasks. Some examples:

* convergence studies of a property of interest with respect to ranging values of different input parameters
* Various optimization and minimization problems
* Plotting and data processing

QECalc is a set of Quantum Espresso launchers and input/ouput parsers available
organized  under single API.
Its primary goal is to use its classes to streamline users workflow,
offer new functionality and provide the machinery  to build new  features using
numpy, scipy, and matplotlib. One such example can be the class Converger from
qecalc/converger.py wich can be  used to converge such
properties as 'total energy', 'geometry', and 'single phonon' with respect to
any iteratable variable of PW config file. More examples can be seen in examples
directory.

Installation
------------
Currently there is no easy_install option for qecalc.
You should have svn client installed
and go through the following steps:

1. Go to your installation dir ($INSDIR), for example, ~/apps and type:
   svn co svn://svn@danse.us/AbInitio/espresso/qecalc
   qecalc project tree will be created

2. Add qecalc to your PYTHONPATH variable:
   export PYTHONPATH=$INSDIR/qecalc:$PYTHONPATH

3. This module also depends on `diffpy.Structure <http://pypi.python.org/pypi/diffpy.Structure>`_  package. Make sure  setuptools is installed and type:
   easy_install diffpy.Structure

Usage
------------
It is essential the user knows how to use Qunatum Espresso for the basic tasks.
Excellent place to start is the `Quantum Espresso wiki <http://www.quantum-espresso.org/wiki>`_ page.
It is important to check that QE input files lead to statisfactory results
before using them in automated manner.

In order to run python scripts with Quantum Espresso, one needs to provide all
the appropriate config files (scf.in for total energy or geometry optimization;
additionally ph.in and dynmat.in for single phonon; or ph.in, q2r.in and matdyn.in
for multi phonon calculation; etc) and place config.ini
into working dir, which specifies parallel environment of your task as well as
all the relevant input and output files. An example of config.ini is located in qecalc directory. All
its sections do not need to be populated, only the parameters needed for a
specific task. If some of the parameters are missing, default values will be used.
The default values are located in qecalc/qecalc/settings.py


Before the run, check all the pseudopotentials from the pw config file
are available and your output dir existsts (e.g. temp/ ). Make sure
Quantum Espresso is in your $PATH environment variable. See examples directory
for more details

Execute your python script which uses qecalc API from your working dir.

