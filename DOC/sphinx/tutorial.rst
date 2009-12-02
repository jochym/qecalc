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
Its primary goal is to use its classes to streamline user's work flow,
offer new functionality and provide the machinery  to build new  features using
numpy, scipy, and matplotlib. One such example can be the class Converger from
qecalc/converger.py which can be  used to converge such
properties as 'total energy', 'geometry', and 'single phonon' with respect to
any iteratable variable of PW config file. More examples can be seen in examples
directory. Sources can be checked out following the installation instructions
or can be looked up here <http://dev.danse.us/trac/AbInitio/browser/espresso/qecalc`_

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
It is essential the user knows how to use Quantum Espresso for the basic tasks.
Excellent place to start is the `Quantum Espresso wiki <http://www.quantum-espresso.org/wiki>`_ page.
It is important to check that QE input files lead to satisfactory results
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


Before the run, check that all the pseudopotentials from the pw config file
are available and your output dir existsts (e.g. temp/ ). Also make sure
Quantum Espresso is in your $PATH environment variable.

Execute your python script which uses qecalc API from your working dir.

See examples directory as well as API documentation for more details

Examples
------------

PWCalc
^^^^^^^

PWCalc consists of one task, launching pw.x. Before running the example, one needs
to create config.ini file in the current dir as well as scf.in input file for pw.x.
Example of config.ini is provided below::

    [Launcher]
    # parallelization parameters
    # if this section is empty - serial mode is used
    paraPrefix:   mpiexec -n 8
    paraPostfix: -npool 8

    #useTorque: True
    #paraPrefix: mpirun --mca btl openib,sm,self
    #paraPostfix: -npool 900

    # this string will be passed to qsub, -d workingDir -V are already there:
    #torqueResourceList: -l nodes=16:ppn=12 -N MyJobName -j oe


    [pw.x]
    # pwscf input/output files
    pwscfInput:  scf.in
    pwscfOutput: scf.out


lookupProperty() goes through the all hte  output files of given qalc::

    # PWCalc
    from qecalc.pwcalc import PWCalc
    pwcalc = PWCalc('config.ini')
    pwcalc.launch()
    pwcalc.lookupProperty('total energy')
    pwcalc.lookupProperty('total energy', withUnits = True)
    pwcalc.lookupProperty('stress', withUnits = True)
    pwcalc.lookupProperty('forces', withUnits = True)


MultiPhononCalc
^^^^^^^^^^^^^^^^

config.ini, pw.x, ph.x, q2r.x, and matdyn.x input files should be in the
current dir. config.ini should have additional sections corresponding to
additional tasks::

    [ph.x]
    #ph.x input/ouput, relevant to all phonon calculations:
    phInput:  ph.in
    phOutput: ph.out


    [dynmat.x]
    #dynmat.x input/output files relevant to single phonon calculation
    dynmatInput:  dynmat.in
    dynmatOutput: dynmat.out


    [q2r.x]
    # input/output files relevant to multiple phonon calculation
    q2rInput:      q2r.in
    q2rOutput:     q2r.out


    [matdyn.x]
    # input/output files relevant to multiple phonon calculation
    matdynInput:   matdyn.in
    matdynOutput:  matdyn.out
    matdynModes:   matdyn.modes
    matdynFreqs:   matdyn.freq
    matdynfldos:   matdyn.phdos

In the following example it is also assumed outputs are laready there
after a succesfull run::

    from qecalc.multiphononcalc import MultiPhononCalc
    mphon = MultiPhononCalc('config.ini')
    for task in mphon.taskList:
        task.output.parse()
    mphon.lookupProperty('total energy', withUnits = True)
    # this will output out qpoints, frequencies and eigen modes
    mphon.lookupProperty('multi phonon', withUnits = True)
    mphon.dispersion.launch('M', 'Gamma', 'A','L', 50, 50, 50)
    mphon.dispersion.plot()
    
Converger
^^^^^^^^^^^

Class converger will converge a value  with respect to k-points or different parameters in 'system'
namelist of pw.x input file. Currently, the value can be 'total energy',
'fermi energy' or 'single phonon'::

    from qecalc.converger import Converger
    opt = Converger('config.ini','total energy', tolerance = 0.1)
    ecut = opt.converge(what = 'ecutwfc', startValue = 18, step = 4)
    conv_thr = opt.converge(what = 'conv_thr', startValue = 1e-4, multiply = 0.1)

    
