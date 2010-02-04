Quick Tutorial
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
properties as 'total energy' and 'single phonon' with respect to
any iteratable variable of PW config file. More examples can be seen in examples
directory. Sources can be checked out following the installation instructions.

Installation
------------

There are two ways of installing qecalc. One is easy_install and another is
checking out current development version from our svn repository.


For easy_install option make sure  setuptools is installed and type::

    easy_install qecalc

or for custom location (which has to be in $PYTHONPATH) type::

    easy_install -d $INSDIR qecalc



For svn option, you should have svn client installed
and go through the following steps:

1. Go to your installation dir ($INSDIR), for example, ~/apps and type::

       svn co svn://svn@danse.us/AbInitio/espresso/qecalc

   qecalc project tree will be created

2. Add qecalc to your PYTHONPATH variable::

       export PYTHONPATH=$INSDIR/qecalc:$PYTHONPATH

  ($INSDIR = ~/apps in this example)

3. This module also depends on `diffpy.Structure <http://pypi.python.org/pypi/diffpy.Structure>`_  package. Make sure  setuptools is installed and type::

    easy_install diffpy.Structure


If you already have .tar.gz source distribution, you may as well just
decompress it and add qecalc to PYTHONPATH

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
The default values are located in qetask.py, pwtask.py, etc in qecalc/qetask
folder.


Before the run, check that all the pseudopotentials from the pw config file
are available.  Also make sure Quantum Espresso is in your $PATH environment variable.

Execute your python script which uses qecalc API from your working dir.

Please, see examples directory as well as API documentation for more details.
IPython is a nice tool to play with the examples.

Examples
------------

PWCalc
^^^^^^^

PWCalc consists of one task, launching pw.x. Before running the example, one needs
to create config.ini file in the current dir as well as scf.in input file for pw.x.
Example of config.ini is provided below::

    # all the relevant input files must be preconfigured for specific tasks
    # before using this class

    [Launcher]
    # parallelization parameters
    # if this section is empty - serial mode is used
    #paraPrefix:   mpiexec -n 8
    #paraPostfix: -npool 8

    #serialPrefix: mpiexec -n 1
    #serialPostfix:

    # default: False
    useTorque: True
    paraPrefix: mpirun --mca btl openib,sm,self
    paraPostfix: -npool 900

    serialPrefix: mpirun
    serialPostfix:

    #Name of a script to execute a command on multiple nodes
    #relevant if outdir is not located on Parallel/Network File system.
    #Default value is empty
    #paraRemoteShell: bpsh -a

    # this string will be passed to qsub, -d workingDir -V are already there:
    paraTorqueParams: -l nodes=4:ppn=12 -N myjob -j oe
    serialTorqueParams: -l nodes=1:ppn=1 -N myjob -j oe

    outdir: temp/

    [pw.x]
    pwfInput: scf.in
    pwOutput: scf.out

[Launcher] section is common for all tasks (but each task has corresponding
variables independantly from other tasks). Some tasks are serial, some are
parallel. para/serial variables specify launching parameters for these two classes
of tasks. If serialPrefix is empty, a serial task will be launched on head node.

Task sections can also contain Quantum Espresso varialbes, corresponding to input/output
files which are usually  specified in QE config files. For example 'flvec' from matdyn.x config
file  or 'fildyn' from ph.x input file. Any file variable, specified in a section
of config.ini will override one in corresponding QE config file. If none is specified,
QECalc will try to resolve their default values internally, so it will not affect parsing.
For example, default value of 'flvec' is 'matdyn.modes' and it does not have
to be specified in matdyn config file nor in config.ini.

if 'outdir' is specified in config.ini, it will override any outdir specified
in QE config input files of tasks containing outdir field.

One does not have to specify all the sections. if a section is ommited, default
values are assumed.


lookupProperty() goes through the all the  output files of a given calc::

    # PWCalc
    from qecalc.pwcalc import PWCalc
    pwcalc = PWCalc('config.ini')
    pwcalc.launch()
    print 'looking for properties in output file ', pwcalc.pw.setting.get('pwOutput')
    pwcalc.lookupProperty('total energy')
    pwcalc.lookupProperty('total energy', withUnits = True)
    pwcalc.lookupProperty('stress', withUnits = True)
    pwcalc.lookupProperty('forces', withUnits = True)

Methods task.setting.get(varName) and task.setting.set(varName, varValue) allow
to read and modify QECalc configuration dynamically for  each task.

Config file can also be passed as a string::

    configString = """
    [pw.x]
    pwfInput: scf.in
    pwOutput: scf.out
    """
    pwcalc = PWCalc(configString = configString)
    pwcalc.launch()


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
    dynmatOutput: dyn.out


    [q2r.x]
    # input/output files relevant to multiple phonon calculation
    q2rInput:      q2r.in
    q2rOutput:     q2r.out


    [matdyn.x]
    # input/output files relevant to multiple phonon calculation
    matdynInput:   matdyn.in
    matdynOutput:  matdyn.out

In the following example it is also assumed outputs are already there
after a successful run::

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

    from qeutils.converger import Converger
    opt = Converger('config.ini','total energy', tolerance = 0.1)
    ecut = opt.converge(what = 'ecutwfc', startValue = 18, step = 4)
    conv_thr = opt.converge(what = 'conv_thr', startValue = 1e-4, multiply = 0.1)

