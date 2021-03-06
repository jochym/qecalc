Quick Start
============

Installation
-------------

QECalc itself is pure Python, but it depends on `NumPy <http://numpy.scipy.org>`_,
`diffpy.Structure <http://pypi.python.org/pypi/diffpy.Structure>`_, 
and `matplotlib <http://matplotlib.sourceforge.net>`_  packages. Matplot lib is used in few examples and is not needed for QECalc to work.  
There are two ways of installing QECalc. One is easy_install and another is checking out current development version from our svn repository.


For easy_install option make sure  setuptools is installed and type::

    easy_install qecalc

or for custom location (which has to be in $PYTHONPATH) type::

    easy_install -d $INSDIR qecalc



svn option allows to get current development snapshot of QECalc. You should have svn client installed
and go through the following steps:

1. Go to your installation dir ($INSDIR), for example, ~/lib-py and type::

       svn co svn://svn@danse.us/AbInitio/espresso/qecalc

   qecalc project tree will be created

2. Add qecalc to your PYTHONPATH variable::

       export PYTHONPATH=$INSDIR/qecalc:$PYTHONPATH

  ($INSDIR = ~/lib-py in this example)

3. This module also depends on `diffpy.Structure <http://pypi.python.org/pypi/diffpy.Structure>`_  packages. Make sure  setuptools is installed and type::

    easy_install diffpy.Structure


If you already have .tar.gz source distribution, you may as well just
decompress it and add qecalc to PYTHONPATH


How to use
-----------
It is essential the user knows how to use Quantum Espresso for the basic tasks.
Excellent place to start is the `Quantum Espresso wiki <http://www.quantum-espresso.org/wiki>`_ page.
It is important to check that QE input files lead to satisfactory results
before using them in automated manner.

In order to run python scripts with Quantum Espresso, one needs to provide all
the appropriate config files (e.g. pw.x input for total energy or geometry optimization;
additionally ph.x input and dynmat.x input for single phonon; etc) as well as create a configuration file,
which specifies parallel environment of your task and
all the relevant input and output files. An example of config.ini is located in qecalc directory. All
its sections do not need to be populated, only the parameters needed for a
specific task. If some of the parameters are missing, default values will be used.
The default values are located in qetask.py, pwtask.py, etc in qecalc/qetask
folder.


Before the run, check that all the pseudopotentials from the pw config file
are available.  Also make sure Quantum Espresso is in your $PATH environment
variable. One does not have to place the python scripts which use qecalc API
into the working dir(where the input files are), but for the sake of simplicity,
this manual assumes everything is in one place.

Please, see examples directory as well as API documentation for more details.
IPython is a nice tool to play with the examples.


The following lines can be used to test QECalc's installation. Just place working pw.x
input file, called 'scf.in' into the directory with the script::

    from qecalc.pwcalc import PWCalc

    configString = """
    [pw.x]
    pwfInput: scf.in
    pwOutput: scf.out
    """
    
    pwcalc = PWCalc(configString = configString)
    pwcalc.launch()


