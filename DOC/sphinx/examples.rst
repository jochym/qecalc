Examples
------------

PWCalc
^^^^^^^

PWCalc consists of one task, launching pw.x. Before running the example, one needs
to create a config file as well as scf.in input file for pw.x.


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