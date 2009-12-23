from qecalc.converger import Converger

#convergence in 'total energy',  'geometry' or 'single phonon' calculations
#can be studied  with respect to the following variables:
#                'nbnd'         : 'system',
#                'degauss'      : 'system',
#                'ecutwfc'      : 'system',
#                'ecutrho'      : 'system',
#                'conv_thr'     : 'electrons',
#                'etot_conv_thr': 'control',
#                'forc_conv_thr': 'control',
#                'path_thr'     : 'ions',
#                'kpoints'      : 'k_points'


configString = """
# all the relevant input files must be preconfigured for specific tasks
# before using this class

[Launcher]
# parallelization parameters
# if this section is empty - serial mode is used
paraPrefix:   mpiexec -n 8
paraPostfix: -npool 8

serialPrefix: mpiexec -n 1
serialPostfix:

# default: False
#useTorque: True
#paraPrefix: mpirun --mca btl openib,sm,self
#paraPostfix: -npool 900

#serialPrefix: mpirun
#serialPostfix:

#Name of a script to execute a command on multiple nodes
#relevant if outdir is not located on Parallel/Network File system.
#Default value is empty
#paraRemoteShell: bpsh -a

# this string will be passed to qsub, -d workingDir -V are already there:
#paraTorqueParams: -l nodes=4:ppn=12 -N myjob -j oe
#serialTorqueParams: -l nodes=1:ppn=1 -N myjob -j oe

outdir: temp/


[pw.x]
# pw input/output files
pwfInput:  scf.in
pwOutput: scf.out


[ph.x]
#ph.x input/ouput, relevant to all phonon calculations:
phInput:  ph.in
phOutput: ph.out

[dynmat.x]
#dynmat.x input/output files relevant to single phonon calculation
dynmatInput:  dynmat.in
dynmatOutput: dyn.out
"""

opt = Converger(configString = configString, taskName = 'total energy', tolerance = 0.1)
conv_thr = opt.converge(what = 'conv_thr', startValue = 1e-4, multiply = 0.1)
ecut = opt.converge(what = 'ecutwfc', startValue = 18, step = 4)
ecutrho = opt.converge('ecutrho', ecut*4, 16)
#opt.converge('kpoints',[12,12,12],[2,2,2])

#print opt.getForces()
#print opt.getStress()

