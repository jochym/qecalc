from converger import Converger

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


opt = Converger('config.ini','total energy', tolerance = 0.1)

ecut = opt.converge(what = 'ecutwfc', startValue = 18, step = 4)
ecutrho = opt.converge('ecutrho', ecut*4, 16)
#opt.converge('kpoints',[12,12,12],[2,2,2])

#print opt.getForces()
#print opt.getStress()

