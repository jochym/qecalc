from qecalc.setting import Setting
from qecalc.qetask.taskmatdyn import MatdynTask

setting = Setting('config.ini')
matdyn = MatdynTask(setting)
matdyn.output.parse()
Pol, Omegas, qPoints = matdyn.output.property('multi phonon', withUnits = True)
print Pol
print Pol[0].shape[0]
#setting.matdynModes = 'blahblah.modes'
#matdyn.output.parse()
#mphonOut, Omegas, qPoints = matdyn.output.property('multi phonon', withUnits = True)


