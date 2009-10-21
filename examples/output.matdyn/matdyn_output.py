from qecalc.setting import Setting
from qecalc.qetask.taskmatdyn import MatdynTask

setting = Setting('config.ini')
matdyn = MatdynTask(setting)
matdyn.output.parse()
Pol, Omegas, qPoints = matdyn.output.property('multi phonon')
print Pol
print Omegas
#setting.matdynModes = 'blahblah.modes'
#matdyn.output.parse()
#mphonOut, Omegas, qPoints = matdyn.output.property('multi phonon', withUnits = True)


