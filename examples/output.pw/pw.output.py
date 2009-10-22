from qecalc.setting import Setting
from qecalc.qetask.taskpw import PWTask

setting = Setting('config.ini')
pw = PWTask(setting)
pw.output.parse()
energy =pw.output.property('total energy')
print 'energy: ', energy
print 'all properties: ', pw.output.properties()
pw.setting.pwscfOutput = 'scf.full.out'
pw.output.parse()
energy =pw.output.property('total energy', withUnits = True )
stress = pw.output.property('stress') 
print 'energy: ', energy
print 'stress: ', stress
print 'all properties: ', pw.output.properties()


