from qecalc.qetask.matdyntask import MatdynTask

matdyn = MatdynTask('config.ini')
matdyn.output.parse()
Pol, Omegas, qPoints = matdyn.output.property('multi phonon')
print Pol
print Omegas
#setting.matdynModes = 'blahblah.modes'
#matdyn.output.parse()
#mphonOut, Omegas, qPoints = matdyn.output.property('multi phonon', withUnits = True)


