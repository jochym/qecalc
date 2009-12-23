from qecalc.qetask.matdyntask import MatdynTask

configString = """
[matdyn.x]
flvec: matdyn.modes
"""

matdyn = MatdynTask(configString = configString)
matdyn.output.parse()
Pol, Omegas, qPoints = matdyn.output.property('multi phonon')
print Pol
print Omegas


