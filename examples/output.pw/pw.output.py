from qecalc.qetask.pwtask import PWTask

configString = """
[pw.x]
pwOutput: scf.out
"""

pw = PWTask(configString = configString)
pw.output.parse()
energy =pw.output.property('total energy')
print 'energy: ', energy
#print 'all properties: ', pw.output.properties()
print
print 'Another output:'
pw.setting.set('pwOutput', 'scf.full.out')
pw.output.parse()
energy =pw.output.property('total energy', withUnits = True )
stress = pw.output.property('stress') 
print 'energy: ', energy
print 'stress: ', stress
#print 'all properties: ', pw.output.properties()


