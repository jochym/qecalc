from qecalc.qetask.matdyntask import MatdynTask
from qecalc.qetask.pwtask import PWTask
from qeutils import kmesh
from qeutils.bandstobxsf import bandstobxsf


configString = """
[pw.x]
pwInput:  scf.in
pwOutput: scf.out

[matdyn.x]
matdynInput:   matdyn.in
matdynOutput:  matdyn.out
"""

matdyn = MatdynTask(configString = configString)
pw = PWTask(configString = configString)

pw.input.parse()
matdyn.input.parse()

qmesh = [17,17,17]

qpoints = kmesh.kMeshCart(qmesh,pw.input.structure.lattice.reciprocalBase())

matdyn.input.qpoints.set(qpoints)
matdyn.input.save()

matdyn.launch()

#matdyn.syncSetting()
#matdyn.output.parse()

modes, freqs, qpts = matdyn.output.property('multi phonon')

bands = kmesh.packBands(qmesh, freqs)
phononEnergy = 720.0
bandstobxsf(phononEnergy, pw.input.structure.lattice.reciprocalBase(), bands,\
'phonons.bxsf')

