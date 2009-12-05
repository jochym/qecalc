from qecalc.qetask.matdyntask import MatdynTask
from qecalc.qetask.pwtask import PWTask
from qeutils import kmesh
from qeutils.bandstobxsf import bandstobxsf


matdyn = MatdynTask('config.ini')
pw = PWTask('config.ini')

pw.input.parse()
matdyn.input.parse()

qmesh = [17,17,17]

qpoints = kmesh.kMeshCart(qmesh,pw.input.structure.lattice.reciprocalBase())

matdyn.input.qpoints.set(qpoints)
matdyn.input.save()

matdyn.launch()

#matdyn.output.parse()

modes, freqs, qpts = matdyn.output.property('multi phonon')

bands = kmesh.packBands(qmesh, freqs)
phononEnergy = 720.0
bandstobxsf(phononEnergy, pw.input.structure.lattice.reciprocalBase(), bands,\
'phonons.bxsf')

