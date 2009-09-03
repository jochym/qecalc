import os
from qephon import QEPhon

a_range = [ 5.7285141671438824, 5.7317877776644526, 5.7350667579521337, 5.7381421125337928, 5.7414536572621246, 5.7445384926803174, 5.7477444904967463, 5.7514802489352217, 5.755333850422792, 5.7589075115531303, 5.7629454999299998, 5.7666501568080868]
c_a_range = [2.2090980266298565, 2.2097169512764472, 2.2103232326527733, 2.2111584911247455, 2.2117142863526338, 2.2125256870194541, 2.213190687501422, 2.2141090078689722, 2.214882104904242, 2.2159692795744084, 2.2165114564527064, 2.2174293230773712]

qe = QEPhon('config.ini')
indexRange = range(5,6,2)
print indexRange
for i in indexRange:
    qe.structure.lattice.a = a_range[i]
    qe.structure.lattice.c = c_a_range[i]*a_range[i]
    qe.structure.lattice.printBase()
    qe.structure.saveStructureToPWSCF()
    qe.cleanOutDir()
    qe.multiPhononTaskLauncher()

    cpMatdynModesCmdStr = "cp "+"matdyn.modes "+ "matdyn" + str((i+1)*2) + ".modes"
    os.system(cpMatdynModesCmdStr)
    backupName = 'backup'+str((i+1)*2)+'.tgz'
    backupCmdStr = "tar -zcf" + backupName + " *  --exclude='*tgz' --no-recursion"
    os.system(backupCmdStr)
