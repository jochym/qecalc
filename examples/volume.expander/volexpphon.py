import os
from qecalc.calcmultiphonon import MultiPhononCalc

# range of parameters starting from 0% expansion
a_range = [5.7883493196, 5.7915203461101745, 5.7943730191592495, 5.797358090068105, 5.8003716998968606, 5.8033186751159649, 5.8061612088663823, 5.8094295051649869, 5.8121388187147254, 5.8152279604208692, 5.8182413994691, 5.821214436566617, 5.8240276956669499]
c_a_range = [1.13667006535, 1.137073615299389, 1.1376612947393401, 1.1381676060959949, 1.1386537405854196, 1.1391757603566126, 1.1397559358795943, 1.1400820726311993, 1.1407340161439217, 1.1411590161587388, 1.1416252954554391, 1.1421120734472419, 1.1426896187884585]

mphon = MultiPhononCalc('config.ini')
indexRange = [0,2,4]
print indexRange
for i in indexRange:
    mphon.pw.input.parse()
    mphon.pw.input.structure.lattice.a = a_range[i]
    mphon.pw.input.structure.lattice.c = c_a_range[i]*a_range[i]
    mphon.pw.input.structure.save()    
    mphon.launch()

    cpMatdynModesCmdStr = "cp "+"matdyn.modes "+ "matdyn" + str(i) + ".modes"
    os.system(cpMatdynModesCmdStr)
    backupName = 'backup'+str(i)+'.tgz'
    backupCmdStr = "tar -zcf" + backupName + " *  --exclude='*tgz' --no-recursion"
    os.system(backupCmdStr)

