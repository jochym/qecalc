#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Xiaoli Tang,  Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from baseoutput import BaseOutput

class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        self.parsers = {
                        #'qpoints'       : self.getQpoints,
                        'qpoints'       : self.parse_qpoints_from_dyn
                        }

    def getQpoints(self, setting):
        """
        Extract ireducible qpoints from  ph output
        """
        #read Espresso output into memory:
        file = open(setting.get('phOutput'))
        phOut = file.readlines()
        posList =  \
        [i for i,line in enumerate(phOut) if 'Dynamical matrices for (' in line]
        i = posList[-1] + 3
        qpoints = []
        while len(phOut[i].split()) == 4:
            qpoints.append([ float(w) for w in phOut[i].split()[1:]])
            i = i + 1
        return [(qpoints, None)]

    def parse_qpoints_from_dyn(self, setting):
        """parse qpoints from dyn file
        for example: si.dyn is specified in ph.input
        then ph output will have si.dyn0, si.dyn1...etc
        this parser reads in si.dyn0 first getting irreducible qpoints
        and then proceed to read all the other files getting full list of qpoints

        """
        file = open(setting.get('fildyn')+'0','r')
        lines = file.readlines()
        file.close()
        [nq1,nq2,nq3]=[int(f) for f in lines[0].split()]
        Nq_indep = int(lines[1].split()[0])
        qpoints_indep = []
        qpoints_full = []
        for i in range(0,Nq_indep):
            qpoints_indep.append([float(f) for f in lines[i+2].split()])
        for i in range(0,Nq_indep):
            file = open(setting.get('fildyn')+str(i+1),'r')
            lines = file.readlines()
            file.close()
            for index, line in enumerate(lines):
                if 'axes' in line:
                    qpoints_full.append([float(f) for f in lines[index+2].split()[3:6]])

        #print [nq1,nq2,nq3]
        #print qpoints_indep
        #print len(qpoints_full)
        return ([nq1,nq2,nq3], None), (qpoints_indep, None), (qpoints_full, None)

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Nov 3, 2009 6:53:12 PM$"
