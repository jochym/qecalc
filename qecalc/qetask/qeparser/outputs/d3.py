#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Xiaoli Tang
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
from baseoutput import BaseOutput

class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        self.parsers = {
                'd3 tensor' : self.getD3Tensor,
                }


    def getD3TensorOld(self, setting):
        '''
        Obtain third order anharmonicity Aqq from QE d3.x
        '''

        file = open(setting.d3Output, 'r')
        for line in file.readlines():
            if 'Number of q in the star =' in line:
                numberq = int(line.split()[7])
        file.close()

        file = open(setting.d3fildyn,'r')
        line = file.readline()
        line = file.readline()
        line = file.readline()
        ntype = int(line.split()[0])
        natom = int(line.split()[1])
        ibrav = int(line.split()[2])
        celldm = line.split()[3:]
        symbol = [0] * ntype
        mass = [0] * ntype
        type = [0] * natom
        position = [[0] * 3] * natom
        qpoints = []
        d3tensor = []

        for i in range(0,ntype):
            line = file.readline()
            symbol[i] = (line.split()[1])
            mass[i] = line.split()[2]
        for i in range(0,natom):
            line = file.readline()
            type[i] = line.split()[1]
            position[i] = line.split()[2:]

        for q in range(0,numberq):
            line = file.readline()
            line = file.readline()
            line = file.readline()
            line = file.readline()
            qpoints.append([float(f) for f in line.split()[3:6]])
            line = file.readline()
            for ii in range(0,3 * natom):
                line = file.readline()
                line = file.readline()
                if ii is not int(line.split()[1]) - 1:
                   return SyntaxError('wrong match')
                line = file.readline()

                for j in range(0,natom):
                    for k in range(0,natom):
                        line = file.readline()
                        for jdir in range(0,3):
                            jj = (j - 1) * 3 + jdir
                            kk = (k - 1) * 3
                            line = file.readline() + file.readline()
                            _val = [float(f) for f in line.split()]
                            d3tensor.append(complex(_val[0],_val[1]))
                            d3tensor.append(complex(_val[2],_val[3]))
                            d3tensor.append(complex(_val[4],_val[5]))

        file.close()
        #print 'qpoints are',qpoints
        #print 'ntype is ', ntype
        #print 'natom is ', natom
        #print 'd3 tensor is ', d3tensor
        return  (np.rollaxis(np.array(d3tensor).reshape(\
                                 numberq,natom,3,natom,natom,3,3),2,5), None), \
                (np.array(qpoints),None),


    def getD3Tensor(self, setting):
        '''
        Obtain third order anharmonic tensor from QE d3.x
        '''

        file = open(setting.d3fildyn, 'r')
        lines = file.readlines()
        ntype = int(lines[2].split()[0])
        natom = int(lines[2].split()[1])
        d3tensor = []
        qpoints = []
        numberq = 0
        for index, line in enumerate(lines):
            if 'q = (' in line:
                numberq = numberq + 1
                qpoints.append([float(f) for f in lines[index].split()[3:6]])
                for mode in range(0,3*natom):
                    for i in range(0,natom):
                        for j in range(0,natom):
                            for idir in range(0,3):
                                _itemp = index+mode*(natom*natom*7 + 3)+\
                                         (i*natom+j)*7 + 6 + idir*natom
                                temp = lines[_itemp]+lines[_itemp + 1]
                                _val = [float(f) for f in temp.split()]
                                d3tensor.append(complex(_val[0],_val[1]))
                                d3tensor.append(complex(_val[2],_val[3]))
                                d3tensor.append(complex(_val[4],_val[5]))

        return  \
        (np.rollaxis(np.array(d3tensor).reshape(\
        numberq,natom,3,natom,natom,3,3),2,5), None), (np.array(qpoints), None),

if __name__ == "__main__":
    print "Hello World";
    qpoints, d3tensor = parse_d3('si.anh_X','si.d3X.out')
    print 'qpoints are', qpoints
    for q in range(0,len(qpoints)):
          print 'd3tensors of qpoint', qpoints[q], 'is:'
          print d3tensor[q]
          print ' '

__author__="Xiaoli Tang"
__date__ ="$Nov 28, 2009 4:57:29 PM$"