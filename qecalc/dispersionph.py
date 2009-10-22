#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class PHDispersion(QEDispersion):
    def __init__(self, structure):
        QEDispersion.__init__(structure)

    def setPhononPath(self,matdynTask, *pathNPoints):
        from parser.qe_io_dict import *
        matdynIn = read_file(matdynTask.setting.matdynInput)
        keyStart = find_key_from_string(matdynIn, '/')
        newDic = {}
        for i in range(keyStart+1)[1:]:
            newDic[i] = matdynIn[i]
        save_dic(newDic, matdynTask.setting.matdynInput)

        self.setPath(*pathNPoints)

        file = open(matdynTask.setting.matdynInput, 'a')
        file.write(self.__toMatdynString())
        file.close()
        matdynTask.launch()
        pol, disp, qpoints =  matdynTask.output.property('multi phonon')
        self._dispersion = disp

    def __toMatdynString(self):
        string = str(len(self.__path)) + '\n'
        for elem, coord in zip(self.__path, self.__axis):
            string = string + \
                   "%f    %f    %f    %f\n" % (elem[0], elem[1], elem[2], coord)
        return string

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 19, 2009 1:18:20 PM$"
