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

from qedispersion import QEDispersion

class PHDispersion(QEDispersion):
    def __init__(self, pwInput, matdynTask):
        QEDispersion.__init__(self,pwInput.structure)
        self.matdynTask = matdynTask

    # temporary hack
    def setPhononPath(self, *pathNPoints):
        from qetask.qeparser.output.parser.qe_io_dict import *
        matdynIn = read_file(self.matdynTask.setting.matdynInput)
        keyStart = find_key_from_string(matdynIn, '/')
        newDic = {}
        for i in range(keyStart+1)[1:]:
            newDic[i] = matdynIn[i]
        save_dic(newDic, self.matdynTask.setting.matdynInput)


        pwInput.parse()
        self.setPath(*pathNPoints)

        file = open(self.matdynTask.setting.matdynInput, 'a')
        file.write(self.__toMatdynString())
        file.close()
        self.matdynTask.launch()
        pol, disp, qpoints =  self.matdynTask.output.property('multi phonon')
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
