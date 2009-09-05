# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="markovsk"
__date__ ="$Sep 4, 2009 2:26:50 PM$"

import numpy

#High symmetry points:
hexSymPoints = {
            'Gamma' : [0.0, 0.0, 0.0],
            'Gamma2': [1.0, 1.0, 1.0],
            'M' :     [0.0, 1.0/2.0, 0.0],
            'K' :     [1./3., 1./3., 0.0],
            'H' :     [1.0/3.0, 1.0/3.0, 1.0/2.0],
            'L' :     [0.0, 0.5, 0.5],
            'A' :     [0.0, 0.0, 1.0/2.0]
}

class Dispersion():
    def __init__(self, qe):
        self.__points = []
        self.__path = []
        self.__qe = qe

    def setPath(self, *pathNPoints):
        import copy

        def addPoints(x,y): return [ \
            x[j] + y[j] for j in range(len(x))]

        def min_gt(seq, val):
            return min([v for v in seq if v > val])

        def min_ge(seq, val):
            return min([v for v in seq if v >= val])

        def max_lt(seq, val):
            return max([v for v in seq if v < val])

        def max_le(seq, val):
            return max([v for v in seq if v <= val])

        nPoints = []
        path = []
        for elem in pathNPoints:
            if isinstance(elem,int):
                nPoints.append(elem)
            else:
                path.append(elem)
        if len(nPoints) != len(path) - 1:
            raise Exception('Dispersion: Numbers of points do not agree')
        self.__nPoints = nPoints
        self.__path = path

        kPoints = []
        if self.__qe.structure.lattice.ibrav == 4:
            numPoints = [0]
            for i, ipnt in enumerate(nPoints):
                numPoints.append(ipnt + numPoints[i])
            for symbol, ipnt in zip(path, numPoints ):
                print hexSymPoints[symbol]
                kPoints.append(hexSymPoints[symbol])
                self.__points.append((symbol, ipnt))
        else:
            raise Exception("Dispersion path generator does not support this lattice symmetry")
    
        for i, kpt in enumerate(kPoints):
            kPoints[i] = list(self.__qe.structure.lattice.recipCartesian(kPoints[i]))

        startPath = copy.deepcopy(kPoints[0])

        startPath.append( 0. )  # add 0. for Quantum Espresso's last column counter
        Path = [ startPath ]

        for iPoint in range(len(kPoints)-1) :
            point  = kPoints[iPoint]
            point1 = kPoints[iPoint + 1]
            dPoint = [ float(point1[i] - point[i])/nPoints[iPoint] \
                for i in range(len(point)) ]
            dPoint.append( min_gt( ( abs( dPoint[i] )  \
                           for i in range(len(dPoint) ) ), 0) )
            lastElem = len(Path) - 1
            for i in range(nPoints[iPoint]):
                Path.append( addPoints( Path[lastElem + i], dPoint ) )

        self.__path = []
        self.__axis = []
        for elem in Path:
            for i in range(len(elem)):
                if abs(elem[i]) < 1.e-14:
                    elem[i] = 0.0e0
            self.__path.append(elem[0:3])
            self.__axis.append(elem[3])
            format = "%1.12f  %1.12f  %1.12f  %1.5f"
            print format % (elem[0], elem[1], elem[2], elem[3])

        print "Number of phonons in list: ", sum(nPoints) + 1


    def printPoints(self):
        for elem in self.__points:
            print elem[0], elem[1]

    def matdyn(self, *pathNPoints):
        from parser.qe_io_dict import *
        matdynIn = read_file(self.__qe.matdynInput)
        keyStart = find_key_from_string(matdynIn, '/')
        newDic = {}
        for i in range(keyStart+1):
            newDic[i] = matdynIn[i]
        save_dic(newDic, self.__qe.matdynInput)
#        ikey = keyStart + 1
#        print matdynIn
#        while ikey < len(matdynIn):
#            print ikey
#            print matdynIn[ikey]
#            del matdynIn[ikey]
#        save_dic(matdynIn, self.__qe.matdynInput)

        self.setPath(*pathNPoints)

        file = open(self.__qe.matdynInput, 'a')
        file.write(self.__toMatdynString())
        self.__qe.matdynLauncher()
        pol, omega, qpoints =  self.getMultiPhonon()
        return omega


    def plot(self, valPath):
        """Takes numpy array of values (nqpoints, nbands/nmodes)"""
        import pylab
        for i in range(valPath.shape[1]):
            pylab.plot(self.__axis, self.valPath[:,i])
        for i in range(len(self.__points)):
            xcoord = self.__axis[ self.__points[i][1] ]
            pylab.axvline(x=xcoord, color='black')
            # need to add greek labels
        pylab.show()

    def __toMatdynString(self):
        string = str(len(self.__path)) + '\n'
        for elem, coord in zip(self.__path, self.__axis):
            for e in elem:
                string = string + str(e) + '    '
            string = string + str(coord) + '\n'
        return string


if __name__ == "__main__":
    print "Hello World";
