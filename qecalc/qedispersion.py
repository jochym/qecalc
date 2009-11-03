#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

#Simple cubic symmetry points:
SCubSymPoints = {
            'Gamma' : [0.0, 0.0, 0.0],
            'Gamma2': [1.0, 1.0, 1.0],
            'X'     : [0.5, 0.0, 0.0],
            'M'     : [0.5, 0.5, 0.0],
            'R'     : [0.5, 0.5, 0.5]
}

#Body Centered cubic symmetry points:
BCCubSymPoints = {
            'Gamma' : [0.0, 0.0, 0.0],
            'Gamma2': [1.0, 1.0, 1.0],
            'H'     : [0.5, 0.5, -0.5],
            'P'     : [0.25, 0.25, 0.25],
            'N'     : [0.0, 0.0, 0.5]
}


#Face Centered cubic symmetry points:
FCCubSymPoints = {
            'Gamma' : [0.0, 0.0, 0.0],
            'Gamma2': [1.0, 1.0, 1.0],
            'X'     : [0.5, 0.5, 0.0],
            'L'     : [0.5, 0.5, 0.5],
            'W'     : [0.5, 0.75, 0.25]
}




class QEdispersion():
    def __init__(self, structure):
        self.points = []
        self.dispersion = []
        self.path = []
        self._structure = structure

    def setValues(self, array):
        self.dispersion = array

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
        self.path = path
        self.points = []

        kPoints = []
        if self._structure.lattice.ibrav == 4:
            numPoints = [0]
            for i, ipnt in enumerate(nPoints):
                numPoints.append(ipnt + numPoints[i])
            for symbol, ipnt in zip(path, numPoints ):
                kPoints.append(hexSymPoints[symbol])
                self.points.append((symbol, ipnt))
        else:
            raise Exception("Dispersion path generator does not support this lattice symmetry")
    
        for i, kpt in enumerate(kPoints):
            kPoints[i] = list(self._structure.lattice.recipCartesian(kPoints[i]))

        startPath = copy.deepcopy(kPoints[0])

        startPath.append( 0. )  # add 0. for Quantum Espresso's last column counter
        Path = [ startPath ]

        for iPoint in range(len(kPoints)-1) :
            point  = kPoints[iPoint]
            point1 = kPoints[iPoint + 1]
            dPoint = [ float(point1[i] - point[i])/nPoints[iPoint] \
                for i in range(len(point)) ]
            dPoint.append( min_gt( [abs(dPoint[i]) for i in range(len(dPoint))], 1.0e-7) ) # 1.0e-7 to deal with ~ 0.0 points, like 1.e-17
            lastElem = len(Path) - 1
            for i in range(nPoints[iPoint]):
                Path.append( addPoints( Path[lastElem + i], dPoint ) )

        self.path = []
        self.axis = []
        for elem in Path:
            for i in range(len(elem)):
                if abs(elem[i]) < 1.e-14:
                    elem[i] = 0.0e0
            self.path.append(elem[0:3])
            self.axis.append(elem[3])

        self.path = numpy.array(self.path)
        self.axis = numpy.array(self.axis)

        return self.path, self.axis
#            format = "%1.12f  %1.12f  %1.12f  %1.5f"
#            print format % (elem[0], elem[1], elem[2], elem[3])

#        print "Number of phonons in list: ", sum(nPoints) + 1


    def printPoints(self):
        """Generates phonon dispersions. Presently requires properly configures
           preconfigured matdyn.in file as well as .fc (force constants) file"""
        for elem in self.points:
            print elem[0], elem[1]

    def plot(self):
        """Takes numpy array of values (nqpoints, nbands/nmodes)"""
        import pylab
        import matplotlib
        legStr = []
        markerStyle = ['b.', 'g.', 'r.', 'c.', 'm.', 'y.', 'k.', 'b,', 'g,', 'r,', 'c,', 'm,', 'y,', 'k,', 'bo', 'go', 'ro', 'co', 'mo', 'yo', 'ko']
#        cmap = matplotlib.colors.ListedColormap(colors, name = 'myPlotCollors', N = self.dispersion.shape[1])
        for i in range(self.dispersion.shape[1]):
            pylab.plot(self.axis, self.dispersion[:,i], markerStyle[i])#, color = colors[i])

            legStr.append(str(i))

#            pylab.ylim([0.0, 3.0])
        pylab.legend(legStr, numpoints = 1)
        for i in range(len(self.points)):
            xcoord = self.axis[ self.points[i][1] ]
            pylab.axvline(x=xcoord, color='black')
            # need to add greek labels
        pylab.show()


if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Sep 4, 2009 2:26:50 PM$"