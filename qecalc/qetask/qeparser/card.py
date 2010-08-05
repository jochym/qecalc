#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
Card - class that corresponds to Card in QE config file
"""

# TODO: Add some more convenience methods?

from block import Block

class Card(Block):

    def __init__(self, name, arg = None):
        self._arg      = arg
        self._name     = name.lower() # keeps lower name
        self._lines    = []


    def arg(self):
        "Returns card argument"
        return self._arg


    def setArg(self, arg):
        """
        Sets card argument

            arg: (str) -- Card argument
        """
        if arg is None:
            self._arg   = None
            return

        self._arg    = arg.lower()


    def line(self, num):
        """
        Returns line in the card with order number num
        
            num: (int) -- Order number
        """
        self._checkRange(num)
        return self._lines[num]


    def lines(self):
        "Returns list of lines"
        return self._lines


    def addLine(self, line):
        """
        Append line to the card

            line: (str) -- Line appended to the end of the card
        """
        line    = line.strip()
        self._lines.append(line)


    def setLines(self, lines):
        """
        Sets card value (list of lines) to a new one

            lines: (list) -- List of lines
        """
        self._lines    = lines


    def removeLine(self, num):
        """
        Removes line with order number num

        Parameters:
            num:        int
                Order number of the line's list
        """
        self._checkRange(num)
        self._lines.pop(num)


    def removeLines(self):
        "Removes cards value"
        self._lines = []


    def _checkRange(self, num):
        """
        Check if num is in range of lines list

        Parameters:
            num:        int
                Order number of the line's list
        """
        assert num >= 0
        assert len(self._lines) > num


    def toString(self, indent=1, br="\n"):
        """
        Dumps card as a sting

        Parameters:
            indent:     int
                Number of spaces in indent for parameters
            br:         str
                Separator between parameters
        """
        ind  = ""
        for i in range(indent):    # Populate indent
            ind += " "

        s = self.name().upper()
        if self._arg is not None:
            s += ' (%s)%s' % (self._arg, br)
        else:
            s += br

        for l in self._lines:
            s += '%s%s%s' % (ind, l.strip().strip('\n').strip(br), br)

        return s


    # DEPRICATED: Use setLines() instead
    def editLines(self, lines):
        self.setLines(lines)


__date__ = "$Aug 27, 2009 7:34:32 AM$"




