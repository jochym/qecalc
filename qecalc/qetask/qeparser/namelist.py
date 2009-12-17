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

from orderedDict import OrderedDict

class Namelist():
    """Namelist class that corresponds to Namelist in QE config file"""

    def __init__(self, name):
        self.__name = name.lower() # keeps lower name
        self.params = OrderedDict() # Replace dictionary by ordered dictionry

    def name(self):
        return self.__name

    def setName(self, name):
        self.__name = name.lower()

    def param(self, param, quotes = True):
        """Returns value of parameter 'param'"""
        if self.__paramExists(param):
            if quotes:
                return self.params[param]
            else:
                return self._unquote(self.params[param])
            return self.params[param]

        return None

    def add(self, param, val, quotes = False):
        # Replaces addParam() Add verification?
        param = param.lower()
        if quotes:
            val     = self._quote(val)

        self.params[param]  = val


    def set(self, param, val, quotes = False):
        #  Replaces editParam() and addParam(). Merge with add()?
        """Edits parameter. If it doesn't exist, it just ignores it """
        if self.__paramExists(param):
            if quotes:
                val     = self._quote(val)

            self.params[param] = val

    def remove(self, param):
        """Deletes parameter"""
        if self.__paramExists(param):
            del(self.params[param])


    def exists(self,param):
        return self.__paramExists(param)


    def _quote(self, val):
        return "'" + val.strip('"').strip("'") + "'"


    def _unquote(self, val):
        return val.strip('"').strip("'")


    def toString(self, indent="    ", br="\n"):
        # Dump namelist
        # Should I use space?
        s = '&%s%s' % (self.name().upper(), br)

        for p in self.params.keys():
            s += '%s%s = %s,%s' % (indent, p, self.params[p], br)

        s += "/%s" % br
        return s


    def __paramExists(self, param):
        try:
            param = param.lower()
            self.params[param]
            return True
        except KeyError:    # parameter is not present
            return False


    # Depricated methods:
    # Depricated
    def addParam(self, param, val):
        self.add(param, val)

    # Depricated
    def editParam(self, param, val):
        self.set(param, val)

    # Depricated
    def removeParam(self, param):
        self.remove()

__date__ = "$Aug 27, 2009 7:30:39 AM$"





