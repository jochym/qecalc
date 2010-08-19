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
Namelist class that corresponds to Namelist in QE configuration file
"""

from orderedDict import OrderedDict
from block import Block

class Namelist(Block):
    
    def __init__(self, name):
        """
            name: (str) -- Name of the namelist in lower case. Example: "control"
        """
        self._name      = name.lower()  # keeps lower name
        self._params    = OrderedDict() # Replace dictionary by ordered dictionry


    def get(self, param, quotes = True):
        """
        Returns paramater value. If no parameter exists, return None.
        When quotes=True quotes are not added to parameter's value.

            param: (str) -- Parameter of the namelist
            quotes: (bool) -- True - if add quotes '' to parameters value,
                              False - otherwise

        Note: replaces param()
        """
        if not self._paramExists(param):
            return None

        param   = param.lower()
        if quotes:
            return self._params[param]
        
        return self._unquote(self._params[param])


    def set(self, param, val, quotes = False):
        """
        Sets existing parameter to the specified value.
        If no parameter exists, create one
        
            param: (str) -- Parameter name
            val: (str) -- Parameter value
            quotes: (bool) -- Add quotes to the value or not
        """
        param   = param.lower()
        if quotes:
            val     = self._quote(val)

        self._params[param] = val


    # TODO: Rename to params() in the future
    def paramlist(self):
        """
        Returns list of parameter names
        """
        return self._params.keys()


    def remove(self, param):
        """
        Deletes parameter

            param: (str) -- Name of the parameter
        """
        if self._paramExists(param):
            del(self._params[param])


    def exists(self,param):
        """
        Checks if parameter exists in the namelist

            param: (str) -- Name of the parameter
        """
        return self._paramExists(param)


    def _quote(self, val):
        """
        Quotes value with "'" quote mark

            val: (str) -- Value to be quoted
        """
        return "'" + val.strip('"').strip("'") + "'"
    
    
    def _unquote(self, val):
        """
        Removes quotes "'" (unquotes) on both sides of the string

            val: (str) -- Value to be unquoted
        """
        return val.strip('"').strip("'")


    def toString(self, indent = 4, br = "\n"):
        """
        Dumps namelist as a sting
        
            indent: (int) -- Number of spaces in indent for parameters
            br: (str) -- Separator between parameters
        """
        ind  = ""
        for i in range(indent):    # Populate indent
            ind += " "

        s = '&%s%s' % (self.name().upper(), br)

        for p in self._params.keys():
            s += '%s%s = %s,%s' % (ind, p, self._params[p], br)

        s += "/%s" % br 
        return s


    def _paramExists(self, param):
        """
        Checks if parameter exists in self._params

            param: (str) -- Name of the parameter
        """
        try:
            param = param.lower()
            self._params[param]
            return True
        except KeyError:    # parameter is not present
            return False


    # DEPRICATED METHODS:
    # DEPRICATED: Use get() instead
    def param(self, param, quotes = True):
        return self.get(param, quotes)

    # DEPRICATED: Use set() instead!
    def add(self, param, val, quotes = False):
        self.set(param, val, quotes)

    # DEPRICATED: Use set() instead!
    def editParam(self, param, val):
        self.set(param, val)

    # DEPRICATED: Use set() instead!
    def addParam(self, param, val):
        self.add(param, val)

    # DEPRICATED: Use remove() instead!
    def removeParam(self, param):
        self.remove(param)

__date__ = "$Aug 27, 2009 7:30:39 AM$"





