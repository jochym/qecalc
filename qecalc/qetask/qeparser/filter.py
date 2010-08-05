# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2010  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
Filter - class for filtering card and namelist parameters in configuration input

Filtering can be positive (add/set parameters): apply(input, type='plus')
or negative (remove/discard parameters): apply(input, type='minus')
"""

from namelist import Namelist
from card import Card

CARD_KEYS       = ("name", "lines", "arg")
NAMELIST_KEYS   = ("name", "params")
CARD_REQ        = "name"
NAMELIST_REQ    = "name"
POSITIVE        = "plus"
NEGATIVE        = "minus"

# Auxiliary functions
# Ternary C operator: '?:' (e.g. 'a ? a: 4')
ifelse  = lambda a,b,c: (b,c)[not a]    #

class Filter(object):

    def __init__(self, name = None):
        """
            name:  (str) -- Name of the filter
        """
        self._name          = name
        self._fnamelists    = []    # Filtered namelists
        self._fcards        = []    # Filtered cards


    def name(self):
        "Returns name of the filter"
        return self._name


    def setParam(self, namelist, param, value=""):
        """
        Set parameter of the namelist. If the parameter or namelist do not exist,
        it creates them.

            namelist: (str) -- Name of the namelist
            param: (str) -- Name of the parameter
            value: (str) -- Valur of the parameter
        """
        for nl in self._fnamelists:
            if nl.name() == namelist:
                nl.set(param, value)
                return

        # namelist doesn't exists
        nl  = Namelist(namelist)
        nl.set(param, value)
        self._fnamelists.append(nl)


    def removeParam(self, namelist, param):
        """
        Removes parameter from the namelist. If the parameter or namelist do
        not exist it will just ignore. If the last parameters if removed from
        the namelist, namelist still remains

            namelist: (str) -- Name of the namelist
            param: (str) -- Name of the parameter
        """
        for nl in self._fnamelists:
            if nl.name() == namelist:   # No name found, just ignore
                nl.remove(param)


    def setNamelist(self, namelist):
        """
        Sets namelist as dictionaty:

            namelist: (dict) -- Dictionary of namelist

        Format:
        {"name":      <namelist name>,
         "params":    {param1:    value1,
                       param2:    value2,
                       ...}}

        Example:
        {"name":      "control",
         "params":    {calculation:   'scf',
                       restart_mode:  'from_scratch',
                       tprnfor:       .true.}}
        """
        if not namelist:    # Ignore empty card
            return

        self._checkDictFormat(namelist, NAMELIST_KEYS, NAMELIST_REQ)

        nl       = Namelist(namelist["name"])
        if namelist.has_key("params"):
            for p in namelist["params"].keys():
                value   = namelist["params"][p]
                nl.set(p, value)

        self._fnamelists.append(nl)



    def removeNamelist(self, name):
        """
        Removes namelist specified by name

            name: (str) -- name of the namelist
        """
        self._remove(name, self._fnamelists)


    def namelists(self):
        """
        Returns namelists
        """
        return self._fnamelists


    def setCard(self, card):
        """
        Set card as dictionary

            card: (dict) -- Dictionary of card

        Format:
        {"name":      <card name>,
         "lines":    (line1,
                      line2,
                       ...),
         "arg":       <argument>}

        Example:
        {"name":      "k_points",
         "lines":     ("4 4 4 1 1 1",),
         "arg":       "automatic"}
        """
        if not card:    # Ignore empty card 
            return

        self._checkDictFormat(card, CARD_KEYS, CARD_REQ)
        
        c       = Card( card["name"],
                        arg = ifelse(card.has_key("arg"), card.get("arg"), None))
        if card.has_key("lines"):
            c.setLines(card["lines"])
            
        self._fcards.append(c)


    def removeCard(self, name):
        """
        Removes card specified by name

            name: (str) -- name of the card
        """
        self._remove(name, self._fcards)
        

    def cards(self):
        """
        Returns filter cards
        """
        return self._fcards
    

    def apply(self, input, type="plus"):
        """
        Applies filter to the input

            input: (object: QEInput) -- Input object
            type: (str: 'plus' or 'minus') -- Type of operation
        """
        if type == POSITIVE:
            self._applyPositive(input)

        if type == NEGATIVE:
            self._applyNegative(input)

        # Otherwise, ignore


    def _applyPositive(self, input):
        for fnl in self._fnamelists:
            nl      = input.namelist(fnl.name())
            params  = fnl.paramlist()
            for p in params:
                nl.set(p, fnl.get(p))

        for fc in self._fcards:
            card    = input.card(fc.name())
            card.setLines(fc.lines())


    def _applyNegative(self, input):
        for fnl in self._fnamelists:
            nl      = input.namelist(fnl.name())
            params  = fnl.paramlist()
            for p in params:
                nl.remove(p)

            self._removeIfEmpty(input, nl, fnl) # If namelist (from filter or input) empty

        for fc in self._fcards:
            input.removeCard(fc.name())
            

    def _removeIfEmpty(self, input, nl, fnl):
        """
        Removes namelist from input if it is empty or filter namelist has no parameters

            input: (object: QEInput) -- Input the namelist is checked from
            nl: (object: Namelist) -- namelist object
            fnl: (object: Namelist) -- filter namelist object
        """
        assert nl.name() == fnl.name()
        # No params (None or []) in filter or input namelist, remove namelist from input
        if not fnl.paramlist() or not nl.paramlist():  
            input.removeNamelist(nl.name())


    def _checkDictFormat(self, item, keys=None, reqkeys=None ):
        """
        Checks dictionary format

            item: (dict) -- Dictionary which format is being checked
            keys: (list) -- Keys allowed in the dictionary
            reqkeys: (str) -- Required keys in the dictionary, single string for now
        """
        if not type(item) == dict:
            raise TypeError("Parameter '%s' must be dictionary" % str(item))

        for k in item.keys():
            if not k in keys:
                raise KeyError("Invalid key: %s" % k)

        if not item.has_key(reqkeys):
            raise KeyError("No key: '%s'", reqkeys)


    def _remove(self, name, list):
        """
        Removes element from list specified by name

            name: (str) -- name of the element
        """
        if not name:    # No name, just ignore
            return

        for e in list:
            if e.name() == name:  # Remove the first element which matches name
                list.remove(e)


__date__ = "$Jul 28, 2010 2:35:31 PM$"


