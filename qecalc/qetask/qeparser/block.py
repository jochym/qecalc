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

class Block(object):
    
    def name(self):
        "Return name of the namelist"
        return self._name


    def setName(self, name):
        "Set name in lower case"
        self._name = name.lower()


    def toString(self):
        return "Block.toString()"


__date__ = "$Jul 29, 2010 6:04:03 PM$"


