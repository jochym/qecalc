#!/usr/bin/env python

# Parses configuration file (for pw.x) and stores it in dictionary.
# It also allows to dump the dictionary to create the configuration file
# The main use case is to parse the existing configuration file, change 
# some values and save it back to the configuration file.

# Input parameters are defined in INPUT_PW.html 

"""
Stability issues:
- Card starts with card title on a separate line and values between card titles.
- Prints both Namelists and Cards in capital 
- Refactoring?  Introduce class relation: Namelist(Block), Card(Block)
"""

from orderedDict import OrderedDict
from namelist import Namelist
from card import Card
from qeparser import QEParser

#Type of the configuration file can be:
#type =
#    'pw'               - (default)
#    'ph'               -
#    'pp'               -
#    'bands'            -
#    'cppp'             -
#    'd3'               -
#    'dos'              -
#    'dynmat'           -
#    'initial_state'    -
#    'gipaw'            -
#    'd1'               -
#    'matdyn'           -
#    'projwfc'          -
#    'pwcond'           -
#    'q2r'              -

class QEInput(object):
    """Quantum Espresso configuration class. It can:
    - Parse existing configuration file
    - Add, Edit or Remove parameters from/to namelist or card
    """

    # Either filename or config (not both) can be specified
    def __init__(self, filename=None, config=None, type='pw'):
        self.header     = None
        self.filename   = filename
        self.config     = config
        self.parser     = QEParser(filename, config, type)
        self.type       = type
        self.namelists  = OrderedDict()
        self.cards      = OrderedDict()
        self.attach     = None          # Specific for 'matdyn'
        self.namelistRef    = None
        self.cardRef        = None
        self.qe         = [self.header, self.namelists, self.cards, self.attach]

    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary """
        (self.header, self.namelists, self.cards, self.attach) = self.parser.parse()

    def createNamelist(self, name):
        """Creates namelist and adds to QEInput. """
        nl  = Namelist(name)
        self.namelists[name] = nl

    def addNamelist(self, namelist):
        """Adds namelist. """
        self.namelists[namelist.name()] = namelist

    def removeNamelist(self, name):
        try:
            del(self.namelists[name])
        except KeyError:    # parameter is not present
            return

    def namelist(self, name):
        # Do I need editNamelist()?
        try:
            return self.namelists[name]
        except KeyError:    # parameter is not present
            raise

    def createCard(self, name):
        """Creates card and adds to QEInput. """
        self.cards[name] = Card(name)

    def addCard(self, card):
        """Adds card. """
        self.cards[card.name()] = card

    def removeCard(self, name):
        try:
            del(self.cards[name])
        except KeyError:    # parameter is not present
            return

    def card(self, name):
        # Do I need editNamelist()?
        try:
            return self.cards[name]
        except KeyError:    # parameter is not present
            raise

    def toString(self):
        (self.namelistRef, self.cardRef)    = self.parser.setReferences()
        s = ''
        if self.header:             # Add header
            s   += self.header

        for name in self.namelistRef:   # Add namelists
            nl  = self.getObject(name, self.namelists)
            if nl is not None:
                s   += nl.toString()

        for name in self.cardRef:   # Add cards
            c  = self.getObject(name, self.cards)
            if c is not None:
                s   += c.toString()

        if self.attach:             # Add attribute (e.g. for type='matdyn')
            s   += self.attach

        return s


    def getObject(self, name, dict):
        """Returns object that corresponds to 'name'"""
        for n in dict.values():
            if n.name() == name:
                return dict[name]

        return None


    def save(self, filename=None):
        """ Saves the QEInput to the configuration file"""
        default = "config.out"

        if filename is None:
            if self.filename is not None:
                filename = self.filename
            else:
                filename = default

        f = open(filename, "w")
        f.write(self.toString())
        f.close()

    def type(self):
        return self.type

        
def _import(package):
    return __import__(package, globals(), locals(), [''], -1)


# Tests
def testCreateConfig():
    print "Testing creation of config file"
    qe  = QEInput()
    nl  = Namelist('control')
    nl.add('title', "'Ni'")
    nl.add('restart_mode', "'from_scratch'")
    print "Adding parameters to namelist:\n%s" % nl.toString()
    nl.set('title', "'Fe'")
    qe.addNamelist(nl)
    print "Adding namelist to QEInput:\n%s" % qe.toString()

    c = Card('atomic_species')
    c.addLine('Ni  26.98  Ni.pbe-nd-rrkjus.UPF')
    print "Adding line to card:\n%s" % c.toString()
    qe.addCard(c)
    print "Adding card to QEInput:\n%s" % qe.toString()
    #qe.save()


def testParseConfig():
    print "Testing parsing config file"
    qe  = QEInput("../tests/ni.scf.in")
    qe.parse()
    print qe.toString()
    nl  = qe.namelist('control')
    nl.add('title', 'Ni')
    nl.remove('restart_mode')
    qe.removeCard('atomic_species')
    nl.set('calculation', "'nscf'")
    c = qe.card('atomic_positions')
    c.editLines(['Say Hi! :)'])
    print qe.toString()
    #qe.save("../tests/ni.scf.in.mod")

def testAttach():
    qe  = QEInput("../tests/si.ph.in", type="ph")
    qe.parse()
    qe.save("../tests/si.ph.in.mod")


if __name__ == "__main__":
    pass
    #testCreateConfig()
    #testParseConfig()
    #testAttach()



