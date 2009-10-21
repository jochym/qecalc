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

# Regular expressions
COMMENT         = '!.*'                 # Comment
NAME            = '([a-zA-Z_]*)[^/]'    # Extracts namelist name ()
SPACES          = '[ \t]*'              # Spaces and tabs
NO_SPACES       = '[^\s]*'              # No spaces
NEWLINE         = '[\n\r]*'             # New line ()
PARAMTER        = '[\w,()]+'            # Parameter characters (space is not allowed)
VALUE           = '[^\s,]+'             # Parameter's value (numerate all possible characters)
EXPRESSION      = '(%s%s=%s%s)' % (PARAMTER, SPACES, SPACES, VALUE)     # Parameter's expression
NAMELIST        = """%s&%s%s([^&]*)/""" % (SPACES, SPACES, NAME)        # Namelist block (handles directory slashes)
OPEN_BRACKET    = '[({]?'               # Open bracket
CLOSE_BRACKET   = '[)}]?'               # Close bracket
CARD            = '(%s[\w]+)%s%s(%s[\w]*%s)%s' % (SPACES, SPACES, OPEN_BRACKET, SPACES, SPACES, CLOSE_BRACKET)  # Card name
EMPTY_LINE  = r'^\s*'                # Empty line

import re
from orderedDict import OrderedDict
from namelist import Namelist
from card import Card

class QEParser:
    def __init__(self, filename=None, configText=None, type='pw'):
        self.namelists  = OrderedDict()
        self.cards      = OrderedDict()
        self.filename   = filename
        self.configText = configText
        self.namelistRef    = None
        self.cardRef        = None
        self.type           = type

    def parse(self):
        self._setReferences()

        if self.configText is not None: # First try use configText
            text = self.configText
        elif self.filename is not None: # ... then from file
            text = self._getText(self.filename)
        else:
            raise NameError('Dude, set config text or filename')  # Compalain
        
        self._parseNamelists(text)
        self._parseCards(text)
        return (self.namelists, self.cards)

    def toString(self):
        for n in self.namelists.keys():
            print self.namelists[n].toString()

        for c in self.cards.keys():
            print self.cards[c].toString()


    def _setReferences(self):
        input   = "input%s" % self.type
        module = __import__("inputs.%s" % input, globals(), \
                            locals(), [''], -1)        
        #module  = _import("inputs.%s" % input)
        self.namelistRef   = getattr(module, "namelists")
        self.cardRef       = getattr(module, "cards")
        

    def _parseNamelists(self, text):
        namelists  = OrderedDict()
        p   = re.compile(COMMENT)
        s1  = re.sub(p, '', text)           # Remove comments
        p2  = re.compile(NAMELIST)
        matches     = p2.findall(s1)        # Finds all namelist blocks
        for m in matches:
            name    = m[0].lower()
            if name in self.namelistRef:
                params  = self._parseParams(m[1])     # Parse parameters from a namelist block
                namelists[name.lower()] = params

        self._convertNamelists(namelists)

    # Converts from dictionary to Namelist
    def _convertNamelists(self, namelists):
        for name in namelists.keys():
            nl      = Namelist(name)
            for p in namelists[name]:
                nl.add(p[0], p[1])
                
            self.namelists[name] = nl

    # Parses parameters
    def _parseParams(self, text):
        params  = []
        p   = re.compile(EXPRESSION)        # Match expression
        matches = p.findall(text)
        for m in matches:
            pl  = self._getParams(m)         # Parameters list
            params.append(pl)

        return params

    def _getParams(self, text):
        """ Takes string like 'a = 2' and returns tuple ('a', 2) """
        s = text.split('=')
        for i in range(len(s)):
            s[i] = s[i].strip()

        param   = s[0]
        val     = s[1]
        # Assume that there are two values only: (variable, value) pair
        assert len(s) == 2

        return (param, val)

    def _parseCards(self, text):
        p   = re.compile(COMMENT)
        s1  = re.sub(p, '', text)       # Remove comments
        p2  = re.compile(NAMELIST)
        s2  = re.sub(p2, '', s1)        # Remove namelists
        rawlist = []

        p   = re.compile(EMPTY_LINE)
        s   = s2.split('\n')
        for line in s:
            line    = line.strip()
            if line != '':
                rawlist.append(line)

        self._convertCards(self._getCards(rawlist))

    def _getCards(self, rawlist):
        cards       = OrderedDict()
        cardName    = None
        for l in rawlist:
            p   = re.compile(CARD)
            m   = p.match(l)
            if m is not None:
                firstPart   = m.group(1).lower()
                secondPart  = m.group(2).strip().lower()    # Catch argument of the card
                if firstPart in self.cardRef:
                    cardName    = firstPart
                    cards[cardName]    = {}

                    if (secondPart != ''):
                        cards[cardName]["args"] = secondPart
                    else:
                        cards[cardName]["args"] = None
                    cards[cardName]["values"]   = []

                elif cardName is not None:
                    cards[cardName]["values"].append(l)

        return cards

    def _convertCards(self, cards):
        for cname in cards.keys():
            c   = Card(cname)
            c.setArg(cards[cname]["args"])
            for l in cards[cname]["values"]:
                c.addLine(l)

            self.cards[cname]    = c


    def _getText(self, filename):
        try:
            f = open(filename)
        except IOError:
            print "I/O error"
        except:
            import sys
            print "Unexpected error:", sys.exc_info()[0]
            raise

        text       = f.read()
        f.close()

        return text

def _import(package):
    return __import__(package, {}, {}, [''])



textPW = """
 &control
    calculation='scf'
    restart_mode='from_scratch',
    tprnfor = .true.
    prefix='ni',
    pseudo_dir = '',
    outdir=''
 /
 &system
    ibrav=2,
    celldm(1) =6.65,
    nat=  1,
    ntyp= 1,
    nspin=2,
    starting_magnetization(1)=0.5,
    degauss=0.02,
    smearing='gauss',
    occupations='smearing',
    ecutwfc =27.0
    ecutrho =300.0
 /
 &electrons
    conv_thr =  1.0d-8
    mixing_beta = 0.7
 /


ATOMIC_SPECIES
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF

ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

textCards = """
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

textProblem = """
&control
calculation='scf'
! restart_mode='from_scratch',
wf_collect = .true.,
tstress = .true. ,
tprnfor = .true. ,
verbosity = 'high',
prefix='mgb2',
pseudo_dir = '/home/markovsk/projects/pslib/espresso/mgb2/',
lkpoint_dir = .false. ,
outdir='temp/'
/
&system
ibrav=4,
celldm(1) = 5.78739785,
celldm(2) = 5.78739785,
celldm(3) = 1.135794331,
! celldm(1) = 5.8260,
! celldm(2) = 5.8260,
! celldm(3) = 1.1420,
nat = 3,
ntyp = 2,
nspin = 1,
nbnd = 12,
occupations='smearing',
! degauss=0.025
degauss=0.025,
smearing = 'methfessel-paxton' ,
ecutwfc =32.0,
ecutrho =256.0,
la2f = .false.
/
&electrons
conv_thr = 1.0d-12
diago_full_acc=.TRUE.
/
ATOMIC_SPECIES
Mg 24.305 mg_6.ncpp
B 11.000 B.pbe-n-van_ak.UPF
!b_rc_1.4_pcc.ncpp
! B 10.811 B.pbe-tmnc.UPF


ATOMIC_POSITIONS alat
Mg 0.000000000 0.0000000000000000 0.000000000
B 0.500000000 0.2886751345948129 0.5678971655
B 0.000000000 0.5773502691896257 0.5678971655


K_POINTS AUTOMATIC
24 24 24 0 0 0
"""

if __name__ == "__main__":
    qeparserText    = QEParser(configText = textProblem)
    qeparserText.parse()
    qeparserText.toString()
    qeparserFile    = QEParser(filename = "../tests/ni.scf.in")
    qeparserFile.parse()
    qeparserFile.toString()


__date__ = "$Oct 9, 2009 4:34:28 PM$"


