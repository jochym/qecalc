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
NEWLINE         = '[\n\r]*'             # New line (for Windows), not used at this point
PARAMETER       = '[\w,()]+'            # Parameter characters (space is not allowed)
VALUE           = '[^\s,]+'             # Parameter's value (numerate all possible characters)
EXPRESSION      = '(%s%s=%s%s)' % (PARAMETER, SPACES, SPACES, VALUE)     # Parameter's expression
NLHEADER        = '%s&%s%s'  % (SPACES, SPACES, NAME) # Namelist header
NAMELIST        = '%s([^&]*)/' % NLHEADER  # Namelist block (handles directory slashes)
NLSCOPE         = '(%s)' % NAMELIST     # Namelist scope
OPEN_BRACKET    = '[({]?'               # Open bracket
CLOSE_BRACKET   = '[)}]?'               # Close bracket
CARD            = '(%s[\w]+)%s%s(%s[\w]*%s)%s' % (SPACES, SPACES, OPEN_BRACKET, SPACES, SPACES, CLOSE_BRACKET)  # Card name
EMPTY_LINE      = r'^\s*'               # Empty line

ATTACHSIM       = ['matdyn', 'ph', 'd3']      # Simulation types that have attachments

import re
from orderedDict import OrderedDict
from namelist import Namelist
from card import Card

class QEParser:
    """
    Flexible Parser for Quantum Espresso (QE) configuration files.
    
    It parses the file specified by filename or configuration string and stores
    parameters in namelists, cards and attachment data
    structures that later on can be used in parameters' manipulations
    """
    
    def __init__(self, filename = None, configText = None, type = 'pw'):
        """
        Parameters:
            filename    -- absolute or relative filename to be parsed
            configText  -- configuration text to be parsed
            type        -- type of the simulation
            os          -- operating system defined for carriage return
        """
        
        self.header     = None
        self.namelists  = OrderedDict()     # Namelist dictionary
        self.cards      = OrderedDict()     # Cards dictionary
        self.attach     = None
        self.filename   = filename
        self.configText = configText
        self.namelistRef    = None
        self.cardRef        = None
        self.type           = type


    def parse(self):
        """Parses string and returns namelists, cards, attachment and header"""
        self.setReferences()

        if self.configText is not None: # First try use configText
            text = self.configText
        elif self.filename is not None: # ... then from file
            text = self._getText(self.filename)
        else:
            raise NameError('Dude, set config text or filename')  # Compalain

        self._parseHeader(text)
        self._parseNamelists(text)
        self._parseAttach(text)
        self._parseCards(text)
        
        return (self.header, self.namelists, self.cards, self.attach)


    def toString(self):
        """Returns string of parsed values"""
        
        s   = ''
        if self.header:
            s   += self.header

        for n in self.namelists.keys():
            s   += self.namelists[n].toString()

        for c in self.cards.keys():
            s   += self.cards[c].toString()

        if self.attach:
            s   += self.attach

        return s


    def setReferences(self):
        """Sets reference names for namelists and cards for specified simulation type"""
        
        input   = "input%s" % self.type
        module  = _import("inputs.%s" % input)
        self.namelistRef   = getattr(module, "namelists")
        self.cardRef       = getattr(module, "cards")

        return (self.namelistRef, self.cardRef)
    

    def _parseHeader(self, text):
        """Cuts the first line if it header"""
        start   = self._namelistStart(text)
        if start is not None and start == 0:
            return  # There is no header 

        lines   = text.splitlines(True)
        if lines:
            self.header   = lines[0]


    def _parseNamelists(self, text):
        """Parses text and populates namelist dictionary"""

        namelists   = OrderedDict()
        s1          = self._removeComments(text)
        p2          = re.compile(NAMELIST)
        matches     = p2.findall(s1)        # Finds all namelist blocks
        
        for m in matches:
            name    = m[0].lower()
            if name in self.namelistRef:
                params  = self._parseParams(m[1])     # Parse parameters from a namelist block
                namelists[name.lower()] = params

        self._convertNamelists(namelists)


    def _removeComments(self, text):
        """Removes comments from the text"""
        p   = re.compile(COMMENT)
        s   = re.sub(p, '', text)
        return self._removeEmptyLines(s)
        

    def _convertNamelists(self, namelists):
        """Converts dictionary to Namelist"""
        for name in namelists.keys():
            nl      = Namelist(name)
            for p in namelists[name]:
                nl.add(p[0], p[1])
                
            self.namelists[name] = nl


    def _parseParams(self, text):
        """Parses parameters"""
        params  = []
        p       = re.compile(EXPRESSION)        # Match expression
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


    def _parseAttach(self, text):
        """Special case for simulations that have attachments"""

        if self.type in ATTACHSIM:
            self.attach = self._cutNamelistText(text)


    def _parseCards(self, text):
        """Parses text and populates cards dictionary"""
        if self.type in ATTACHSIM:  # There should not be cards for simulations with attachment
            return

        s   = self._cutNamelistText(text)
        self._convertCards(self._getCards(self._rawlist(s) ))


    def _cutNamelistText(self, text):
        """Cuts the namelist text"""
        
        s       = self._removeComments(text)
        end     = self._namelistEnd(s)

        if end is not None:
            s   = s[end + 1:]  # Suppose that cards and attachmet starts with new line

        return s


    def _rawlist(self, text):
        """Removes empty lines or lines with white spaces only"""
        rawlist = []

        s   = text.splitlines(True)
        for line in s:
            stripped    = line.strip()
            if stripped != '':
                rawlist.append(line)    # rawlist contains both cardnames and card lines in order

        return rawlist


    def _removeEmptyLines(self, text):
        """Joins non empty lines or lines with white spaces only in the list to one string"""
        return ''.join(self._rawlist(text))


    def _getCards(self, rawlist):
        cards       = OrderedDict()
        cardName    = None
        for l in rawlist:
            isCardName  = False
            p   = re.compile(CARD)
            m   = p.match(l)
            if m is not None:       # If card name matches
                firstPart   = m.group(1).lower()
                secondPart  = m.group(2).strip().lower()    # Catch argument of the card
                if firstPart in self.cardRef:
                    isCardName  = True
                    cardName    = firstPart
                    cards[cardName]    = {}

                    if (secondPart != ''):
                        cards[cardName]["args"] = secondPart
                    else:
                        cards[cardName]["args"] = None
                    cards[cardName]["values"]   = []

            if cardName is not None and not isCardName:
                cards[cardName]["values"].append(l)

        return cards


    def _convertCards(self, cards):
        for cname in cards.keys():
            c   = Card(cname)
            c.setArg(cards[cname]["args"])
            for l in cards[cname]["values"]:
                c.addLine(l)

            self.cards[cname]    = c


    def _namelistStart(self, text):
        """Returns the start character position of the first namelist in the text"""
        s           = self._removeComments(text)
        p           = re.compile(NAMELIST)
        matches     = p.finditer(s)        # Finds all namelist blocks
        starts      = []

        for m in matches:
            starts.append(m.start())

        if len(starts):
            return starts[0]    # Get first value

        return None


    def _namelistEnd(self, text):
        """
        Returns the end character position of the last namelist in the text
        Notes:
            - text should be clear from comments (filtered by _removeComments(text)).
              Otherwise the end of the last namelist will be incorrect
        """
        s           = self._removeComments(text)

        p           = re.compile(NAMELIST)
        matches     = p.finditer(s)        # Finds all namelist blocks
        ends      = []

        for m in matches:
            ends.append(m.end())

        size    = len(ends)
        if size:
            return ends[size-1]    # Get last position value

        return None


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
    return __import__(package, globals(), locals(), [''], -1)



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
CELL_PARAMETERS
   0.993162743  -0.000000000   0.000000000
  -0.496581371   0.860104165  -0.000000000
  -0.000000000  -0.000000000   4.345938530
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

textMatdyn = """
 &input
    asr='crystal',
    amass(1)=24.305, amass(2)=11.000,
    flfrc='mgb2666.fc'
 /
176
0.000000    0.000000    0.456392    0.000000
0.000000    0.000000    0.447264    0.009128
0.000000    0.000000    0.438137    0.018256
0.000000    0.000000    0.429009    0.027384
0.000000    0.000000    0.419881    0.036511
"""

textDynmat = """
&input  fildyn='mgb2.dynG', asr='simple',
        q(1)=0.0, q(2)=0.0, q(3)=0.0 /
"""

textPh  = """
 &inputph
  tr2_ph=1.0d-10,
  amass(1)=24.305,
  amass(2)=11.000,
  prefix='mgb2',
  outdir='/scratch/markovsk/mgb2'
  fildyn='mgb2.dynG',
 /

"""

textHeader  = """
&INPUTPH
   tr2_ph = 1.0d-12,
   prefix = 'si',
   epsil = .false.,
   trans = .true.,
   zue = .false.,
   outdir = '/scratch/si',
   amass(1) = 28.0855,
   fildyn = 'si.dyn_G',
   fildrho = 'si.drho_G',
/
0.0 0.0 0.0
"""

# This is not a problem text (just add spaces between commas)
textComma   = """&input
   asr='crystal',  dos=.true.
   amass(1)=26.982538, amass(2)=11.000,
   flfrc='mgalb4666.fc', fldos='mgalb4.666.phdos', nk1=28,nk2=28,nk3=28
/
"""

def testMatdyn():
    parser    = QEParser(configText = textMatdyn, type="matdyn")
    parser.parse()
    print parser.toString()


def testDynmat():
    parser    = QEParser(configText = textDynmat, type="dynmat")
    parser.parse()
    print parser.toString()


def testFile():
    parser    = QEParser(filename = "../tests/ni.scf.in")
    parser.parse()
    print parser.toString()

def testCards():
    parser    = QEParser(configText = textCards)
    parser.parse()
    print parser.toString()


def testComma():
    parser          = QEParser(configText = textComma, type="matdyn")
    parser.parse()
    print parser.toString()


def testHeader():
    parser          = QEParser(configText = textHeader, type="ph")
    parser.parse()
    print parser.toString()

def testMgB2():
    parser          = QEParser(filename = "../tests/ph.mgb2.in", type="ph")
    parser.parse()
    print parser.toString()


if __name__ == "__main__":
    pass
    #testMatdyn()
    #testDynmat()
    #testFile()
    #testCards()
    #testComma()
    #testHeader()
    #testMgB2()

__date__ = "$Oct 9, 2009 4:34:28 PM$"




