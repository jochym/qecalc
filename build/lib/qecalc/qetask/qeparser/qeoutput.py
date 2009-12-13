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

class QEOutput(object):
    def __init__(self, setting, type='pw'):
        self.setting   = setting
        self.type       = type
        self._properties = None
        outModule = __import__("outputs." + self.type, globals(), \
                                locals(), ['Output'], -1)
        self.output = outModule.Output()

    def listParsers(self):
        parserList = []
        for parserName in self.output.parsers:
            parserList.append(parserName)
        return parserList

    def parse(self, parserList = 'all'):
        if parserList == 'all':
            parserList = self.listParsers()
        properties = {}
        for parserName in parserList:
            try:
                properties[parserName] = self.output.parse(parserName, self.setting)
            except KeyError: pass
            except IOError: pass
            except TypeError: pass
            except IndexError: pass
            except ValueError:
                properties[parserName] = [(None, None)]
        self._properties = properties

    def property(self, name, withUnits = False):
        """Will output the property of interest as a set of tuples if withUnits is True.
           Will output property itself otherwise"""
        prop = []
        #print self._properties[name]
        #if len(self._properties[name]) > 1:
        #    print 'preved'
        #    properties = self._properties[name]
        #else:
        #    properties = self._properties[name]
        #print self._properties[name]
        if name not in self._properties:
            if withUnits:
                return (None, None)
            else:
                return None
        for tpl in self._properties[name]:
            #print tpl
            if withUnits:
                prop.append(tpl)
            else:
                prop.append(tpl[0])
               # print "hello",  tpl[0]
        #print prop[0]
        if len(self._properties[name]) > 1:
            return prop
        else:
            return prop[0]

    def properties(self, withUnits = False):
        propDic = {}
        for propertyName in self._properties:
            propDic[propertyName] = self.property(propertyName, withUnits)
        return propDic

def test():
    from setting import Setting
    qeOut = QEOutput(Setting('config.ini'), 'pw')
    print qeOut.listParsers()
    qeOut.parse()
    print qeOut.properties()
    qeOut.parse(['total energy', 'stress'])
    print qeOut.properties()
    print qeOut.property('total energy')
    print qeOut.property('stress')
    qeOut = QEOutput(Setting('config.ini'), 'matdyn')
    print qeOut.listParsers()
    qeOut.parse()
    print qeOut.properties()

if __name__ == "__main__":
    test()
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 18, 2009 7:51:00 PM$"