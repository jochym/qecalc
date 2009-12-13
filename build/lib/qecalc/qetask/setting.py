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
import os.path
import StringIO
class Setting:
    def __init__(self, filename = None, configString = None):
        try:
            if filename == None:
                if configString == None:
                    raise NameError("Config should be initialized either with a \
                                                     filename or configString")
                else:
                    self.filename = StringIO.StringIO(configString)
            else:
                self.filename = filename
        except NameError:
            raise        


        
    def section(self, sectionName, configDic = {}):
        """
        will parse self.configFileName values not found in the file will be
        initialized from configDic
        """
        import ConfigParser

        config = ConfigParser.SafeConfigParser(configDic)
        config.read(self.filename)

        if not config.has_section(sectionName):
            config.add_section(sectionName)

        for varName in configDic:
            varValue = config.get(sectionName, varName)
            setattr(self, varName, varValue)
            if 'input' in varName or 'Input' in varName:
                if os.path.isfile(varValue):
                    file = open(varValue, 'r')
                    string = file.read()
                    setattr(self, varName+'Str', string)

