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
class Setting:
    def __init__(self, fname):        
        try:
          if fname == None:
             raise NameError("Config should be initialized with a filename")
        except NameError:
            raise

        self.filename = fname


        
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
            setattr(self, varName, config.get(sectionName, varName))
