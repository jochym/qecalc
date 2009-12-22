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

class SubSetting:
    def __init__(self): pass

class Setting:
    def __init__(self, filename = None, configString = None):
        self._paths = SubSetting()
        try:
            if filename == None:
                if configString == None:
                    raise NameError("Config should be initialized either with a \
                                                     filename or configString")
                else:
                    self.filename = None #StringIO.StringIO(configString)
                    self.configString = configString
            else:
                self.filename = filename
                self.configString = open(filename,'r').read()
        except NameError:
            raise        

        
    def section(self, sectionName, configDic = {}):
        """
        will parse self.configFileName values not found in the file will be
        initialized from configDic
        """
        import ConfigParser
        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        
        file = StringIO.StringIO(self.configString)
        config.readfp(file)
                
        if not config.has_section(sectionName):
            config.add_section(sectionName)

        for varName in configDic.keys():
            if varName not in dir(self):
                setattr(self, varName, configDic[varName])

        for option in config.options(sectionName):
            varValue = config.get(sectionName, option)
            setattr(self, option, varValue)

        Vars = dir(self)
        for varName in Vars:
            if 'input' in varName or 'Input' in varName:
                if os.path.isfile(varName):
                    file = open(varName, 'r')
                    string = file.read()
                    setattr(self, '_'+varName+'Str', string)

    def set(self, name, value):
        setattr(self, name, value)

    def get(self, name):
        if name in dir(self) and getattr(self,name) != None:
            return getattr(self,name)
        else:
            if name in dir(self._paths):
                return getattr(self._paths, name)
            else:
                return None

    def syncAllPathsInNamelist(self, param, namelist, varName, input, defaults = None):
        """
        Syncs path attribute in namelist with setting variable varName
        if varName was not set in Setting it will be initialized from QE
        input. If it is not in QE input it will be initialized from QE default
        values
        """
        var = getattr(self, varName, None)
        if var != None:
            input.namelist(namelist).add(param, var, quotes = True)
            setattr(self._paths, varName, var)
        else:
            if input.namelist(namelist).exists(param):
                inputVar = input.namelist(namelist).param(param,  quotes = False)
                setattr(self, varName, inputVar)
                setattr(self._paths, varName, inputVar)
            else:
                setattr(self, varName, defaults[varName])
                setattr(self.setting._paths, varName, defaults[varName])

    def getAllPathsInNamelist(self, param, namelist, varName, input, defaults = None):
        """
        Retrieves all the filenames relevant to given namelist. Variables
        from class Setting override ones from QE input files. If both are
        empty, default values are used
        """
        var = getattr(self, varName, None)
        fileDict = {}
        if var != None:
            fileDict[param] = var
        else:
            if input.namelist(namelist).exists(param):
                fileDict[param] = input.namelist(namelist).param(param,  \
                                                                 quotes = False)
            else:
                fileDict[param] = defaults[varName]


    def syncPathInNamelist(self, param, namelist, varName, input, defaults = None):
        """
        Syncs path attribute in namelist with setting variable varName
        """
        var = getattr(self, varName, None)
        if var != None:
            input.namelist(namelist).add(param, var, quotes = True)
            setattr(self._paths, varName, var)
        else:
            if input.namelist(namelist).exists(param):
                inputVar = input.namelist(namelist).param(param,  quotes = False)
                setattr(self._paths, varName, inputVar)
            else:
                setattr(self._paths, varName, defaults[varName])