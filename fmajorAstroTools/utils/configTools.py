# -*- coding: utf-8 -*-
from __future__ import (print_function, absolute_import, division)
# Copyright 2015 Wu Jin <wujin_astro@pku.edu.cn>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
''' configTools
to save and load config files
vim:set foldmethod=indent
'''
'''Update history
    2015-11-08 add license
    2015-11-09 add function fieldLoad
    2015-11-09 make fieldLoad to a class
    2015-11-11 change the settingName prefix from "#" to "@"
               change the behavior of toNumpyArray
    2016-01-27 mv str to a demo file
    2016-07-22 support multiline input for list and json
               fix a bug caused by json package
    2016-08-11 add support for custom delimiter
'''
'''Bugs or Cautions:
    1. if you input json format like {"a":123.0}
        you can NOT mark 123.0 as 123.
        or you will have a Exception
            this seems like a bug in the json package....
'''
__author__           = 'Wu Jin <wujin_astro@pku.edu.cn>'
__lastModifiedDate__ = '2016-07-22'
__date__             = __lastModifiedDate__
__version__          = '0.2'

import os
import pdb
from collections import OrderedDict
import numpy as np
import json
import warnings

def safe_mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

class DataLoadingWarning(UserWarning):
    pass

def equalLoad(filename):
    '''
        load config file into a dict, the config file must be in this format
        config file example begin
        # some commments
        # can have empty lines, like below

        # below are data
        a = 2
        b =  fasdfa
        d  =  asdf  # we can have comments after the equation
        config file example end
        data = equalLoad(filename)
        data is like
        {'a':'2', 'b':'fasdfa', 'd':'asdf', 'comments':[a list for all the comments line]}
    '''
    lineNumber = 0
    config=OrderedDict()
    config['comments'] = []
    with open(filename, 'r') as f:
        for eachLine in f:
            lineNumber += 1
            eachLine = eachLine.strip()
            if not eachLine:
                continue
            if eachLine[0] in "#":
                config['comments'].append(eachLine[1:])
            eachPiece = eachLine.split('=')
            try:
                assert len(eachPiece)==2
                config[eachPiece[0].strip()]=eachPiece[1].split('#')[0].strip()
            except:
                raise Exception('config file:{}\nline:{}  error config format!!'.format(filename, lineNumber))
    return config

# functions and class for fieldLoad class
def equalDump(config, filename):
    '''write dict to a config file that have format of {} = {}'''
    with open(filename, 'w') as f:
        keys = config.keys()
        for eachKey in keys:
            f.write("{} = {}\n".format(eachKey, config[eachKey]))

def testNameValidity(name):
    for eachSym in '%$?,':
        if eachSym in name:
            return 0
        return 1

def parseComma(eachLine):
    pieces = eachLine.split(',')
    pieces = [each.split() for each in pieces]
    return pieces

def str2rightType(value, filename="?", lineNumber="?"):
    if value[0]=="{" and value[-1]=="}":
        try:
            result = json.loads(value)
            return result
        except:
            raise Exception("parse error at file:{} Line:{}".format(filename, lineNumber))
    if value[0]=="[" and value[-1]=="]":
        members = value[1:-1].split(",")
        goodMemeber = [each for each in members if each]
        result = [str2rightType(each) for each in goodMemeber]
        return result
    tests = [int, float]
    for testFunc in tests:
        try:
            return testFunc(value)
        except ValueError:
            continue
    # No match
    value = value.strip()
    if value[0] in "'\"" and value[-1] in "'\"":
        return value[1:-1]
    else:
        return value

class myList(list):
    def __init__(self, *a, **b):
        list.__init__(self, *a, **b)
        self.dtypeListType = None
    def setDtype(self, dtypeList):
        '''example: self.setDtype(['S10']*13)'''
        if len(self.struct['type']['header']) != len(dtypeList):
            raise Exception("dtype list and struct not matched!")
        self.dtypeList = dtypeList
        if type(dtypeList[0])==tuple:
            self.dtypeListType = 0
        else:
            self.dtypeListType = 1
    def toNumpyArray(self, startIndex = None):
        if 'dtypeList' not in self.__dict__:
            raise Exception("dtype list not set!!")
        if self.dtypeListType == 1:
            dtypeArray = [(eachName, eachType)\
                        for eachName,eachType in\
                            zip(self.struct['type']['header'], self.dtypeList)]
        else:
            dtypeArray = self.dtypeList
        if startIndex is None:
            startIndex = self.dtypeListType
        return np.array(self[startIndex:], dtype=dtypeArray)

class fieldLoad:
    def __init__(self, filename, delimiter=""):
        self.delimiter = delimiter
        self.result, self.struct, self.settingNames = self.parseFile(filename)

    def __getitem__(self, index):
        return self.result[index]

    def __str__(self):
        return str(self.result)

    def parseFile(self, filename):
        '''
            load config file into a dict, try parseFile(configToolsDemo.txt) and see the results
        '''
        newSettingNameFlag       = 0
        result                   = OrderedDict()
        resultStruct             = OrderedDict()
        configFormat             = ""
        lastArrayLen             = -1
        settingName              = None
        fieldName                = 'default'
        listMultiLineValue       = ""
        dictMultiLineValue       = ""
        theDict                  = result
        theDictStruct            = resultStruct
        theDict[fieldName]       = OrderedDict()
        theDictStruct[fieldName] = OrderedDict()
        settingNameList          = []
        newSettingName           = ''

        with open(filename, 'r') as f:
            for lineNumber, eachLine in enumerate(f):
                lineNumber+=1
                eachLine = eachLine.strip()
                if not eachLine or eachLine[0]=="#":
                    continue
                if eachLine[0] in "@":
                    default = theDict.pop('default')
                    for eachKey in default:
                        theDict[eachKey] = default[eachKey]
                    #theDict.update(default)
                    structDefault = theDictStruct.pop('default')
                    for eachKey in structDefault:
                        theDictStruct[eachKey] = structDefault[eachKey]
                    #theDictStruct.update(theDictStruct['default'])
                    newSettingName = eachLine[1:].strip()
                    fieldName   = 'default'
                    # test validity of setting name
                    good = testNameValidity(newSettingName)
                    if not good:
                        raise Exception("config file:{}\nLine:{} not a good setting name".format(filename, lineNumber))
                    # test ok, now begin a new setting
                    settingName               = newSettingName
                    result[settingName]       = OrderedDict()
                    resultStruct[settingName] = OrderedDict()
                    theDict                   = result[settingName]
                    theDictStruct             = resultStruct[settingName]
                    settingNameList.append(settingName)
                    configFormat = None
                    lastArrayLen = -1
                    theDict[fieldName]       = OrderedDict()
                    theDictStruct[fieldName] = OrderedDict()
                    continue
                # remove possible comment at end of the line
                eachLine = eachLine.split('#')[0].strip()
                # begin a field
                if eachLine[0]=="[" and eachLine[-1]=="]":
                    if settingName is None:
                        theDict       = result
                        theDictStruct = resultStruct
                    fieldName=eachLine[1:-1]
                    # test validity of field name
                    good = testNameValidity(fieldName)
                    if not good:
                        raise Exception("config file:{}\nLine:{} not a good field name".format(filename, lineNumber))
                    configFormat = None
                    lastArrayLen = -1
                    theDict[fieldName]       = OrderedDict()
                    theDictStruct[fieldName] = OrderedDict()
                    continue
                # not blank or comment or [filed] defination, so it must be an item
                if listMultiLineValue: # in multiLine list mode
                    value = eachLine.split('#')[0].strip()
                    if value[-1]!="]":
                        listMultiLineValue += value
                        continue
                    else:
                        listMultiLineValue += value
                        theDict[fieldName][multiLineName] =\
                                str2rightType(listMultiLineValue, filename=filename, lineNumber=lineNumber)
                        theDictStruct[fieldName]['count'] += 1
                        listMultiLineValue=""
                        continue
                if dictMultiLineValue: # in multiLine list mode
                    value = eachLine.split('#')[0].strip()
                    if value[-1]!="}":
                        dictMultiLineValue += value
                        continue
                    else:
                        dictMultiLineValue += value
                        theDict[fieldName][multiLineName] =\
                                str2rightType(dictMultiLineValue, filename=filename, lineNumber=lineNumber)
                        theDictStruct[fieldName]['count'] += 1
                        dictMultiLineValue=""
                        continue
                #print settingName, fieldName, lineNumber
                #<== in case of equation config
                if '=' in eachLine and configFormat!="," and configFormat!=" ":# it should be a equation config
                    # test for mix config type error
                    if configFormat!='=':
                        if not configFormat:
                            configFormat='='
                            theDictStruct[fieldName]['type'] = '='
                            theDictStruct[fieldName]['count'] = 0
                        else:
                            raise Exception("config file:{}\nLine:{} mix config type in one field!!".format(filename, lineNumber))
                    # try to record the term
                    #eachPiece = eachLine.split('=')
                    firstEq = eachLine.index("=")
                    eachPiece = [eachLine[:firstEq], eachLine[firstEq+1:]]
                    value = eachPiece[1].split('#')[0].strip()
                    if not value:
                        theDict[fieldName][eachPiece[0].strip()] = ''
                    elif value[0]=="[" and value[-1]!="]":
                        listMultiLineValue = value
                        multiLineName  = eachPiece[0].strip()
                    elif value[0]=="{" and value[-1]!="}":
                        dictMultiLineValue = value
                        multiLineName  = eachPiece[0].strip()
                    else:
                        theDict[fieldName][eachPiece[0].strip()] =\
                                str2rightType(value, filename=filename, lineNumber=lineNumber)
                    theDictStruct[fieldName]['count'] += 1
                    continue
                #==> incase of equation config
                # in case of arrays
                #<== in case of arrays
                #!! it should be a array with delimiter of ,
                if (',' in eachLine and not self.delimiter) or self.delimiter==",":
                    # test for mix config type error
                    if configFormat!=',':
                        if configFormat is None:
                            configFormat=','
                            theDict[fieldName] = myList()
                            theDict[fieldName].struct = theDictStruct[fieldName]
                            theDictStruct[fieldName]['type'] = {'delimiter':','}
                            theDictStruct[fieldName]['count'] = 0
                        else:
                            raise Exception("config file:{}\nLine:{} mix config type in one field!!".format(filename, lineNumber))
                    thisArray = parseComma(eachLine)
                    # test for array length error
                    if lastArrayLen!=len(thisArray):
                        if lastArrayLen == -1:
                            lastArrayLen = len(thisArray)
                            theDictStruct[fieldName]['type']['header'] = thisArray
                        else:
                            raise Exception("config file:{}\nLine:{} array length error!!".format(filename, lineNumber))
                    # try to record the term
                    if lastArrayLen==1:
                        thisArray = thisArray[0]
                    else:
                        thisArray = tuple(thisArray)
                    theDict[fieldName].append(thisArray)
                    theDictStruct[fieldName]['count'] += 1
                    continue
                #!! costom delimiter
                if self.delimiter and self.delimiter!=" ":
                    # test for mix config type error
                    if configFormat!=self.delimiter:
                        if configFormat is None:
                            configFormat=self.delimiter
                            theDict[fieldName] = myList()
                            theDict[fieldName].struct = theDictStruct[fieldName]
                            theDictStruct[fieldName]['type'] = {'delimiter':self.delimiter}
                            theDictStruct[fieldName]['count'] = 0
                        else:
                            raise Exception("config file:{}\nLine:{} mix config type in one field!!".format(filename, lineNumber))
                    thisArray = eachLine.split(self.delimiter)
                    # test for array length error
                    if lastArrayLen!=len(thisArray):
                        if lastArrayLen == -1:
                            lastArrayLen = len(thisArray)
                            theDictStruct[fieldName]['type']['header'] = thisArray
                        else:
                            raise Exception("config file:{}\nLine:{} array length error!!".format(filename, lineNumber))
                    # try to record the term
                    if lastArrayLen==1:
                        thisArray = thisArray[0]
                    else:
                        thisArray = tuple(thisArray)
                    theDict[fieldName].append(thisArray)
                    theDictStruct[fieldName]['count'] += 1
                    continue
                #!! it now can be a 1 by N array or a array with delimiter " "
                thisArray = eachLine.split()
                # test for mix config type error
                if configFormat!=' ':
                    if not configFormat:
                        configFormat=' '
                        theDict[fieldName] = myList()
                        theDict[fieldName].struct = theDictStruct[fieldName]
                        theDictStruct[fieldName]['type'] = {'delimiter':' '}
                        theDictStruct[fieldName]['count'] = 0
                    else:
                        raise Exception("config file:{}\nLine:{} mix config type in one field!!".format(filename, lineNumber))
                # test for array length error
                if lastArrayLen!=len(thisArray):
                    if lastArrayLen == -1:
                        lastArrayLen = len(thisArray)
                        theDictStruct[fieldName]['type']['header'] = thisArray
                    else:
                        #raise Exception("config file:{}\nLine:{} array length error!!".format(filename, lineNumber))
                        warnings.warn("config file:{}\nLine:{} array length error!!".format(filename, lineNumber), DataLoadingWarning)
                # try to record the term
                if lastArrayLen==1:
                    thisArray = thisArray[0]
                else:
                    thisArray = tuple(thisArray)
                theDict[fieldName].append(thisArray)
                theDictStruct[fieldName]['count'] += 1
                #==> in case of arrays

        theDict.update(theDict['default'])
        theDict.pop('default')
        theDictStruct.update(theDictStruct['default'])
        theDictStruct.pop('default')
        return result, resultStruct, settingNameList

if __name__=="__main__":
    result = fieldLoad(os.path.join(os.getcwd(), 'configToolsDemo.txt'))
    r = result.result

