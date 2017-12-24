# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function, unicode_literals)
# Copyright 2016 Wu Jin <wujin_astro@pku.edu.cn>
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
''' Notes
    vim:set foldmethod=indent
    use vimrc config 'customFolding.vim' from https://gist.github.com/Fmajor to
        make reading the code much more easier!
'''
''' Update history
    2016-08-08 add the license
'''
''' Bugs or Cautions:
    None
'''

import logging
import warnings
import random
from colorama import Fore, Back, Style
import sys
import traceback
import time
import os
import pdb

initLoggerMacro="""
global logger
auxLogger = configs.get("logger", None)
if auxLogger is not None:
    if isinstance(auxLogger, dict):
        logger = messager.Messager(**auxLogger)
    else:
        logger = auxLogger
else:
    try:
        logger = globalLogger
    except NameError as e:
        logger = messager.Messager(name="tempLogger",
                         logFile="log_default.log",
                         warningFile="warning_default.log",
                         errorFile="error_default.log")
        message = "init default log in function: {}".format(sys._getframe().f_back.f_code.co_name)
        logger.p(message)
"""

#region: color print functions
Fcolor = {
    "r": Fore.RED,
    "k": Fore.BLACK,
    "b": Fore.BLUE,
    "c": Fore.CYAN,
    "w": Fore.WHITE,
    "y": Fore.YELLOW,
    "m": Fore.MAGENTA,
    "g": Fore.GREEN,
    "lr": Fore.LIGHTRED_EX,
    "lk": Fore.LIGHTBLACK_EX,
    "lb": Fore.LIGHTBLUE_EX,
    "lc": Fore.LIGHTCYAN_EX,
    "lw": Fore.LIGHTWHITE_EX,
    "ly": Fore.LIGHTYELLOW_EX,
    "lm": Fore.LIGHTMAGENTA_EX,
    "lg": Fore.LIGHTGREEN_EX,
    " ": ""
}
Bcolor = {
    "r": Back.RED,
    "k": Back.BLACK,
    "b": Back.BLUE,
    "c": Back.CYAN,
    "w": Back.WHITE,
    "y": Back.YELLOW,
    "m": Back.MAGENTA,
    "g": Back.GREEN,
    "lr": Back.LIGHTRED_EX,
    "lk": Back.LIGHTBLACK_EX,
    "lb": Back.LIGHTBLUE_EX,
    "lc": Back.LIGHTCYAN_EX,
    "lw": Back.LIGHTWHITE_EX,
    "ly": Back.LIGHTYELLOW_EX,
    "lm": Back.LIGHTMAGENTA_EX,
    "lg": Back.LIGHTGREEN_EX,
    " ": ""
}
Scolor = {
    "b": Style.BRIGHT,
    "d": Style.DIM,
    "n": Style.NORMAL,
    "r": Style.RESET_ALL,
    " ": ""
}

def cPrint(message, colorStr="   "):
    assert len(colorStr)>=3
    index = 0
    #!! processs fore
    if colorStr[index]=="l":
        fore = Fcolor[colorStr[index:index+2]]
        index+=2
    else:
        fore = Fcolor[colorStr[index:index+1]]
        index+=1
    #!! processs back
    if colorStr[index]=="l":
        back = Bcolor[colorStr[index:index+2]]
        index+=2
    else:
        back = Bcolor[colorStr[index:index+1]]
        index+=1
    #!! processs style
    style = Scolor[colorStr[index:index+1]]
    color = "".join([fore, back, style])
    coloredMessage = "{}{}{}".format(color, message, Style.RESET_ALL)
    print(coloredMessage)

#endregion: color print functions

#region: messagers
stderr = sys.stderr
stdout = sys.stdout

class TeeStdout(object):
    def __init__(self, logFile=None,
                 renew=False, bufferSize=1):
        if logFile is not None:
            if renew:
                self.f = open(logFile, "w", bufferSize)
            else:
                self.f = open(logFile, "a", bufferSize)
            if renew and os.path.exists(logFile):
                open(logFile, "w").close()
        else:
            now = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())
            logger = Messager(logFile="log{time}.txt".format(time=now),
                    logFormat=logFormats["raw"])

        self.stdout = stdout
        sys.stdout = self
        self.count=0

    def __del__(self, *arg, **kwarg):
        sys.stdout = self.stdout
        if not self.f.closed:
            self.f.close()

    def write(self, data):
        self.f.write(data)
        self.stdout.write(data)

    def flush(self):
        self.f.flush()
        self.stdout.flush()

    def fileno(self):
        return self.stdout.fileno()

    def isatty(self):
        return self.stdout.isatty()

    def __exit__(self, _type, _value, _traceback):
        sys.stdout = self.stdout
        if not self.f.closed:
            self.f.close()

    def __enter__(self):
        pass

class TeeStderr(object):
    def __init__(self, logFile=None,
                 renew=False, bufferSize=1):

        if logFile is None:
            now = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())
            logFile="log{time}.txt".format(time=now)

        if renew:
            self.f = open(logFile, "w", bufferSize)
        else:
            self.f = open(logFile, "a", bufferSize)
        if renew and os.path.exists(logFile):
            open(logFile, "w").close()

        self.stderr = stderr
        sys.stderr = self

    def __del__(self, *arg, **kwarg):
        try:
            sys.stderr = self.stderr
        except Exception as e:
            pass
            #print(e)
        if not self.f.closed:
            self.f.close()

    def write(self, data):
        self.f.write(data)
        self.stderr.write(data)

    def flush(self):
        self.f.flush()
        self.stderr.flush()

    def fileno(self):
        return self.stderr.fileno()

    def isatty(self):
        return self.stderr.isatty()

    def __exit__(self, _type, _value, _traceback):
        sys.stderr = self.stderr
        if not self.f.closed:
            self.f.close()

    def __enter__(self):
        pass

class fmajorWarning(UserWarning):
    pass
warnings.simplefilter("always", fmajorWarning)

class fmajorError(Exception):
    pass

def formatWarning(message, category, filename, lineno, line):
    message = str(message)
    message = message.split("FROM:LOGGER:")
    if len(message)==1:
        message = message[0]
        logger = ""
    elif len(message)==2:
        logger = "from logger: " + message[1]
        message = message[0]
    else:
        raise Exception("you should not have FROM:LOGGER: in the warning message")

    coloredMessage = "{fore}{back}{style}{message}{reset}".format(
        fore=Fore.YELLOW, back=Back.BLACK, style=Style.BRIGHT,
        reset=Style.RESET_ALL, message=message)

    #return "{filename}:{lineno}😱 😱 😱\n\t{logger}\n\t{message}\n".format(filename=filename,
            #lineno=lineno, message=coloredMessage, logger=logger).encode("utf8")
    return "{filename}:{lineno}😱 😱 😱\n\t{logger}\n\t{message}\n".format(filename=filename,
            lineno=lineno, message=coloredMessage, logger=logger)

def myShowWarning(message, category, filename, lineno, file=None, line=None):
    if file is None:
        file = stderr
        if file is None:
            # sys.stderr is None - warnings get lost
            return
    try:
        file.write(formatWarning(message, category, filename, lineno, line))
    except IOError:
        pass # the file (probably stderr) is invalid - this warning gets lost.

showwarning_orig = warnings.showwarning

def getEnv(stacklevel=2):
    """copy from the python std lib, just want to learn how to get the linenum and filename"""
    # Get context information
    try:
        caller = sys._getframe(stacklevel)
    except ValueError:
        globals = sys.__dict__
        lineno = 1
    else:
        globals = caller.f_globals
        lineno = caller.f_lineno
    if '__name__' in globals:
        module = globals['__name__']
    else:
        module = "<string>"
    filename = globals.get('__file__')
    if filename:
        fnl = filename.lower()
        if fnl.endswith((".pyc", ".pyo")):
            filename = filename[:-1]
    else:
        if module == "__main__":
            try:
                filename = sys.argv[0]
            except AttributeError:
                filename = '__main__'
        if not filename:
            filename = module
    return (filename, lineno)

logFormats = {
        "externalWarning": ["%(asctime)s %(process)d: %(message)s"]*3,
        "raw": ["%(message)s"]*3
}

class DefaultMessager(object):
    """my logger"""
    def __init__(self):
        pass
    def p(self, message):
        print(message)
    def print(self, message):
        print(message)
    def debug(self, message):
        print(message)
    def warn(self, message, warningClass=None):
        coloredMessage = "{fore}{back}{style}{message}{reset}".format(
            fore=Fore.YELLOW, back=Back.BLACK, style=Style.BRIGHT,
            reset=Style.RESET_ALL, message=message)
        if warningClass is None:
            warnings.warn(coloredMessage, self.warningClass)
        else:
            warnings.warn(coloredMessage, warningClass)
    def error(self, message, exceptionClass=None):
        coloredMessage = "\n\t{fore}{back}{style}{message}{reset}".format(
            fore=Fore.RED, back=Back.BLACK, style=Style.BRIGHT,
            reset=Style.RESET_ALL, message=message)
        if exceptionClass is None:
            exceptionClass = Exception
        return exceptionClass(coloredMessage)

    def log(self, message):
        messageNew = "can not log this message, because logger have no log file:\n"+message
        coloredMessage = "{fore}{back}{style}{message}{reset}".format(
            fore=Fore.YELLOW, back=Back.BLACK, style=Style.BRIGHT,
            reset=Style.RESET_ALL, message=message)
        warnings.warn(coloredMessage, self.warningClass)

    def warn2file(self, message):
        messageNew = "can not log this message, because logger have no log file:\n"+message
        coloredMessage = "{fore}{back}{style}{message}{reset}".format(
            fore=Fore.YELLOW, back=Back.BLACK, style=Style.BRIGHT,
            reset=Style.RESET_ALL, message=message)
        warnings.warn(coloredMessage, self.warningClass)

    def error2file(self, message):
        messageNew = "can not log this message, because logger have no log file:\n"+message
        coloredMessage = "{fore}{back}{style}{message}{reset}".format(
            fore=Fore.YELLOW, back=Back.BLACK, style=Style.BRIGHT,
            reset=Style.RESET_ALL, message=message)
        warnings.warn(coloredMessage, self.warningClass)

class Messager(object):
    """my logger"""
    def __init__(self, name=None, logFile=None,
                 warningFile=None, errorFile=None,
                 initLog=False, initWarningLog=False,
                 initErrorLog=False, logFormat=None,
                 renew=False, warningClass=None,
                 exceptionClass=None,
                 renewLog=None, renewWarning=None, renewError=None):
        if name is not None:
            self.loggerName = name
        else:
            now = time.strftime("%Y-%m-%dT%H:%M:%S"+str(random.random())[2:], time.gmtime())
            self.loggerName = "Messager-" + now

        if renewLog is None:
            renewLog = renew
        if renewWarning is None:
            renewWarning = renew
        if renewError is None:
            renewError = renew

        if logFormat is not None:
            formatters = logFormat
        else:
            formatters = ["%(asctime)s %(process)d %(message)s"] +\
                ["%(asctime)s %(process)d %(levelname)8s %(funcName)s\n\t%(message)s"]*2
        # add logger
        if logFile is not None:
            self.logger = logging.getLogger(self.loggerName)
            self.logger.setLevel(logging.DEBUG)
            formatter = logging.Formatter(formatters[0])
            formatter.datefmt = "%Y-%m-%d %H:%M:%S"
            if isinstance(logFile, list):
                for each in logFile:
                    fh = logging.FileHandler(each)
                    fh.setFormatter(formatter)
                    self.logger.addHandler(fh)
            else:
                fh = logging.FileHandler(logFile)
                fh.setFormatter(formatter)
                self.logger.addHandler(fh)
            self.haveLogFile = True
            if isinstance(renewLog, list):
                for each, eachRenew in zip(logFile, renewLog):
                    if eachRenew:
                        open(each,"w").close()
            else:
                if renewLog:
                    open(logFile,"w").close()
            if initLog:
                self.logger.critical("=========init logger {}===========".format(self.loggerName))
        else:
            self.haveLogFile = False
        # add logger for Warning file
        if warningFile is not None:
            self.warningLogger = logging.getLogger(self.loggerName+".warning")
            self.warningLogger.setLevel(logging.DEBUG)
            formatter = logging.Formatter(formatters[1])
            formatter.datefmt = "%Y-%m-%d %H:%M:%S"
            fh = logging.FileHandler(warningFile)
            fh.setFormatter(formatter)
            self.warningLogger.addHandler(fh)
            self.haveWarningLogFile = True
            if renewWarning:
                open(warningFile,"w").close()
            if initWarningLog:
                self.warningLogger.critical("=========init warningLogger {}===========".format(self.loggerName))
        else:
            self.haveWarningLogFile = False
        # add logger for Error file
        if errorFile is not None:
            self.errorLogger = logging.getLogger(self.loggerName+".warning.error")
            self.errorLogger.setLevel(logging.DEBUG)
            formatter = logging.Formatter(formatters[2])
            formatter.datefmt = "%Y-%m-%d %H:%M:%S"
            fh = logging.FileHandler(errorFile)
            fh.setFormatter(formatter)
            self.errorLogger.addHandler(fh)
            self.haveErrorLogFile = True
            if renewError:
                open(errorFile,"w").close()
            if initErrorLog:
                self.errorLogger.critical("=========init errorLoggedr {}===========".format(self.loggerName))
        else:
            self.haveErrorLogFile = False

        self.Exception = self.error
        if warningClass is None:
            self.warningClass = fmajorWarning
        else:
            self.warningClass = warningClass

        if exceptionClass is None:
            self.exceptionClass = fmajorError
        else:
            self.exceptionClass = exceptionClass

    def p(self, message, colorStr="   "):
        cPrint(message, colorStr)
        if self.haveLogFile:
            self.logger.info(message)
    def print(self, message, colorStr="   "):
        cPrint(message, colorStr)
    def log(self, message):
        if self.haveLogFile:
            self.logger.info(message)
    def debug(self, message):
        print(message)
        if self.haveLogFile:
            self.logger.debug(message)
    def warn(self, message, warningClass=None):
        """warnings are colored in woarnings.showwarinigs"""
        filename, lineno = getEnv()
        showwarning_old = warnings.showwarning
        warnings.showwarning = myShowWarning
        if warningClass is None:
            warnings.warn(message+"FROM:LOGGER:"+self.loggerName, self.warningClass, stacklevel=2)
        else:
            warnings.warn(message+"FROM:LOGGER:"+self.loggerName, warningClass, stacklevel=2)
        warnings.showwarning = showwarning_old
        if self.haveWarningLogFile:
            logMessage = "{}:{} from logger: {}\n\t {}".format(filename, lineno, self.loggerName, message)
            self.warningLogger.warn(logMessage)
    def warn2file(self, message):
        filename, lineno = getEnv()
        if self.haveWarningLogFile:
            logMessage = "{}:{} from logger: {}\n\t{}".format(filename, lineno, self.loggerName, message)
            self.warningLogger.warn(logMessage)
    def error(self, message, exceptionClass=None):
        message += "\n\t\t from logger: {}".format(self.loggerName)
        coloredMessage = "\n\t{fore}{back}{style}{message}{reset}".format(
            fore=Fore.RED, back=Back.BLACK, style=Style.BRIGHT,
            reset=Style.RESET_ALL, message=message)
        if exceptionClass is None:
            exceptionClass = self.exceptionClass
        if self.haveErrorLogFile:
            try:
                raise exceptionClass(message)
            except exceptionClass as e:
                logMessage = str(traceback.format_exc())
                logMessageListAux = logMessage.split("\n")
                logMessageList = [each for each in logMessageListAux if each]
                logMessage = "\n\t".join(logMessageList)
                self.errorLogger.error(logMessage)
        return exceptionClass(coloredMessage)
    def error2file(self, message, exceptionClass=None):
        message += "\n\t\t from logger: {}".format(self.loggerName)
        if exceptionClass is None:
            exceptionClass = self.exceptionClass
        else:
            exceptionClass = exceptionClass
        if self.haveErrorLogFile:
            try:
                raise exceptionClass(message)
            except exceptionClass as e:
                logMessage = str(traceback.format_exc())
                self.errorLogger.error(logMessage)
#endregion: messagers

if __name__=="__main__":
    l1=Messager(name=None, logFile="test/log1.log", warningFile="test/warnings.log", errorFile="test/errors.log", initErrorLog=True, renew=True)
    l2=Messager(name=None, logFile="test/log1.log", warningFile="test/warnings.log", errorFile="test/errors.log", initWarningLog=False, initErrorLog=False, renew=True)
    l3=Messager(name="l3:logger", logFile="test/log1.log", warningFile="test/warnings.log", errorFile="test/errors.log", initWarningLog=False, initErrorLog=False, renew=True)
    l4=Messager(name="l4:logger", logFile="test/log1.log", warningFile="test/warnings.log", errorFile="test/errors.log", initWarningLog=False, initErrorLog=False, renew=True)

    l1.p("p from l1")
    l1.print("print from l1")
    l1.log("log from l1")
    l1.warn("warning from l1")
    l1.debug("debug from l1")

    l2.p("p from l2")
    l2.print("print from l2")
    l2.log("log from l2")
    l2.warn("warning from l2")
    l2.debug("debug from l2")

    l3.p("p from l3")
    l3.print("print from l3")
    l3.log("log from l3")
    l3.warn("warning from l3")
    l3.debug("debug from l3")

    l4.p("p from l4")
    l4.print("print from l4")
    l4.log("log from l4")
    l4.warn("warning from l4")
    l4.debug("debug from l4")

    #logger = Messager(logFile="test/externalWarningsAndErrors.log", logFormat=logFormats["externalWarning"])
    ##TeeStderr(logger)
    ##import testLog2
    #l4.p("p from l1")
    #l4.print("print from l4")
    #l4.log("log from l4")
    #l4.warn("warning from l4, balabla")
    #l4.debug("debug from l4")
    #raise l1.error("error from l1")

