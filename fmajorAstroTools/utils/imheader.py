from __future__ import division, absolute_import, print_function
try:
    input = raw_input
except NameError:
    pass

import colorama
from colorama import Fore, Back, Style
#colorama.init(autoreset=True)

import numpy as np
import sys
from astropy.io import fits as pyfits
#import pyfits
import sys, getopt
import pprint
import re
import os
import warnings
import pdb
warnings.filterwarnings("ignore")

#<== fucntions
def isFitsFile(name, fitsExts):
    for eachExt in fitsExts:
        N = len(eachExt)
        # for *.fits[frame] version
        if name.endswith(eachExt) or "{}[".format(eachExt) in name:
            return True
    return False

def frameFileSep(name, fitsExts):
    for eachExt in fitsExts:
        if "{}[".format(eachExt) in name:
            return eachExt
    return ""

def alphaBetaPrint(l):
    ''' print list of word alphabetaly, case sensitive '''
    toPrint = sorted(l)
    last = ""
    for each in toPrint:
        if each[0]!=last:
            last = each[0]
            if last:
                print()
        print(each,end=" ")
    if len(toPrint):
        print()

def safe_mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def testListFile(name, maxListSize):
    with open(name, "rb") as f:
        text = f.read()
        if len(text) > maxListSize:
            print("file {} is too large as a list, do you want to continue use it as a list??\n\t Press Enter to continue, Ctrl+c to excape")
            aux = input()

def getFilePrintName(name,fileNameType):
    # [short|relative|full]
    if fileNameType=="relative":
        return name
    elif fileNameType=="full":
        return os.path.abspath(name)
    elif fileNameType=="short":
        return os.path.basename(name)
    else:
        raise Exception("unknown fileNameType:{}".format(fileNameType))

def getDimStr(header):
    extName = header.get("extname", "")
    dim = header.get("NAXIS", 0)
    dims = ", ".join([str(header.get("NAXIS"+str(each),"?")) for each in range(1, dim+1)])
    obj = header.get("object", "")
    if not dim:
        return "{} {} {}".format(extName, dim, obj)
    else:
        return "{} {} [{}] {}".format(extName, dim, dims, obj)

inRe    = re.compile(r"(\(in\))")
notInRe = re.compile(r"(\(notin\))")
andRe   = re.compile(r"(\(and\))")
orRe    = re.compile(r"(\(or\))")
notRe   = re.compile(r"(\(not\))")
#==> functions

def main(args, configs={}): #<==
    """
    Parameters
    ------------
    args : list
        all input file name (could be ascii list file or fits file)

    Parameters in configs
    -------------
    fileNameType : {"relative", ["short"|"relative"|"full"|"only"] }, optional
        short: only filename
        relative: relative name
        full:full path
    fileNameExp : {"", str}
        do expression for fileName in short print mode, example:
            "{}[:-3]"
            "os.path.abspath({})"
    only : {False, bool}
        only print things for the first ext of the first file
    mode : {"", [""|"f"|"e"|"fe"|"feq"]}, optional
        f : only print file name, return file name
        e : print file name plus extension info, return file name
        fe : print file name plus extension info, return filename[ext]
        feq : no print, return filename[ext]
        not set: print the headers
    keys : {"", str separated by ","}
        keys to showup
    headerMode : {"long", ["long"|"short"|"dict"|"keys"|"diff"]}, optional
        how to print header
    frames : {"0", str separated by ","}, optional
        default frames to print
    allFrames : bool, optional
        when allFrames is on, default frames is all
    boolExp : {"", str}, optional
        use bool expression to filter fits file
            no blanks in the expression!!
        examples:
            -b"np.abs({exptime})>100"
            -b"({obstype}=='object')(and)(({object}=='type')(or)({object}=='flat'))"
            -b"{obstype}(in)['object','flat']"
            -b"{obstype}(notin)['object','flat']"
    extBoolExp: boolExp for each ext
    boolFrames : {"", str seperated by ","}, optional
        which frame to be used if boolExp is set
    listColumn : {0, int}, optional
        which column in the list file to get the fits file
    listDelimiter : {"", str}, optional
        what delimiter to use in the list file to get the column number
    fitsExts : {"fits", str}, optional
        define costum fits file extension(e.g. fits.gz)
    output: {False, bool}, optional
        print to *.fits.headers instead of to the screen
    onlyList: {False, bool}, optional
        only print the key list
    withoutExt: {False, bool}, optional
        do not show ext num in short print
    """#==>

    #<== setups
    fileNameType = configs.get("fileNameType", "relative")
    fileNameExp = configs.get("fileNameExp", "")
    headerMode = configs.get("headerMode", "long")
    mode = configs.get("mode", "")
    only = configs.get("only", False)
    keys = configs.get("keys", "")
    frames = configs.get("frames", "0")
    allFrames = configs.get("allFrames", False)
    boolFrames = configs.get("boolFrames", "")
    boolExp = configs.get("boolExp", "")
    extBoolExp = configs.get("extBoolExp", "")
    listColumn = configs.get("listColumn", 0)
    listDelimiter = configs.get("listDelimiter", "")
    fitsExts = configs.get("fitsExts", "")
    output = configs.get("output", False)
    onlyList = configs.get("onlyList", False)
    debug = configs.get("debug", False)
    withoutExt = configs.get("withoutExt", False)
    #==> setups

    fileSeperator = "="*50
    separator = "-"*50
    maxListSize=1*1024*1024
    # after config process
    if not fitsExts:
        fitsExts="fits"
    else:
        fitsExts+=",fits"
    fitsExts = [each for each in fitsExts.split(",") if each]

    if boolFrames:
        boolFramesAux = boolFrames.split(",")
        boolFrames = [int(each) for each in boolFrames]

    if keys:
        keysaux = keys.split(",")
        keys = [each.strip().upper() for each in keysaux]
    aux = frames.split(",")
    auxFrames = []
    for each in aux:
        if '~' in each:
            _ = each.split('~')
            if len(_) != 2:
                raise ValueError('error format in "-f" argument: {}'.format(frames))
            auxFrames.extend(list(range(int(_[0]), int(_[1])+1)))
        elif ':' in each:
            _ = each.split(':')
            try:
                if len(_)==2:
                    __ = list(range(int(_[0]), int(_[1])))
                elif len(_)==3:
                    __ = list(range(int(_[0]), int(_[1]), int(_[2])))
                else:
                    raise ValueError()
            except ValueError:
                raise ValueError("if you use ':' in -f, it should be like 1:10 or 1:10:2, slices like 1::2 is not allown")
            auxFrames.extend(__)
        else:
            auxFrames.append(int(each))
    checkFrame = False
    if frames!="0":
        checkFrame = True
    frames = [int(each) for each in auxFrames]

    # convert bool exp string to python code
    if boolExp:
        boolExp = inRe.sub(" in ", boolExp)
        boolExp = notInRe.sub(" not in ", boolExp)
        boolExp = andRe.sub(" and ", boolExp)
        boolExp = orRe.sub(" or ", boolExp)
        boolExp = notRe.sub(" not ", boolExp)
    if extBoolExp:
        extBoolExp = inRe.sub(" in ", extBoolExp)
        extBoolExp = notInRe.sub(" not in ", extBoolExp)
        extBoolExp = andRe.sub(" and ", extBoolExp)
        extBoolExp = orRe.sub(" or ", extBoolExp)
        extBoolExp = notRe.sub(" not ", extBoolExp)

    if onlyList:
        only=True
        headerMode="keys"

    #<== setup todos
    #!! set todo files
    # make all file to be a list
    # each have format *.fits or *.fits[*]
    todoFiles = []
    skipExts = ["headers"]
    #auxArgs = [each[2:] if each[:2]=="./" else each for each in args]
    #for eachFile in auxArgs:
    for eachFile in args:
        if eachFile.split(".")[-1] in skipExts:
            continue
        if isFitsFile(eachFile, fitsExts):
            todoFiles.append(eachFile)
        else:
            testListFile(eachFile, maxListSize)
            f = open(eachFile)
            if listDelimiter:
                thisTodoFiles = [each.split(listDelimiter)[listColumn] for each in f]
            else:
                thisTodoFiles = [each.split()[listColumn] for each in f]
            f.close()
            todoFiles.extend(thisTodoFiles)

    #!! set todo exts
    # now to get extensions for each file
    reExpFrame = re.compile(r"\[(.*)\]")
    todoFrames = ["?"] * len(todoFiles)
    for i,eachFile in enumerate(todoFiles):
        thisExt = frameFileSep(eachFile, fitsExts)
        if not thisExt:# not have [*]
            todoFrames[i] = [frames, 'all'][allFrames]
        else:
            checkFrame = True
            aux = eachFile.split("."+thisExt)
            assert len(aux)==2
            todoFiles[i] = aux[0]+"."+thisExt
            thisFrame = reExpFrame.findall(aux[1])
            assert len(thisFrame)==1
            todoFrames[i] = [int(each.strip()) for each in thisFrame[0].split(",")]

    #!! try to get file obj length eiter in bool test or do it seperately
    fileSize = []
    #!! use boolExp to filter file if it is set
    reExp = re.compile(r"\{(.*?)\}")
    keyFilterList=reExp.findall(boolExp)
    filesToUse = []
    if boolExp:
        print("filter: {}, boolFrames: {}".format(boolExp, boolFrames))
        lastFile=""
        for eachFile, eachFrames in zip(todoFiles, todoFrames): # filter
            if lastFile!=eachFile:
                thisObjs=pyfits.open(eachFile)
                N = len(thisObjs)
            fileSize.append(N)
            if eachFrames=="all":# try all frames for the bool operation
                eachFrames=range(len(thisObjs))
            if len(boolFrames):# only use the boolFrames to do bool operation
                eachFrames=boolFrames
            for eachFrame in eachFrames:
                if eachFrame>=N:
                    continue
                eachHeader = thisObjs[eachFrame].header
                thisKeyParis = [(eachKey, eachHeader.get(eachKey,"")) if not isinstance(eachHeader.get(eachKey,""), str)
                                else (eachKey, "'{}'".format(eachHeader.get(eachKey,"")))
                                    for eachKey in keyFilterList]
                todoDict = dict(thisKeyParis)
                toExec = boolExp.format(**todoDict)
                l=locals()
                try:
                    exec("aux="+toExec, globals(), l)
                except ValueError:
                    pass
                except Exception as e:
                    print(e)
                    raise Exception('error when exec boolExp\naux={}'.format(toExec)) from None
                aux=l['aux']
                if aux:
                    if eachFile not in filesToUse:
                        filesToUse.append(eachFile)
                        if only:# if only, only use the first file
                            break

        # update todoFiles
        aux1=[]; aux2=[]; aux3=[]
        for eachFile, eachFrame, eachSize in zip(todoFiles, todoFrames, fileSize):
            if eachFile in filesToUse:
                aux1.append(eachFile)
                aux2.append(eachFrame)
                aux3.append(eachSize)
        todoFiles = aux1; todoFrames = aux2; fileSize = aux3
    elif checkFrame:
        for eachFile in todoFiles:
            fileSize.append(len(pyfits.open(eachFile)))

    #!! drop files and frames by obj length
    if checkFrame:
        aux1=[]; aux2=[]
        for eachFile, eachFrames, eachSize in zip(todoFiles, todoFrames, fileSize):
            aux = [each for each in eachFrames if each<eachSize]
            auxDrop = [str(each) for each in eachFrames if each>=eachSize]
            if len(auxDrop):
                print("does not exists: {}[{}]".format(eachFile, ", ".join(auxDrop)), file=sys.stderr)
            if len(aux):
                aux1.append(eachFile)
                aux2.append(aux)
            else:
                print("file: {} all dropped".format(eachFile), file=sys.stderr)

        todoFiles = aux1; todoFrames = aux2

    if len(todoFiles)<1:
        print("no fits to do!")
        return {"names":[], "todoFiles":[], "todoFrames":[]}
    #!! if only
    if only:
        todoFiles  =  todoFiles[:1]
        todoFrames = todoFrames[:1]
        if not allFrames:
            todoFrames[0] = todoFrames[0][:1]

    #==> setup todos

    #print(fitsExts, boolExp, boolFrames, keys, frames)
    if debug:
        print("todoFiles:", todoFiles, "  todoFrames:", todoFrames)
    #<== begin to run
    # [f|e|fe|feq|""]
    expandFrames = []
    if extBoolExp:
        keyFilterList=reExp.findall(extBoolExp)
    if mode: # simple print form
        fileNameReturns = []
        printed = False
        for i, (eachFile, eachFrames) in enumerate(zip(todoFiles, todoFrames)):
            try:
                obj = pyfits.open(eachFile)
            except Exception as e:
                print('open {} error'.format(eachFile))
                raise
            objN = len(obj)
            filePrintName = getFilePrintName(eachFile, fileNameType)
            if eachFrames=="all":
                eachFrames=range(len(obj))
            if extBoolExp:
                _eachFrames = []
                if not printed:
                    print("ext filter: {}".format(extBoolExp))
                    printed = True
                lastFile=""
                thisObjs=pyfits.open(eachFile)
                N = len(thisObjs)
                fileSize.append(N)
                for eachFrame in eachFrames:
                    if eachFrame>N:
                        continue
                    eachHeader = thisObjs[eachFrame].header
                    thisKeyParis = [(eachKey, eachHeader.get(eachKey,"")) if not isinstance(eachHeader.get(eachKey,""), str)
                                    else (eachKey, "'{}'".format(eachHeader.get(eachKey,"")))
                                        for eachKey in keyFilterList]
                    todoDict = dict(thisKeyParis)
                    toExec = extBoolExp.format(**todoDict)
                    l=locals()
                    try:
                        exec("aux="+toExec, globals(), l)
                    except ValueError:
                        continue
                    except Exception:
                        raise
                    aux=l['aux']
                    #print("toExec:{} result:{}".format(toExec, aux))
                    if aux:
                        _eachFrames.append(eachFrame)
                eachFrames = _eachFrames
            expandFrames.append(eachFrames)

            fnWithExt = ["{}[{}]".format(filePrintName, eachExt)
                                  for eachExt in eachFrames]
            fnAll = [" ".join([
                        eachPrint,
                        obj[eachExt].__module__.split("astropy.io.fits.")[-1],
                        getDimStr(obj[eachExt].header)])
                            for eachPrint, eachExt in zip(fnWithExt, eachFrames)]
            if "feq" in mode:# no print, return filename[ext]
                fileNameReturns += fnAll
            elif "fe" in mode:# print file name plus extension info, return filename[ext]
                for eachPrint, eachExt in zip(fnWithExt, eachFrames):
                    eachObj = obj[eachExt]
                    dimsStr = getDimStr(eachObj.header)
                    eachObjStr = eachObj.__module__.split("astropy.io.fits.")[-1]
                    print(eachPrint, eachObjStr, dimsStr)
                fileNameReturns += fnWithExt
            elif "f" in mode:# only print file name, return file name
                print(filePrintName)
                fileNameReturns.append(filePrintName)
            elif "e" in mode:# print file name plus extension info, return file name
                print(fileSeperator)
                print(filePrintName)
                fileNameReturns.append(filePrintName)
                print(separator)
                for eachPrint, eachExt in zip(fnWithExt, eachFrames):
                    if eachExt<objN:
                        eachObj = obj[eachExt]
                        dimsStr = getDimStr(eachObj.header)
                        eachObjStr = eachObj.__module__.split("astropy.io.fits.")[-1]
                        print(eachPrint, eachObjStr, dimsStr)
            else:
                raise Exception("unknown mode:{}".format(mode))
        result = {"names":fileNameReturns, "todoFiles":todoFiles, "todoFrames":expandFrames}
        return result

    diffDir = "imheader_diff"
    if headerMode=="diff": # only make command file here
        safe_mkdir(diffDir)
        if len(todoFiles)%2>0: # the first half is A file, the last half is B file
            raise Exception("todoFiles:{} should be even number".format(todoFiles))
        nTotal = len(todoFiles)
        todoFramesStr = [[[str(each) for each in thisFrameStr], ["all"]][thisFrameStr=="all"] for thisFrameStr in todoFrames]
        allVimFiles = ["{}-{}-{}_{}_.txt".format(i%int(nTotal/2), nTotal, eachFits, ",".join(todoFramesStr[i]))
                       for i, eachFits in enumerate(todoFiles)]
        allCommands = ["vimdiff {} {}\n".format(allVimFiles[i], allVimFiles[i+int(nTotal/2)]) for i in range(int(nTotal/2))]
        commandF = open(os.path.join(diffDir, "command"), "w")
        commandF.writelines(allCommands)
        commandF.close()

    #!! prepare for the short print array
    if headerMode=="short":
        if len(keys)==0:
            raise Exception("you must give some keys to make a short print!")
        dtype=[("FileName", "U128")] + [(str(each), "U64") for each in keys]
        arrayList = np.empty(0, dtype=dtype)

    lastFits = ""
    stdout = sys.stdout
    #!! main loop
    for i, eachFits in enumerate(todoFiles):
        if eachFits!=lastFits:
            thisObjs = pyfits.open(eachFits)
            lastFits = eachFits
        thisFrames = todoFrames[i]
        thisPrintName = getFilePrintName(eachFits, fileNameType)
        if thisFrames=="all":
            thisFrames = range(len(thisObjs))

        if headerMode in ["long", "diff", "dict", "keys"]:
            todoObjs    = [thisObjs[eachIndex] for eachIndex in thisFrames]
            todoNames   = ["{}[{}]".format(thisPrintName, eachIndex) for eachIndex in thisFrames]
            todoHeaders = [thisObjs[eachIndex].header for eachIndex in thisFrames]
            todoExtNames= [each.get("extname", "") for each in todoHeaders]
            maxFileLen = np.max([len(each) for each in todoNames])
            maxExtNameLen = np.max([len(each) for each in todoExtNames])
            #!! header info for long print (name, extName, type(obj)) 
            toPrintStr = "{{:>{}}} {{:>{}}} {{}} {{}}".format(maxFileLen, maxExtNameLen)
            if headerMode=="diff":
                todoFramesStr = [[[str(each) for each in thisFrameStr],
                                  ["all"]][thisFrameStr=="all"]
                                      for thisFrameStr in todoFrames]
                thisFile = os.path.join(diffDir, "{}-{}-{}_{}_.txt".\
                        format(i%int(nTotal/2), nTotal, eachFits, ",".join(todoFramesStr[i])))
                thisF = open(thisFile, "w")
                sys.stdout = thisF
            elif output:# output to *.fits.header, only work in long mode
                sys.stdout = stdout
                basename = os.path.basename(eachFits)
                dirname = os.path.dirname(eachFile)
                thisFile = os.path.join(dirname, "{}.headers".format(basename))
                print("output to {}".format(thisFile))
                thisF = open(thisFile, "w")
                sys.stdout = thisF

            #!! print header before long print (to screen or to file or to diff file)
            print(fileSeperator)
            for obj,name,extName in zip(todoObjs, todoNames, todoExtNames):
                dimStr = getDimStr(obj.header)
                eachObjStr = obj.__module__.split("astropy.io.fits.")[-1]
                print(toPrintStr.format(name, extName, eachObjStr, dimStr))
            #!! main work: print header content
            for thisFrame in thisFrames:
                eachHeader = thisObjs[thisFrame].header
                if len(keys)>0:
                    auxHeader = pyfits.Header()
                    allThisKeys = list(eachHeader) #!! can this work in python3?
                    for eachKey in keys:
                        if eachKey in allThisKeys:
                            index = eachHeader.index(eachKey)
                            auxHeader.append(eachHeader.cards[index])
                    eachHeader=auxHeader
                else:#!! print all keys
                    try:
                        while not eachHeader[-1]:#!! remove blank comment lines at last
                            eachHeader.pop()
                    except Exception as e:
                        pass

                print(separator)
                print("{}[{}]".format(eachFits, thisFrame))
                if len(list(eachHeader)):
                    if headerMode in ["long", "diff"]:
                        pprint.pprint(eachHeader)
                    elif headerMode=="dict": # make dict print
                        thisKeys = sorted(list(eachHeader))
                        thisKeys = list(filter(lambda _:_, thisKeys))
                        lastKey=""
                        for eachKey in thisKeys:
                            if eachKey[0]!=lastKey:
                                coloredKey = Fore.RED + Style.BRIGHT + eachKey + Style.RESET_ALL
                                lastKey = eachKey[0]
                            else:
                                coloredKey = eachKey
                            try:
                                coloredContent = Fore.GREEN + Style.BRIGHT + str(eachHeader[eachKey]) + Style.RESET_ALL
                            except:
                                coloredContent = Fore.GREEN + Style.BRIGHT + "????" + Style.RESET_ALL
                            print("{}:{}, ".format(coloredKey,coloredContent),end="")
                        print()
                    elif headerMode=="keys":
                        thisKeys = list(eachHeader)
                        alphaBetaPrint(thisKeys)
            if headerMode == "diff" or output:
                thisF.close()
            sys.stdout = stdout
        else: # short or dict print
            for thisFrame in thisFrames:
                eachHeader = thisObjs[thisFrame].header
                if withoutExt:
                    if fileNameExp:
                        todoStr = "thisPrintName = {}".format(fileNameExp.format("thisPrintName"))
                        try:
                            exec(todoStr)
                        except Exception as e:
                            raise Exception("Error when exec({})".format(todoStr))
                    toAppend = tuple([thisPrintName] +
                                 [[str(eachHeader.get(eachKey, "-")),
                                     "undefined"][isinstance(eachHeader.get(eachKey), pyfits.card.Undefined)]
                                     for eachKey in keys]) # to be test
                else:
                    toAppend = tuple(["{}[{}]".format(thisPrintName, thisFrame)] +
                                 [[str(eachHeader.get(eachKey, "-")),
                                     "undefined"][isinstance(eachHeader.get(eachKey), pyfits.card.Undefined)]
                                     for eachKey in keys]) # to be test
                #toAppend = tuple(["{}[{}]".format(thisPrintName, thisFrame)] +
                                 #[str(eachHeader.get(eachKey, "-")) for eachKey in keys])
                arrayList = np.append(arrayList, np.array(toAppend, dtype=dtype))

    if headerMode=="diff": # if not short print, return
        with open(os.path.join(diffDir, "command")) as f:
            print(f.read())
        return
    elif headerMode != "short":
        return
    else: # last work for short print
        keys = ["FileName"] + keys
        nColumn = len(keys)
        formatStr="{{:>{}}} "*(nColumn) # add file list
        allWidth = [0] * nColumn
        for i,eachKey in enumerate(keys):
            thisValueMaxWidth = np.max([len(eachValue) for eachValue in arrayList[eachKey]])
            thisMaxWidth = np.max([len(eachKey), thisValueMaxWidth])
            allWidth[i] = thisMaxWidth
        formatStr = formatStr.format(*allWidth)
        print(formatStr.format(*keys))
        for i,each in enumerate(arrayList.tolist()):
            print(formatStr.format(*each))
    #==> end of run 


