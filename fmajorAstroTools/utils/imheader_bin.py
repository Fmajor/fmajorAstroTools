from __future__ import division, absolute_import, print_function
from fmajorAstroTools.utils import imheader
import getopt
import sys

usageStr="""imheader filename [options]
    filename: fits files or a file list(with any suffix except for '.fits')

    without options:
        print the full header of the 0 extension
    options:
        -d: print the header using dict form for the first file and first extension
        -e: only print the extension info

        -k key1,key2...: only print key1,key2 in the header
        -s: use short print
            if use short print, must also set the keys using -k
        -a: use all frame
        -f num: use frame num in each file, Examples -f1,2,3 -f1~3,5~9 -f1:10:2
        -b boolExp: filtering the files using bool expression
           examples:    ({EXPTIME}>10(and){OBJTYPE}='object')(or)({RA}>100(and){DEC}<10)
                        ('acq'(notin){OBJECT})
               all the keys must be inside "{}"
        -o: only print for the first extension of the first file
        -v: have debug output
        --sn, --shortName: remove dirname
        --onlyKeys: only list the keywords in some header alphybetically
        --onlyFile: only print filename
        --onlyFileExt: print filename.fits[ext]
        --boolFrames: specify which frame to be use for bool filter (default is all frames)
        --diff: output and give command to vimdiff two pairs of headers
        --output: output the header to *.fits.header
        --withoutExt: not show ext num in the short print mode
        --fileNameExp: do operation on the output filename in short print mode, example:
            --fileNameExp='os.path.abspath({})'
            --fileNameExp='{}[:-2]'
        --show: equal -a --onlyFileExt
        --debug: print debug info
"""
def usage():
    print(usageStr)
    sys.exit(0)

c = {
    "fileNameType": "relative",
    "headerMode": "long",
    "mode": "",
    "only": False,
    "keys": "",
    "frames": "0",
    "allFrames": False,
    "boolFrames": "",
    "boolExp": "",
    "extBoolExp": "",
    "listColumn": 0,
    "listDelimite": "",
    "fitsExts": "",
    "output": False,
    "onlyList": False,
    "debug": False,
    "withoutExt": False,
}

def main(argv):
    opts, args = getopt.gnu_getopt(argv[1:],
            "aB:b:def:k:sov", ["diff", "boolFrames=", "output", "sn",
                           "onlyKeys", "onlyFile", "onlyFileExt", "show",
                           "exp=", "withoutExt", "debug", "fileNameExp="])

    debug = False
    if len(args)==0:
        usage()
    for op, value in opts:
        if op == "-d":
            c["headerMode"] = "dict"
        elif op == "-e":
            c["mode"] = "e"
        elif op == "-k":
            c["keys"] = value
        elif op == "-s":
            c["headerMode"] = "short"
        elif op == "-a":
            c["allFrames"] = True
        elif op == "-f":
            c["frames"] = value
        elif op == "-b":
            c["boolExp"] = value
        elif op == "-B":
            c["extBoolExp"] = value
        elif op == "-o":
            c["only"] = True
        elif op == "-v":
            c["debug"] = True
        elif op == "--exp":
            c["boolExp"] = value
        elif op == "--sn":
            c["fileNameType"] = "short"
        elif op == "--onlyKeys":
            c["headerMode"] = "keys"
        elif op == "--onlyFile":
            c["mode"] = "f"
        elif op == "--onlyFileExt":
            c["mode"] = "fe"
        elif op == "--show":
            c["mode"] = "fe"
            c["allFrames"] = True
        elif op == "--boolFrames":
            c["boolFrames"] = value
        elif op == "--diff":
            c["headerMode"] = "diff"
        elif op == "--output":
            c["output"] = True
        elif op == "--debug":
            debug = True
        elif op == "--fileNameExp":
            c["fileNameExp"] = value
        elif op == "--withoutExt":
            c["withoutExt"] = True
        else:
            print("unknow args: {}:{}".format(op, value))
            usage()

    if debug:
        print("args: {}".format(args))
        print("configs:{}".format(c))
    result = imheader.main(args, configs=c)
    return result

if __name__=="__main__":
    main(sys.argv)
