#!/bin/bash

RSTDOC=ds10.doc.sh.rst
echo "" > $RSTDOC
title () {
    tput bold; tput setaf 6; echo "${title}"; tput setaf 7; tput sgr0
    echo "${title}::" >> $RSTDOC
    echo "" >> $RSTDOC
    read -p"press ENTER to continue"
}

ndoing () {
    echo 'skip'
}
nsdoing () {
    echo 'skip'
}

doing () {
    title
    tput setaf 2; echo "\$ ${todo}"; tput setaf 7;
    read -p"press ENTER to run this command"
    eval ${todo}
    echo "    \$ ${todo}" >> $RSTDOC
}
sdoing () {
    tput setaf 2; echo "\$ ${todo}"; tput setaf 7;
    read -p"press ENTER to run this command"
    eval ${todo} || exit 1
    echo "    \$ ${todo}" >> $RSTDOC
    read -p"press ENTER to continue"
}

title="Open an empty ds9 window. You must exit it manually to type 'exit()' in the ipython console or use Ctrl-D."
todo="ds10"
doing

title="Open a fits in existing ds9 window. '-a' means open all extension in the fits file, it is equivalent to use '-f0,1', or '-f0~1', but extension 0 have no data, so we only add 1 new frame into the ds9"
todo="ds10 N20150904S0254.fits -a --exit"
doing

title="Open fits from file list 'test_list.txt' in a new ds9 window with name 'test' using the '--name=' argument. Only *.txt files are identified as file list, other file will be opened as fits file."
todo="ds10 test_list.txt -f1 --name=test --exit"
doing

title="ds10 will open 2d image by ds9, 1d array by matplotlib and table in ipython interactive console, use d.t to see all tables and d.tn to see names of all tables, in ipython"
todo="ds10 N20150903S0254.fits bbody* spec-* -a"
doing

title="Filter files and extensions to open using -b and -B arguments, here we also use --notOpen and --exit argument to not really open this files(only show their filter results) and exit the interactive mode"
todo="imheader N2015090* bbody*  -k\"object,date,exptime,FRMNAME\" -s -a"
doing
todo="ds10 N2015090* bbody* -a --notOpen --exit"
sdoing
todo="ds10 N2015090* bbody* -b\"{exptime}\" -a --notOpen --exit"
sdoing
todo="ds10 N2015090* bbody* -b\"(not){exptime}\" -a --notOpen --exit"
sdoing
todo="ds10 N2015090* bbody* -b\"(not){exptime}\" --boolFrames=0 -a --notOpen --exit"
sdoing
todo="ds10 N2015090* bbody* -b\"(not)(float({exptime})<1)(and){DATE}>'2015-09-03'\" -a --notOpen --exit"
sdoing
todo="ds10 N2015090* bbody* -B\"(not)(float({exptime})<1)(and){DATE}>'2015-09-03'\" -a --notOpen --exit"
sdoing

title="Adding new fits image into the 'test' ds9 window. We stay at the interactive mode and give you a list of command to playwith. the '--exec' options allow you to run some ds10 command after all the images are loaded"
todo="ds10 N2015*.fits -a 2Trncsci-N20160126S0150-153.fits -f0~8 --name=test_command --exec='d.h;d.hb;d.hpro'"
echo $commandList
doing

title="Here are the most interest functions of ds10: image debug tools using region, we open a new window to test it"
todo="ds10 N20150903*.fits -a 2Trncsci-N20160126S0150-153.fits -f0~4 m2fs_flat* --name=test_regions --exec='d.h;d.hr;d.hc;d.hf'"
doing
