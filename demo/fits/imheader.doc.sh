#!/bin/bash

title () {
    tput bold; tput setaf 6; echo "${title}"; tput setaf 7; tput sgr0
}

doing () {
    echo 'skip'
}
sdoing () {
    echo 'skip'
}

doing () {
    title
    tput setaf 2; echo "\$ ${todo}"; tput setaf 7;
    eval ${todo}
}
sdoing () {
    tput setaf 2; echo "\$ ${todo}"; tput setaf 7;
    eval ${todo} || exit 1
}

title="with no options: output the full header of extension 0 for this file"
todo="imheader N20150903S0254.fits"
doing

title="work for file list:output the full header of extension 0 for 𝔼 file in a file list"
todo="cat test_list.txt"
doing
todo="imheader test_list.txt"
sdoing

title="Overview fits struct using '--show'"
todo="imheader N20150903* vstandard_comb.fits --show"
doing

title="select key names to print using '-k'"
todo="imheader N20150903*.fits -k\"object,date,exptime"""
doing

title="print keys in 'short' mode (the prints above are called 'long' mode) using '-s', you use use '-k' at the same time to select keys for printing. If some key not exists in some extensions, they are markd as '-'"
todo="imheader *.fits -k\"object,date,exptime\" -s"
doing

title="print more extensions also ⊇ '-a' or '-f'"
todo="imheader N20150903S0254.fits vstandard_comb.fits -k\"object,date,exptime\" -s -a"
doing

title="or you can add the extension name directly to a filename like this (you can not do this in zsh, since zsh parse '[]' itself)"
todo="imheader N20150903S0254.fits[0] -k\"object,date,exptime\" -s"
doing
todo="cat test_list_ext.txt"
sdoing
todo="imheader test_list_ext.txt -k\"object,date,exptime\" -s"
sdoing

title="print the header as colorful dict using '-d' (you can see the color in the terminal instead in the doc here, the first key of ∀ character are colored red)"
todo="imheader vstandard_comb.fits -d"
doing

title="only print key names using '--onlyKeys'"
todo="imheader vstandard_comb.fits --onlyKeys"
doing

title="filter filer using bool operation, '-b', there should be no blank between '(and)', '(or)', '(in)' and '(notin)'"
todo="imheader N2015090* bbody*  -k\"object,date,exptime,FRMNAME\" -s -a"
doing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"{exptime}\""
sdoing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"(not){exptime}\""
sdoing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"float({exptime})>1\""
sdoing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"(not)(float({exptime})<1)(and){DATE}>'2015-09-03'\""
sdoing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"'HIP'(in){OBJECT}\""
sdoing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"'078'(in){FRMNAME}\""
sdoing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"'078'(in){FRMNAME}\" --boolFrames=1"
sdoing

title="you can output only the filename or filename with extension number using '--onlyFile' and '--onlyFileExt'"
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"'HIP'(in){OBJECT}\" --onlyFile"
doing
todo="imheader N2015090* bbody*  -k\"object,date,exptime\" -s -b\"'HIP'(in){OBJECT}\" -a --onlyFileExt"
sdoing

