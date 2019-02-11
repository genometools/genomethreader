#!/bin/sh -ex

clear
cd $HOME/work/genometools
make 64bit=yes cairo=no opt=no $*
cd $HOME/work/genomethreader
make 64bit=yes licensemanager=no opt=no $*
