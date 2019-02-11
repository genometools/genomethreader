#!/bin/sh -e
#
# Copyright (c) 2005-2006 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

range=
force=
verbose=
gzip=
bzip2=
range_value=5;

usage() {
  echo "Usage: $0 [-r range] [-f] [-v] [-g | -b] intermediate_files" >&2
  exit 1
}

while getopts fvgbr: name
do
    case $name in
    r)    range=1
          range_value="$OPTARG";;
    f)    force=1;;
    v)    verbose=1;;
    g)    gzip=1;;
    b)    bzip2=1;;
    ?)    usage;;
    esac
done

range_opt=
force_opt=
verbose_opt=
gzip_opt=
bzip2_opt=

if [ -n "$force" ]
then
  force_opt="-force"
fi

if [ -n "$range" ]
then
  range_opt="-range $range_value"
fi

if [ -n "$verbose" ]
then
  verbose_opt="-v"
fi

if [ -n "$gzip" ]
then
  gzip_opt="-gzip"
fi

if [ -n "$bzip2" ]
then
  bzip2_opt="-bzip2"
fi

shift $(($OPTIND - 1))

if [ $# -eq 0 ]
then
  usage
fi

if [ -n "$verbose" ]
then
  echo '$ split according to alignment score'
fi

gthsplit -alignmentscore $range_opt $force_opt $verbose_opt $gzip_opt $bzip2_opt $*

if [ -n "$verbose" ]
then
  echo '$ split according to coverage'
fi

for filenamebase in $*
do
  upperrange=100
  lowerrange=`expr 100 - $range_value`
  while [ $lowerrange -ge 0 ]
  do
    filename="$filenamebase.scr$lowerrange-$upperrange"
    if [ -r "$filename" ]
    then
      gthsplit -coverage $range_opt $force_opt $verbose_opt $gzip_opt $bzip2_opt $filename
      rm -f $filename
      if [ -n "$verbose" ]
      then
        echo "\$ split file removed: $filename"
      fi
    fi
    upperrange=$lowerrange
    lowerrange=`expr $lowerrange - $range_value` || true
  done
  shift
done
