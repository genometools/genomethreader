#!/usr/bin/env ruby
#
# Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg
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

# the GenomeThreader test suite (employs ``stest'').

if $0 == __FILE__
  $:<< "."            # favor the local stest version
  require 'stest'
  at_exit do
    OnError do exit 1 end
  end
end

# set some global variables
if $arguments["path"] then
  $path=File.join($arguments["path"], "")
else
  $path=""
end

if $arguments["testdata"] then
  $testdata=File.join($arguments["testdata"], "")
else
  $testdata=File.join(Dir.pwd, "..", "testdata", "")
end

if $arguments["bin"] then
  $bin=File.join($arguments["bin"], "")
else
  $bin=File.join(Dir.pwd, "..", "bin", "")
end

if $arguments["cur"] then
  $cur=$arguments["cur"]
else
  $cur=File.join(Dir.pwd, "..", "")
end

$transdir=File.join(Dir.pwd, "..", "gtdata" , "trans", "")
$obodir=File.join(Dir.pwd, "..", "gtdata" , "obo_files", "")

$scriptsdir=File.join(Dir.pwd, "..", "scripts", "")

if $arguments["gthtestdata"] then
  $gthtestdata=File.join($arguments["gthtestdata"], "")
end

$systemname=`uname -s`
$systemname.chomp!

# define helper functions
def run_test(str, opts = {})
  if $arguments["memcheck"] then
    if $systemname == "Linux" then
      $memcheck = "valgrind --tool=memcheck --suppressions="+
                  File.join($testdata, "gth.supp")+
                  " --leak-check=yes --error-exitcode=1 -q"
    elsif $systemname == "OpenBSD" then
      $memcheck = "env MALLOC_OPTIONS='GJ'"
    end
  else
    $memcheck = ""
  end
  run("#{$memcheck} #{$path}#{str}", opts)
end

def file_exists(filename, mtimes)
  if not File.exist?(filename) then
    raise TestFailed, "file #{filename} does not exist"
  end
  if mtimes then
    mtimes[filename] = File.mtime(filename)
  end
end

def cmp_mtimes(filename, mtimes)
  if (File.mtime(filename) != mtimes[filename]) then
    raise TestFailed, "mtimes differ for file #{filename}"
  end
end

# include the actual test modules
require 'gth_include'
require 'gthconsensus_include'
require 'gthbssmfileinfo_include'
require 'gthbssmprint_include'
require 'gthbssmrmsd_include'
require 'gthbssmtrain_include'
require 'gthunit_include'
require 'fastdp_include'
require 'ngasp_include'
require 'ngasp_large_include'
