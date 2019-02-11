Name "gthconsensus -help"
Keywords "gthconsensus"
Test do
  run_test "#{$bin}gthconsensus -help"
  run "grep -v ^Usage: #{$last_stdout}"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthconsensus/gthconsensus_-help.out"
end

Name "gthconsensus -help+"
Keywords "gthconsensus"
Test do
  run_test "#{$bin}gthconsensus -help+"
  run "grep -v ^Usage: #{$last_stdout}"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthconsensus/gthconsensus_-help+.out"
end

Name "gthconsensus empty_file"
Keywords "gthconsensus"
Test do
  run_test("#{$bin}gthconsensus #{$testdata}empty_file", :retval => 1)
  grep $last_stderr, "no element found"
end

Name "gthconsensus empty_file.gz"
Keywords "gthconsensus"
Test do
  run_test("#{$bin}gthconsensus #{$testdata}empty_file.gz", :retval => 1)
  grep $last_stderr, "no element found"
end

Name "gthconsensus empty_file.bz2"
Keywords "gthconsensus"
Test do
  run_test("#{$bin}gthconsensus #{$testdata}empty_file.bz2", :retval => 1)
  grep $last_stderr, "no element found"
end
