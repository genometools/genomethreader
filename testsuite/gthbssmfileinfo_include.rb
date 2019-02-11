Name "gthbssmfileinfo (arabidopsis)"
Keywords "gthbssmfileinfo"
Test do
  run_test "#{$bin}gthbssmfileinfo arabidopsis"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthbssmfileinfo/arabidopsis.out"
end

Name "gthbssmfileinfo (human)"
Keywords "gthbssmfileinfo"
Test do
  run_test "#{$bin}gthbssmfileinfo human"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthbssmfileinfo/human.out"
end
