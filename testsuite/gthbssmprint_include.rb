Name "gthbssmprint (arabidopsis)"
Keywords "gthbssmprint"
Test do
  run_test "#{$bin}gthbssmprint arabidopsis.bssm"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthbssmprint/arabidopsis.out"
end

Name "gthbssmprint (human)"
Keywords "gthbssmprint"
Test do
  run_test "#{$bin}gthbssmprint human.bssm"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthbssmprint/human.out"
end
