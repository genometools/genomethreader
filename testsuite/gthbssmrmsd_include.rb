Name "gthbssmrmsd (arabidopsis)"
Keywords "gthbssmrmsd"
Test do
  run_test "#{$bin}gthbssmrmsd arabidopsis.bssm arabidopsis.bssm"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gthbssmrmsd/arabidopsis.out"
end
