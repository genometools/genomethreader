Name "gthbssmtrain (examples from manual)"
Keywords "gthbssmtrain"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}U89959_genomic.fas " +
           "-cdna #{$testdata}U89959_ests.fas -gff3out -skipalignmentout " +
           "-md5ids -o arab.gff3", :maxtime => 600
  run "diff --strip-trailing-cr arab.gff3 #{$testdata}arab.gff3"
  run_test "#{$bin}gthbssmtrain -seed 1238225656 " +
           "-seqfile #{$testdata}U89959_genomic.fas arab.gff3"
  run_test "#{$bin}gthbssmbuild -gtdonor -agacceptor -datapath training_data " +
           "-bssmfile arab.bssm"
  run "diff --strip-trailing-cr arab.bssm #{$testdata}arab.bssm"
  run_test "#{$bin}gth -bssm arab -genomic #{$testdata}U89959_genomic.fas " +
           "-cdna #{$testdata}U89959_ests.fas", :maxtime => 300
end

Name "gthbssmtrain (file, no introns)"
Keywords "gthbssmtrain"
Test do
  run_test "#{$bin}gthbssmtrain -seqfile #{$testdata}U89959_genomic.fas " +
           "-matchdesc -seed 913609837 #{$testdata}U89959_csas.gff3"
  run "diff --strip-trailing-cr    #{$last_stdout} #{$testdata}gthbssmtrain.out"
  run "rm -f training_data/gthbssmtrain.run"
  run "diff --strip-trailing-cr -r training_data #{$testdata}training_data"
end

Name "gthbssmtrain (file, introns)"
Keywords "gthbssmtrain"
Test do
  run_test "#{$bin}gthbssmtrain -seqfile #{$testdata}U89959_genomic.fas " +
           "-matchdesc -seed 913609837 #{$testdata}U89959_csas_introns.gff3"
  run "diff --strip-trailing-cr    #{$last_stdout} #{$testdata}gthbssmtrain.out"
  run "rm -f training_data/gthbssmtrain.run"
  run "diff --strip-trailing-cr -r training_data #{$testdata}training_data"
end

Name "gthbssmtrain (stdin, no introns)"
Keywords "gthbssmtrain"
Test do
  run_test "#{$bin}../../genometools/bin/gt gff3 " +
           "#{$testdata}U89959_csas.gff3 | " +
           "#{$bin}gthbssmtrain -seqfile #{$testdata}U89959_genomic.fas " +
           "-matchdesc -seed 913609837 -"
  run "diff --strip-trailing-cr    #{$last_stdout} #{$testdata}gthbssmtrain.out"
  run "rm -f training_data/gthbssmtrain.run"
  run "diff --strip-trailing-cr -r training_data #{$testdata}training_data"
end

Name "gthbssmtrain (stdin, introns)"
Keywords "gthbssmtrain"
Test do
  run_test "#{$bin}../../genometools/bin/gt gff3 " +
           "#{$testdata}U89959_csas_introns.gff3 | " +
           "#{$bin}gthbssmtrain -seqfile #{$testdata}U89959_genomic.fas " +
           "-matchdesc -seed 913609837 -"
  run "diff --strip-trailing-cr    #{$last_stdout} #{$testdata}gthbssmtrain.out"
  run "rm -f training_data/gthbssmtrain.run"
  run "diff --strip-trailing-cr -r training_data #{$testdata}training_data"
end

Name "gthbssmtrain (short intron)"
Keywords "gthbssmtrain"
Test do
  run_test "#{$bin}gthbssmtrain -seqfile #{$testdata}U89959_genomic.fas " +
           "-matchdesc #{$testdata}short_intron.gff3"
  grep $last_stderr, "ignoring intron of length < 2 for sequence ID " +
                     "'md5:063b1024d68e26716b7f38caf958316f:1877523'"
end

if $gthtestdata then
  Name "gthbssmtrain (arabidopsis)"
  Keywords "gthbssmtrain"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -sort " +
             "#{$gthtestdata}gthbssmtrain/22329272.gff3 " +
             "#{$gthtestdata}gthbssmtrain/22330780.gff3"
    run_test "#{$bin}../../genometools/bin/gt csa #{$last_stdout}"
    run_test "#{$bin}../../genometools/bin/gt cds -regionmapping " +
             "#{$gthtestdata}gthbssmtrain/mapping.lua -minorflen 1 " +
             "-startcodon yes -finalstopcodon yes #{$last_stdout}",
             :maxtime => 300
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}gthbssmtrain/bssmtrain_input.gff3"
    FileUtils.mkdir_p("training_data")
    run_test "#{$bin}gthbssmtrain " +
             "-gcdonor no -goodexoncount 3 " +
             "-seed 4067582214 -intermediate -regionmapping " +
             "#{$gthtestdata}gthbssmtrain/mapping.lua "+
             "#{$gthtestdata}gthbssmtrain/bssmtrain_input.gff3"
    run "rm -f training_data/gthbssmtrain.run"
    run "diff --strip-trailing-cr -r training_data #{$gthtestdata}gthbssmtrain/training_data"
  end

  Name "gthbssmtrain (22329272 input)"
  Keywords "gthbssmtrain"
  Test do
    run_test "#{$bin}gth -species arabidopsis -introncutout " +
             "-genomic #{$gthtestdata}gthbssmtrain/22329272.fas " +
             "-cdna #{$testdata}U89959_ests.fas -intermediate -gff3out",
             :maxtime => 300
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}gthbssmtrain/22329272.gff3"
  end

  Name "gthbssmtrain (22330780 input)"
  Keywords "gthbssmtrain"
  Test do
    run_test "#{$bin}gth -species arabidopsis -introncutout " +
             "-genomic #{$gthtestdata}gthbssmtrain/22330780.fas " +
             "-cdna #{$testdata}U89959_ests.fas -intermediate -gff3out",
             :maxtime => 300
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}gthbssmtrain/22330780.gff3"
  end
end
