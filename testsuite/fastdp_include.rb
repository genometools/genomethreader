Name "fastdp 01"
Keywords "gth fastdp"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}U89959_genomic.fas " +
           "-cdna #{$testdata}gth/fastdp_01.fas -fastdp"
end

Name "fastdp (U89959)"
Keywords "gth fastdp"
Test do
  run_test "#{$bin}gth -species arabidopsis " +
           "-genomic #{$testdata}U89959_genomic.fas " +
           "-cdna #{$testdata}U89959_ests.fas -gff3out -intermediate " +
           "-fastdp", :maxtime => 300
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}U89959_sas_fast.gff3"
end

Name "fastdp (U89959, introncutout)"
Keywords "gth fastdp"
Test do
  run_test "#{$bin}gth -species arabidopsis " +
           "-genomic #{$testdata}U89959_genomic.fas " +
           "-cdna #{$testdata}U89959_ests.fas -gff3out -intermediate " +
           "-fastdp -introncutout", :maxtime => 300
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}U89959_sas_fast_ic.gff3"
end

if $gthtestdata then
  Name "fastdp (unused)"
  Keywords "gth fastdp diss"
  Test do
    run_test "#{$bin}gth -bssm #{$gthtestdata}ngasp/ngasp -enrichchains " +
             "-scoreminexonlen 1 -maskpolyatails -dpminintronlen 10 " +
             "-gcmaxgapwidth 22000 -gff3out -intermediate " +
             "-genomic #{$gthtestdata}ngasp32/training_regions.fa " +
             "-cdna #{$gthtestdata}ngasp32/rna_small.fa", :maxtime => 300
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/rna_small.gff3"
  end

  Name "fastdp (-fastdp)"
  Keywords "gth fastdp diss"
  Test do
    run_test "#{$bin}gth -bssm #{$gthtestdata}ngasp/ngasp -enrichchains " +
             "-scoreminexonlen 1 -maskpolyatails -dpminintronlen 10 " +
             "-gcmaxgapwidth 22000 -gff3out -intermediate " +
             "-genomic #{$gthtestdata}ngasp32/training_regions.fa " +
             "-cdna #{$gthtestdata}ngasp32/rna_small.fa -fastdp",
             :maxtime => 120
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/rna_small_fast.gff3"
  end

  Name "fastdp (-introncutout)"
  Keywords "gth fastdp diss"
  Test do
    run_test "#{$bin}gth -bssm #{$gthtestdata}ngasp/ngasp -enrichchains " +
             "-scoreminexonlen 1 -maskpolyatails -dpminintronlen 10 " +
             "-gcmaxgapwidth 22000 -gff3out -intermediate " +
             "-genomic #{$gthtestdata}ngasp32/training_regions.fa " +
             "-cdna #{$gthtestdata}ngasp32/rna_small.fa -introncutout",
             :maxtime => 120
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/rna_small_ic.gff3"
  end

  Name "fastdp (-fastdp -introncutout)"
  Keywords "gth fastdp diss"
  Test do
    run_test "#{$bin}gth -bssm #{$gthtestdata}ngasp/ngasp -enrichchains " +
             "-scoreminexonlen 1 -maskpolyatails -dpminintronlen 10 " +
             "-gcmaxgapwidth 22000 -gff3out -intermediate " +
             "-genomic #{$gthtestdata}ngasp32/training_regions.fa " +
             "-cdna #{$gthtestdata}ngasp32/rna_small.fa -fastdp -introncutout"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/rna_small_fast_ic.gff3"
  end
end
