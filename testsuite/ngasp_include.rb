if $gthtestdata then
  Name "nGASP (check source files)"
  Keywords "gth ngasp diss"
  Test do
    run "cd #{$gthtestdata}ngasp && md5sum -c ngasp.md5sum"
  end

  Name "nGASP (training_confirmed_tidy.gff)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/confirmed.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_confirmed_tidy.gff3"
  end

  Name "nGASP (training_unconfirmed_tidy.gff)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/unconfirmed.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_unconfirmed_tidy.gff3"
  end

  Name "nGASP (test_confirmed_tidy.gff)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/test_regions_confirmed.gff3.gz"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_confirmed_tidy.gff3"
  end

  Name "nGASP (test_unconfirmed_tidy.gff)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/test_regions_unconfirmed.gff3.gz"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_unconfirmed_tidy.gff3"
  end

  Name "nGASP (training_confirmed_md5.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/training_regions.fa " +
             "#{$gthtestdata}ngasp/training_confirmed_tidy.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_confirmed_md5.gff3"
  end

  Name "nGASP (training_unconfirmed_md5.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/training_regions.fa " +
             "#{$gthtestdata}ngasp/training_unconfirmed_tidy.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_unconfirmed_md5.gff3"
  end

  Name "nGASP (test_confirmed_md5.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/test_regions.fa " +
             "#{$gthtestdata}ngasp/test_confirmed_tidy.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_confirmed_md5.gff3"
  end

  Name "nGASP (test_unconfirmed_md5.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/test_regions.fa " +
             "#{$gthtestdata}ngasp/test_unconfirmed_tidy.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_unconfirmed_md5.gff3"
  end

  Name "nGASP (training_confirmed_sorted.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -sort " +
             "#{$gthtestdata}ngasp/training_confirmed_md5.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_confirmed_sorted.gff3"
  end

  Name "nGASP (training_unconfirmed_sorted.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -sort " +
             "#{$gthtestdata}ngasp/training_unconfirmed_md5.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_unconfirmed_sorted.gff3"
  end

  Name "nGASP (test_confirmed_sorted.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -sort " +
             "#{$gthtestdata}ngasp/test_confirmed_md5.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_confirmed_sorted.gff3"
  end

  Name "nGASP (test_unconfirmed_sorted.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -sort " +
             "#{$gthtestdata}ngasp/test_unconfirmed_md5.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_unconfirmed_sorted.gff3"
  end

  Name "nGASP (training_ref2.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt merge " +
             "#{$gthtestdata}ngasp/training_confirmed_sorted.gff3 " +
             "#{$gthtestdata}ngasp/training_unconfirmed_sorted.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_ref2.gff3"
  end

  Name "nGASP (test_ref2.gff3)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt merge " +
             "#{$gthtestdata}ngasp/test_confirmed_sorted.gff3 " +
             "#{$gthtestdata}ngasp/test_unconfirmed_sorted.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_ref2.gff3"
  end

  Name "nGASP (training_ref1.gff3, combined)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/confirmed.gff3 | " +
             "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/training_regions.fa | " +
             "#{$bin}../../genometools/bin/gt gff3 -sort"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_ref1.gff3"
  end

  Name "nGASP (training_ref2.gff3, combined)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/confirmed.gff3 " +
             "#{$gthtestdata}ngasp/unconfirmed.gff3 | " +
             "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/training_regions.fa | " +
             "#{$bin}../../genometools/bin/gt gff3 -sort"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_ref2.gff3"
  end

  Name "nGASP (test_ref1.gff3, combined)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/test_regions_confirmed.gff3.gz | " +
             "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/test_regions.fa | " +
             "#{$bin}../../genometools/bin/gt gff3 -sort"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_ref1.gff3"
  end

  Name "nGASP (test_ref2.gff3, combined)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt gff3 -tidy -addids no " +
             "#{$gthtestdata}ngasp/test_regions_confirmed.gff3.gz " +
             "#{$gthtestdata}ngasp/test_regions_unconfirmed.gff3.gz | " +
             "#{$bin}../../genometools/bin/gt id_to_md5 -usedesc -seqfile " +
             "#{$gthtestdata}ngasp32/test_regions.fa | " +
             "#{$bin}../../genometools/bin/gt gff3 -sort"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/test_ref2.gff3"
  end

  Name "nGASP (protein_matches.fa)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt seqtransform -addstopaminos " +
             "#{$gthtestdata}ngasp32/protein_matches_complete_filtered.fa",
             :maxtime => 300
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp32/protein_matches.fa"
  end

  Name "nGASP (splice site info)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt splicesiteinfo -seqfile " +
             "#{$gthtestdata}ngasp32/training_regions.fa " +
             "#{$gthtestdata}ngasp/training_ref2.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/training_ref2.splicesiteinfo"
  end

  Name "nGASP (gthbssmtrain & gthbssmbuild)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}gthbssmtrain -filtertype CDS -seed 2040532791 -seqfile " +
             "#{$gthtestdata}ngasp32/training_regions.fa " +
             "#{$gthtestdata}ngasp/training_ref2.gff3", :maxtime => 600
    run_test "#{$bin}gthbssmbuild -bssmfile ngasp.bssm -datapath " +
             "training_data -gtdonor -agacceptor"
    run "diff --strip-trailing-cr -r ngasp.bssm #{$gthtestdata}ngasp/ngasp.bssm"
  end

  Name "nGASP (intron length distribution)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt stat -intronlengthdistri " +
             "#{$gthtestdata}ngasp/training_ref2.gff3"
    run "diff --strip-trailing-cr #{$last_stdout} " +
        "#{$gthtestdata}ngasp/training_ref2.intronlengthdistri"
  end

# XXX: this test never worked
=begin
  Name "nGASP (genometools triad)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt filter -mingenescore 0.85 " +
             "#{$gthtestdata}ngasp/ngasp.inter.gff3 | " +
             "#{$bin}../../genometools/bin/gt csa | " +
             "#{$bin}../../genometools/bin/gt cds -seqfile " +
             "#{$gthtestdata}ngasp32/test_regions.fa " +
             "-startcodon yes -o ngasp.gff3", :maxtime => 300
    run "diff --strip-trailing-cr ngasp.gff3 #{$gthtestdata}ngasp/ngasp.gff3"
  end
=end

  Name "nGASP (sensitivity)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt eval " +
             "#{$gthtestdata}ngasp/test_ref1.gff3 " +
             "#{$gthtestdata}ngasp/ngasp.gff3", :maxtime => 120
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/ngasp.sensitivity"
  end

  Name "nGASP (specificity)"
  Keywords "gth ngasp diss"
  Test do
    run_test "#{$bin}../../genometools/bin/gt eval " +
             "#{$gthtestdata}ngasp/test_ref2.gff3 " +
             "#{$gthtestdata}ngasp/ngasp.gff3", :maxtime => 300
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/ngasp.specificity"
  end
end
