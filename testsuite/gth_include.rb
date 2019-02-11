Name "gth"
Keywords "gth"
Test do
  run_test("#{$bin}gth", :retval => 1)
  grep($last_stderr, /option "-genomic" is mandatory/);
end

Name "gth -help"
Keywords "gth"
Test do
  run_test "#{$bin}gth -help"
  run "grep -v ^Usage: #{$last_stdout}"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/gth_-help.out"
end

Name "gth -help+"
Keywords "gth"
Test do
  run_test "#{$bin}gth -help+"
  run "grep -v ^Usage: #{$last_stdout} | tr -d '\r'"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/gth_-help+.out"
end

Name "gth -version"
Keywords "gth"
Test do
  run_test "#{$bin}gth -help"
  grep($last_stdout, /gth/);
end

Name "gth regression test (assertion in gthsafilter)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/sctg_21447.fasta -cdna #{$testdata}gth/sctg_21447.est -maskpolyatails -minalignmentscore 0.95 -undetcharweight -2.0"
end

Name "gth regression test (-gff3out -intermediate)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -gff3out -intermediate"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/standard_inter.gff3"
end

Name "gth regression test (-gff3out -intermediate -fastdp)"
Keywords "gth fastdp"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -gff3out -intermediate -fastdp"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/standard_inter.gff3"
end

Name "gth regression test (-gff3out -skipalignmentout)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -gff3out -skipalignmentout"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/standard_skip.gff3"
end

Name "gth regression test (-gff3out -skipalignmentout -fastdp)"
Keywords "gth fastdp"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -gff3out -skipalignmentout -fastdp"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/standard_skip.gff3"
end

Name "gth regression test (-gff3out -intermediate -o)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -gff3out -intermediate -o output.gff3"
  run "diff --strip-trailing-cr output.gff3 #{$testdata}gth/standard_inter.gff3"
end

Name "gth regression test (-gff3out -skipalignmentout -o)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -gff3out -skipalignmentout -o output.gff3"
  run "diff --strip-trailing-cr output.gff3 #{$testdata}gth/standard_skip.gff3"
end

Name "gthconsensus regression test (-gff3out -intermediate)"
Keywords "gth gthconsensus"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -intermediate -xmlout -o test.inter"
  run_test "#{$bin}gthconsensus -intermediate -gff3out test.inter"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/standard_inter.gff3"
end

Name "gthconsensus regression test (-gff3out -skipalignmentout)"
Keywords "gth gthconsensus"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -intermediate -xmlout -o test.inter"
  run_test "#{$bin}gthconsensus -skipalignmentout -gff3out test.inter"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/standard_skip.gff3"
end

Name "gthconsensus regression test (-gff3out -intermediate -o)"
Keywords "gth gthconsensus"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -intermediate -xmlout -o test.inter"
  run_test "#{$bin}gthconsensus -intermediate -gff3out -o output.gff3 test.inter"
  run "diff --strip-trailing-cr output.gff3 #{$testdata}gth/standard_inter.gff3"
end

Name "gthconsensus regression test (-gff3out -skipalignmentout -o)"
Keywords "gth gthconsensus"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -intermediate -xmlout -o test.inter"
  run_test "#{$bin}gthconsensus -skipalignmentout -gff3out -o output.gff3 test.inter"
  run "diff --strip-trailing-cr output.gff3 #{$testdata}gth/standard_skip.gff3"
end

Name "gth regression test (-exondistri)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -species human -exondistri"
  grep($last_stdout, /length distribution of all exons/)
end

Name "gth regression test (-introndistri)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -species human -introndistri"
  grep($last_stdout, /length distribution of all introns/)
end

Name "gth regression test (-refseqcovdistri)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -species human -refseqcovdistri"
  grep($last_stdout, /reference sequence coverage distribution/)
end

Name "gth regression test (all distributions)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -species human -exondistri -introndistri -refseqcovdistri"
  grep($last_stdout, /length distribution of all exons/)
  grep($last_stdout, /length distribution of all introns/)
  grep($last_stdout, /reference sequence coverage distribution/)
end

Name "gth regression test (-skipindexcheck)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -species human -createindicesonly"
  run_test "#{$bin}gth -genomic #{$testdata}gth/genomic_sequence_single_exon.fas -cdna #{$testdata}gth/cdna_sequence_single_exon.fas -species human -skipindexcheck | grep -v -e '^\\$' -e 'file='"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/single_exon.out"
end

Name "gth regression test (index reconstruction)"
Keywords "gth"
Test do
  # copy input files to test directory
  FileUtils.cp("#{$testdata}gth/genomic_sequence_single_exon.fas", ".")
  FileUtils.cp("#{$testdata}gth/cdna_sequence_single_exon.fas", ".")
  # run gth
  run_test "#{$bin}gth -genomic genomic_sequence_single_exon.fas -cdna cdna_sequence_single_exon.fas -species human"
  # make sure index files have been created & save mtimes
  mtimes = {}
  file_exists("cdna_sequence_single_exon.fas.dna.al1", mtimes)
  file_exists("cdna_sequence_single_exon.fas.dna.des", mtimes)
  file_exists("cdna_sequence_single_exon.fas.dna.ois", mtimes)
  file_exists("cdna_sequence_single_exon.fas.dna.prj", mtimes)
  file_exists("cdna_sequence_single_exon.fas.dna.sds", mtimes)
  file_exists("cdna_sequence_single_exon.fas.dna.tis", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.al1", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.bck", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.bwt", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.des", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.lcp", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.llv", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.ois", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.prj", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.sds", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.sti1", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.suf", mtimes)
  file_exists("genomic_sequence_single_exon.fas.dna.tis", mtimes)
  # run gth again
  run_test "#{$bin}gth -genomic genomic_sequence_single_exon.fas -cdna cdna_sequence_single_exon.fas -species human"
  # make sure index files have have not changed
  cmp_mtimes("cdna_sequence_single_exon.fas.dna.al1", mtimes)
  cmp_mtimes("cdna_sequence_single_exon.fas.dna.des", mtimes)
  cmp_mtimes("cdna_sequence_single_exon.fas.dna.ois", mtimes)
  cmp_mtimes("cdna_sequence_single_exon.fas.dna.prj", mtimes)
  cmp_mtimes("cdna_sequence_single_exon.fas.dna.sds", mtimes)
  cmp_mtimes("cdna_sequence_single_exon.fas.dna.tis", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.al1", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.bck", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.bwt", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.des", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.lcp", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.llv", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.ois", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.prj", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.sds", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.sti1", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.suf", mtimes)
  cmp_mtimes("genomic_sequence_single_exon.fas.dna.tis", mtimes)
end

Name "gth regression test (-createindicesonly & -noautoindex)"
Keywords "gth"
Test do
  # copy input files to test directory
  FileUtils.cp("#{$testdata}gth/genomic_sequence_single_exon.fas", ".")
  FileUtils.cp("#{$testdata}gth/cdna_sequence_single_exon.fas", ".")
  # run gth to construct index
  run_test "#{$bin}gth -createindicesonly " +
           "-genomic genomic_sequence_single_exon.fas " +
           "-cdna cdna_sequence_single_exon.fas " +
           "-inverse" # option is necessary to use -frompos later on!
  # run gth without automatic index construction
  run_test "#{$bin}gth -skipindexcheck " +
           "-genomic genomic_sequence_single_exon.fas " +
           "-cdna cdna_sequence_single_exon.fas " +
           "-frompos 8000 -topos 10000"
end

Name "gth regression test (coverage assertion)"
Keywords "gth"
Test do
  run_test "#{$bin}/gth -genomic #{$testdata}gth/mydb -cdna #{$testdata}gth/myqy -minmatchlen 10 -seedlength 8 -mincutoffs"
end

Name "gth -gff3out -intermediate (long intron files)"
Keywords "gth"
Test do
  FileUtils.cp("#{$testdata}gth/long_intron_genomic_1.fas", ".")
  FileUtils.cp("#{$testdata}gth/long_intron_genomic_2.fas", ".")
  FileUtils.cp("#{$testdata}gth/long_intron_cdna.fas", ".")
  run_test "#{$bin}/gth -genomic long_intron_genomic_1.fas " +
           "long_intron_genomic_2.fas -cdna long_intron_cdna.fas " +
           "-introncutout -gff3out -intermediate", :maxtime => 300
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/long_intron_inter.gff3"
end

Name "gth -gff3out -skipalignmentout (long intron files)"
Keywords "gth"
Test do
  FileUtils.cp("#{$testdata}gth/long_intron_genomic_1.fas", ".")
  FileUtils.cp("#{$testdata}gth/long_intron_genomic_2.fas", ".")
  FileUtils.cp("#{$testdata}gth/long_intron_cdna.fas", ".")
  run_test "#{$bin}/gth -genomic long_intron_genomic_1.fas " +
           "long_intron_genomic_2.fas -cdna long_intron_cdna.fas " +
           "-introncutout -gff3out -skipalignmentout", :maxtime => 300
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/long_intron_final.gff3"
end

Name "gth -gff3out -frompos (correct positions in sequence region)"
Keywords "gth"
Test do
  run_test "#{$bin}/gth -genomic " +
           "#{$testdata}gth/genomic_sequence_single_exon.fas " +
           "-cdna #{$testdata}gth/cdna_sequence_single_exon.fas " +
           "-gff3out -skipalignmentout -frompos 8001 -topos 10000"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/frompos.gff3"
end

Name "gth -gff3out -intermediate (correct reverse positions)"
Keywords "gth"
Test do
  run_test "#{$bin}gth -gff3out -intermediate -genomic " +
           "#{$testdata}gth/genomic_sequence_single_exon.fas -cdna " +
           "#{$testdata}gth/cdna_sequence_single_exon.fas -r"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/reverse.gff3"
end

Name "gth arabidopsis"
Keywords "gth arabidopsis"
Test do
  run_test("#{$bin}gth -species arabidopsis -genomic " +
           "#{$testdata}U89959_genomic.fas -cdna " +
           "#{$testdata}U89959_ests.fas -intermediate -xmlout -o test.inter",
           :maxtime => 300)
  run_test "#{$bin}gthconsensus -o test_sas.gff3 -gff3out -intermediate " +
           "test.inter"
  run      "diff --strip-trailing-cr test_sas.gff3 #{$testdata}U89959_sas.gff3"
  run_test "#{$bin}gthconsensus -o test_csas.gff3 -gff3out -skipalignmentout " +
           "test.inter"
  run_test "#{$bin}../../genometools/bin/gt eval test_csas.gff3 #{$testdata}U89959_csas.gff3"
  run      "diff --strip-trailing-cr #{$last_stdout} #{$testdata}U89959_csas.eval"
end

Name "gth illegal character in protein file"
Keywords "gth"
Test do
  run_test("#{$bin}gth -genomic #{$testdata}U89959_genomic.fas " +
           "-protein #{$testdata}gth/illegal_protein.fas", :retval => 1)
  grep($last_stderr, /Illegal character 'x' in file/);
end

Name "align_dna"
Keywords "gth align_dna"
Test do
  run_test "#{$bin}align_dna " +
           "#{$testdata}gth/genomic_sequence_single_exon.fas " +
           "#{$testdata}gth/cdna_sequence_single_exon.fas"
  run "grep -v -e '^Genomic Template:' -e 'EST Sequence:' #{$last_stdout}"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}gth/align_dna.out"
end

Name "gth -gff3descranges (single exon)"
Keywords "gth gff3descranges"
Test do
  run_test "#{$bin}gth -gff3descranges -gff3out -intermediate " +
           "-genomic #{$testdata}genomic_sequence_single_exon_descrange.fas " +
           "-cdna #{$testdata}cdna_sequence_single_exon.fas"
  run "diff --strip-trailing-cr #{$last_stdout} #{$testdata}single_exon_descrange.gff3"
end

if $gthtestdata then
  Name "gth regression test (-proteinsmap TransProt11)"
  Keywords "gth"
  Test do
    run_test "#{$bin}gth -species maize -genomic #{$gthtestdata}156255364.fsa.raw -protein #{$gthtestdata}OSpep_multi -proteinsmap TransProt11"
  end

  Name "gth regression test (multi vs. single protein file)"
  Keywords "gth"
  Test do
    run_test "#{$bin}gth -species maize -genomic #{$gthtestdata}156255364.fsa.raw -protein #{$gthtestdata}OSpep_multi"
    run "grep -v -e '^\\\$' -e '^Protein Sequence:' #{$last_stdout}"
    FileUtils.mv($last_stdout, "multi.txt")
    run_test "#{$bin}gth -species maize -genomic #{$gthtestdata}156255364.fsa.raw -protein #{$gthtestdata}OSpep_single"
    run "grep -v -e '^\\\$' -e '^Protein Sequence:' #{$last_stdout}"
    FileUtils.mv($last_stdout, "single.txt")
    run "diff --strip-trailing-cr multi.txt single.txt"
  end

  Name "gth (detect small initial exon)"
  Keywords "gth gencode"
  Test do
    run_test "#{$bin}gth -species human -exdrop 1 -minmatchlen 12 " +
             "-seedlength 12 -fragweightfactor 0.5 -introncutout " +
             "-enrichchains -scoreminexonlen 1 -detectsmallexons -gff3out " +
             "-skipalignmentout -genomic #{$gthtestdata}gencode/ENr233.fa.gz " +
             "-cdna #{$gthtestdata}gencode/small_initial_exon_cdna.fas"
    run_test "#{$bin}../../genometools/bin/gt eval " +
             "#{$gthtestdata}gencode/small_initial_exon.gff3 #{$last_stdout}"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}gencode/small_initial_exon.eval"
  end

  Name "gth (detect small terminal exon)"
  Keywords "gth gencode"
  Test do
    run_test "#{$bin}gth -species human -exdrop 1 -minmatchlen 12 " +
             "-seedlength 12 -fragweightfactor 0.5 -introncutout " +
             "-enrichchains -scoreminexonlen 1 -detectsmallexons -gff3out " +
             "-skipalignmentout -genomic #{$gthtestdata}gencode/ENm005.fa.gz " +
             "-cdna #{$gthtestdata}gencode/small_terminal_exon_cdna.fas"
    run_test "#{$bin}../../genometools/bin/gt eval " +
             "#{$gthtestdata}gencode/small_terminal_exon.gff3 #{$last_stdout}"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}gencode/small_terminal_exon.eval"
  end

  Name "gth (detect small exon bug 1)"
  Keywords "gth gencode"
  Test do
    run_test "#{$bin}gth -species human -exdrop 2 -minmatchlen 12 " +
             "-maskpolyatails -noicinintroncheck " +
             "-seedlength 12 -fragweightfactor 0.5 -introncutout " +
             "-enrichchains -scoreminexonlen 1 -detectsmallexons -gff3out " +
             "-skipalignmentout -genomic #{$gthtestdata}gencode/ENm001.fa.gz " +
             "-cdna #{$gthtestdata}gencode/detect_small_exon_bug_1_cdna.fas"
  end

  Name "gth (detect small exon bug 2)"
  Keywords "gth gencode bug"
  Test do
    run_test("#{$bin}gth -fragweightfactor 0.5 -exdrop 1 -seedlength 12 " +
             "-minmatchlen 12 -enrichchains -scoreminexonlen 1 " +
             "-noicinintroncheck -detectsmallexons -introncutout " +
             "-maskpolyatails -species human -genomic " +
             "#{$gthtestdata}gencode/ENm003.fa.gz " +
             "-cdna #{$gthtestdata}gencode/detect_small_exon_bug_2_cdna.fas",
             :maxtime => 900)
  end

  Name "gth (detect small exon bug 3)"
  Keywords "gth gencode"
  Test do
    run_test("#{$bin}gth -fragweightfactor 0.5 -exdrop 1 -seedlength 12 " +
             "-minmatchlen 12 -enrichchains -scoreminexonlen 1 " +
             "-noicinintroncheck -detectsmallexons -introncutout " +
             "-maskpolyatails -species human -genomic " +
             "#{$gthtestdata}gencode/ENm004.fa.gz " +
             "-cdna #{$gthtestdata}gencode/detect_small_exon_bug_3_cdna.fas",
             :maxtime => 300)
  end

  Name "gth (detect small exon bug 4)"
  Keywords "gth gencode"
  Test do
    run_test("#{$bin}gth -fragweightfactor 0.5 -exdrop 1 -seedlength 12 " +
             "-minmatchlen 12 -enrichchains -scoreminexonlen 1 " +
             "-noicinintroncheck -detectsmallexons -introncutout " +
             "-maskpolyatails -species human -genomic " +
             "#{$gthtestdata}gencode/ENm004.fa.gz " +
             "-cdna #{$gthtestdata}gencode/detect_small_exon_bug_4_cdna.fas",
             :maxtime => 300)
  end

  Name "gth (detect small exon bug 5)"
  Keywords "gth gencode"
  Test do
    run_test("#{$bin}gth -fragweightfactor 0.5 -exdrop 1 -seedlength 12 " +
             "-minmatchlen 12 -enrichchains -scoreminexonlen 1 " +
             "-noicinintroncheck -detectsmallexons -introncutout " +
             "-maskpolyatails -species human -genomic " +
             "#{$gthtestdata}gencode/ENm014.fa.gz " +
             "-cdna #{$gthtestdata}gencode/detect_small_exon_bug_5_cdna.fas",
             :maxtime => 600)
  end

  Name "gth -gff3descranges (nGASP)"
  Keywords "gth gff3descranges ngasp"
  Test do
    run_test("#{$bin}gth -genomic " +
             "#{$gthtestdata}ngasp/training_regions.fa.gz " +
             "-protein #{$gthtestdata}ngasp32/protein_small.fa -intermediate " +
             "-xmlout -o out.gz -gzip", :maxtime => 900)
    grep $last_stderr, /does not end with a stop amino acid/
    run_test"#{$bin}gthconsensus -gff3out -intermediate out.gz"
    run "diff --strip-trailing-cr #{$last_stdout} #{$gthtestdata}ngasp/protein_small_default.gff3"
    run_test"#{$bin}gthconsensus -gff3descranges -gff3out -intermediate out.gz"
    run "diff --strip-trailing-cr #{$last_stdout} " +
        "#{$gthtestdata}ngasp/protein_small_descranges.gff3"
  end

  1.upto(14) do |i|
    Name "gth -fastdp protein test #{i}"
    Keywords "gth fastdp"
    Test do
      run_test "#{$bin}gth -fastdp " +
               "-genomic #{$gthtestdata}protein/proteintest#{"%02d" % i}_dna.fna " +
               "-protein #{$gthtestdata}protein/proteintest#{"%02d" % i}_protein.fna"
    end
  end
end
