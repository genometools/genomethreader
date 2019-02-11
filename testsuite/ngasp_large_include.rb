if $gthtestdata then
  files="-genomic #{$gthtestdata}ngasp32/test_regions.fa " +
        "-cdna #{$gthtestdata}ngasp32/rna_matches.fa " +
        "-protein #{$gthtestdata}ngasp32/protein_matches.fa"

  opts="-bssm #{$gthtestdata}ngasp/ngasp -dpminintronlen 10 " +
       "-gcmaxgapwidth 22000 -enrichchains -scoreminexonlen 1 -maskpolyatails"

  Name "nGASP (gene prediction, complete)"
  Keywords "gth ngasp diss md5ids"
  Test do
    run "#{$bin}gth #{files} #{opts} -minalignmentscore 0.85 -startcodon yes " +
        "-skipalignmentout -gff3out -md5ids -o ngasp.gff3", :maxtime => 3600
    run "diff --strip-trailing-cr ngasp.gff3 #{$gthtestdata}ngasp/ngasp.gff3"
  end

  Name "nGASP (gene prediction, intermediate)"
  Keywords "gth ngasp diss md5ids"
  Test do
    run "#{$bin}gth #{files} #{opts} -intermediate " +
        "-xmlout -gzip -o ngasp.inter.gz", :maxtime => 3600
    run "#{$bin}gthconsensus -minalignmentscore 0.85 -startcodon yes " +
        "-skipalignmentout  -gff3out -md5ids " +
        "-o ngasp.inter.gff3 ngasp.inter.gz ", :maxtime => 1200
    run "diff --strip-trailing-cr ngasp.inter.gff3 #{$gthtestdata}ngasp32/ngasp.inter.gff3"
  end
end
