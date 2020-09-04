#!/usr/bin/env ruby

require_relative 'citation'
require 'date'

def expand_macro(text)
  ["Gene","Species","Program"].each do |w|
    text = text.gsub(/\\#{w}\{([^\}]*)\}/,"<i>\\1</i>")
  end
  return text
end

def expand_author(author)
  return author.gsub(/(et\. al)\.$/,"<i>\\1</i>.")
end

journalnames = [
  "Nature Genetics",
  "Nature",
  "Nature Biotechnology",
  "Nature Communications",
  "Science",
  "Genome Research",
  "Genome Biology",
  "Nucl. Acids Res.",
  "Bioinformatics",
  "PLOS one",
  "Gene",
  "Genetics",
  "Plant Biology",
  "Plant Molecular Biology",
  "Journal of Plant Physiology",
  "Molecular Plant Pathology",
  "Rice",
  "The Plant Cell",
  "The Plant Journal",
  "Front Plant Sci.",
  "Genomics",
  "BMC Genomics",
  "BMC Genetics",
  "BMC Plant Biology",
  "BMC Bioinformatics",
  "Animal Genetics",
  "Current Bioinformatics",
  "Database--the journal of biological databases and curation",
  "Microbial Cell Factories",
  "Molecular Genetics and Genomics",
  "Philos Trans A Math Phys Eng Sci.",
  "International Journal of Genomics"]

journalrank = Hash.new()
journalnames.each_with_index do |name,idx|
  journalrank[name] = idx
end

def compare_rank(a,b,journalrank)
  if not journalrank.has_key?(a.journal)
    STDERR.puts "rank of \"#{a.journal}\" not defined"
    exit 1
  end
  if not journalrank.has_key?(b.journal)
    STDERR.puts "rank of \"#{b.journal}\" not defined"
    exit 1
  end
  ra = journalrank[a.journal]
  rb = journalrank[b.journal]
  if ra < rb
    return -1
  elsif ra > rb
    return 1
  elsif a.year < b.year
    return 1
  elsif a.year > b.year
    return -1
  else
    return 0
  end
end

puts <<'HEADER'
<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<!--http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd-->
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>GenomeThreader Gene Prediction Software</title>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="description" content="The GenomeThreader gene prediction software computes gene structure predictions using a similarity-based approach where additional cDNA/EST and/or protein sequences are used to predict gene structures via spliced alignments." />
<meta name="keywords" content="gene structure prediction, gene prediction, cDNA alignment, EST alignment, protein alignment, cDNA mapping, EST mapping, protein mapping, spliced alignment, consensus spliced alignment, genome annotation, bioinformatics, computational biology" />
<link rel="stylesheet" type="text/css" href="style.css" />
</head>
<body>

<h1><i>GenomeThreader</i> Gene Prediction Software</h1>

<div id="trialbox">
<ul>
  <li><a href="download.html">Download <i>GenomeThreader</i>!</a></li>
</ul>
</div>

<p>
<i>GenomeThreader</i> is a software tool to compute gene structure predictions.
The gene structure predictions are calculated using a similarity-based approach
where additional cDNA/EST and/or protein sequences are used to predict gene
structures via spliced alignments.

<i>GenomeThreader</i> was motivated by disabling limitations in
<a href="http://bioinformatics.iastate.edu/cgi-bin/gs.cgi">
<i>GeneSeqer</i></a>, a popular gene prediction program which is widely used
for plant genome annotation.
</p>
<h2>Features</h2>

<ul>
  <li>
    Intron Cutout Technique: <br/>
    The intron cutout technique allows to overcome the time and space
    limitations of the dynamic programming (DP) algorithms used in
    <a href="http://bioinformatics.iastate.edu/cgi-bin/gs.cgi"><i>GeneSeqer</i></a>,
    in particular, when applied to organisms containing long introns.
  </li>
  <li>
    Baysian Splice Site Models (BSSMs): <br/>
    With BSSMs it is possible to assign probabilities to GT donor, GC donor,
    and AG acceptor sites. This information is used in the DP to get the exact
    exon/intron boundaries right.
  </li>
  <li>
    Combination of cDNA/EST Based Spliced Alignments with Protein Based Spliced
    Alignments: <br/>
    After (spliced) aligning the supplied cDNAs/ESTs and protein sequences onto
    the genomic template, <i>GenomeThreader</i> computes consensus spliced
    alignments. Consensus spliced alignments combine several spliced alignments
    to resolve the complete gene structure and to uncover alternative splicing.
  </li>
  <li>
    Incremental Updates: <br/>
    When the used cDNA/EST or protein database is updated, a common approach
    was to redo the complete mapping. With <i>GenomeThreader</i>, you can combine
    newly computed spliced alignments with precomputed spliced alignments to
    quickly recompute consensus spliced alignments.
  </li>
  <li>
    XML: <br/>
    The additional <i>GenomeThreader</i> XML output conforms to our gthXML
    standard <a href="GenomeThreader.rng.txt">GenomeThreader.rng.txt</a>. With
    the included script XML2GFF.py, it is possible to convert gthXML output to the
    <a href="http://www.sanger.ac.uk/software/GFF">GFF</a> format.
    A variety of gthXML-specific tools can be found
    <a href="http://brendelgroup.org/mespar1/gthxml/">here</a>.
  </li>
  <li>
    gthDB: <br/>
    We also <a href="http://brendelgroup.org/mespar1/gthxml/">provide</a>
    a schema and load script for gthDB, which permits storage
    and query of <i>GenomeThreader</i> output in a relational format.
  </li>
</ul>
<p>
References have been omitted for brevity; you can find them and more details on
the implementation in the <i>GenomeThreader</i>
<a href="doc/GreBreSpaKur2005.pdf">paper</a>.

How to take advantage of these features and many more is described in depth in
the <i>GenomeThreader</i> <a href="doc/gthmanual.pdf">manual</a>.
Please consult the <a href="faq.html">FAQ</a> page for frequently asked
questions.

All mentioned files and scripts are also part of the <i>GenomeThreader</i>
distribution (see below).
</p>
<h2>Availability</h2>
<p>
<i>GenomeThreader</i> is available free of charge.
You can <a href="download.html">download</a> a copy.
</p>
<h2>Examples</h2>
<ul>
  <li>
    Evaluation <a href="gthcases/softeng.html">cases</a> described in Gremme et
    al. 2005 (see below)
  </li>
  <li>
    A 16.6Kb rice gene structure tractable with <i>GenomeThreader</i> (using
    both an <a href="gthcases/biggene.gth.cut.txt">intron cutout</a> technique
    and <a href="gthcases/biggene.gth.nocut.txt">without</a>), but beyond
    <a href="gthcases/biggene.gsq.txt"><i>GeneSeqer</i></a>'s limitations.
  </li>
  <li>
    A 125Kb intron-containing human
    <a href="gthcases/bigintron.gth.txt">gene structure</a>.
  </li>
  <li>
    Small samples of gzip'ed
    <a href="gthcases/small_demo.gth.out.gz">plain text</a> and
    <a href="gthcases/small_demo.gthxml.out.gz">XML</a>
    GenomeThreader output.
  </li>
</ul>

<h2>Users</h2>
<p>
The following sites use <i>GenomeThreader</i>. This list is not intended to be
comprehensive.
</p>
<ul>
  <li>
    MIPS (Munich Information Center for Protein Sequences),
    <a href="http://www.helmholtz-muenchen.de/en/ibis">Institute of Bioinformatics and
    Systems Biology</a>
    (for plant genome annotation)
  </li>
  <li>
    <a href="http://www.cosmoss.org/">University of Freiburg,
    Plant Biotechnology</a> (to annotate <i>Physcomitrella patens</i>)
  </li>
  <li>
    <a href="http://www.plantgdb.org/">PlantGDB</a>
  </li>
  <li>
    <a href="http://waksman.rutgers.edu/">Waksman Institute</a>, Rutgers
    University
  </li>
  <li>
    <a href="http://sgn.cornell.edu/">SOL Genomics Network (SGN)</a>, Cornell
    University
  </li>
  <li>
    <a href="http://bioinformatics.psb.ugent.be/">Bioinformatics and
      evolutionary genomics division, VIB</a>, Gent University
  </li>
</ul>

<h2>Citations</h2>
<p>
<a name="citations"></a>
Here are the most important publications citing <i>GenomeThreader</i> (sorted by Journal)
</p>
<ol>
HEADER

citation_list = Array.new()
citation_iterator("index.xml") do |citation|
  citation_list.push(citation)
end
cite_stat = Hash.new() {0}

citation_list.sort {|a,b| compare_rank(a,b,journalrank)}.each do |citation|
  puts "<li>"
  puts "  #{expand_author(citation.author)}"
  puts "  <a href=\"#{citation.url}\">"
  puts "  #{expand_macro(citation.title)}</a>,"
  puts "  <i>#{citation.journal}</i>"
  print "  <b>#{citation.volume}</b>"
  if not citation.number.nil?
    if citation.number.match(/^[0-9]*$/)
      print "(#{citation.number})"
    else
      print " #{citation.number}"
    end
  end
  if not citation.pages.nil?
    print ":#{citation.pages}"
  end
  puts ", #{citation.year}."
  puts "</li>"
  cite_stat[citation.year] += 1
end

puts <<'FOOTER'
</ol>
<p>
If I missed a publication which cites <i>GenomeThreader</i>, please contact
<a href="mailto:gordon@gremme.org">me</a>.
</p>

<h2>Developers</h2>
<p>
<i>GenomeThreader</i> is being actively developed by the following individuals:
</p>
<ul>
  <li>
    <a href="http://gremme.org/">Gordon Gremme</a> (primary developer)
  </li>
  <li>
    <a href="http://www.zbh.uni-hamburg.de/en/prof-dr-stefan-kurtz.html">Stefan Kurtz</a> (<a href="http://www.vmatch.de"><i>Vmatch</i></a>,
    <a href="http://www.vmatch.de"><i>Mkvtree</i></a>, libkurtz)
  </li>
  <li>
    <a href="http://brendelgroup.org/group/volker.php">Volker Brendel</a> (BSSMs, conceptual ideas)
  </li>
  <li>
    <a href="http://brendelgroup.org/mespar1/">Michael E Sparks</a>
    (BSSMs, gthXML, gthDB)
  </li>
</ul>

<h2>Publications</h2>
<p>
Please cite the following article in publications about research using
<i>GenomeThreader</i>:
</p>
<ul>
  <li>
    G. Gremme, V. Brendel, M.E. Sparks, and S. Kurtz.
    <a href="doc/GreBreSpaKur2005.pdf">Engineering a software tool for gene
    structure prediction in higher organisms.</a> <i>Information and Software
    Technology</i>, <b>47</b>(15):965-978, 2005
  </li>
</ul>
<p>
For in-depth information about <i>GenomeThreader</i> please refer to the
following dissertation:
</p>
<ul>
  <li>
  G. Gremme.  <a href="http://ediss.sub.uni-hamburg.de/volltexte/2013/6237/pdf/Dissertation.pdf">Computational Gene Structure Prediction.</a> Ph.D. thesis, University of Hamburg, 2012
  </li>
</ul>
FOOTER

today = Date.today
puts <<FOOTER
<div id="footer">
Copyright &copy; 2003-#{today.year} <a href="mailto:gordon@gremme.org">
Gordon Gremme</a>. Last update: #{today}
</div>
FOOTER

totalcitations = 0
cite_stat.sort.each do |k,v|
  totalcitations += v
  STDERR.puts "#{k}\t#{v}\t#{totalcitations}"
end
