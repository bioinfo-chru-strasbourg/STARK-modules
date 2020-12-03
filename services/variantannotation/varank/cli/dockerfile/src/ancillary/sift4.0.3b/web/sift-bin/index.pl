#!/usr/local/bin/perl -w

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

use strict;
use Template;
use warnings; 
print "Content-type: text/html\n\n";

my $title = 'SIFT Home';
my $site = 'SIFT';
my $project_name = 'SIFT';
my @top_menu = (
                 {'link' => 'www/SIFT_help.html', 'menu_name' => 'Help'},
                 {'link' => 'sift-bin/contact.pl', 'menu_name' => 'Contact us'});

my @side_menu = ({'link' => '/', 'menu_name' => 'SIFT Home'},
                 {'link' => 'www/SIFT_help.html', 'menu_name' => 'Help'},
                 {'link' => 'sift-bin/contact.pl', 'menu_name' => 'Contact us'});

my $left_content = '<p><strong>Code &amp; Data</strong></p>'."\n";
$left_content .= '<p>SIFT for in-house use (version 3.0, released March 21, 2008):</p>'."\n";
$left_content .= '<p><a href="/blocks/license.html">Copyright</a> <a href="www/sift3.0.tar">code & exe (Sun, Linux)</a><br /><a href="/sift-bin/contact.pl">Report bugs</a></p>'."\n";
$left_content .= '<p><a href="www/SIFTing_databases.html">Prediction on human SNPs</a></p>';

my @right_content = ({'header' => 'Literature', 'content' => '<p>Our review "Predicting the Effects of Amino Acid Substitutions on Protein Function" in <em>Annual Review of Genomics and Human Genetics</em><br /><a href="www/annurev.genom.7.Ng_Henikoff.pdf">Chapter</a><br /><a href="www/annurev.genom.7.080505_SuppTable1.pdf">Supplementary Table 1</a></p>'},
                     {'header' => 'Referencing SIFT', 'content' => '<p>Predicting Deleterious Amino Acid Substitutions, Genome Res. 2001 May; 11(5): 863.874. <a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11337480&dopt=Abstract">(Link)</a></p><p>Accounting for Human Polymorphisms Predicted to Affect Protein Function, Genome Res. 2002 December; 436-446 <a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11875032&dopt=Abstract">(Link)</a></p><p>SIFT: predicting amino acid changes that affect protein function, Nucleic Acids Research, 2003, Vol. 31, No. 13 3812-3814 <a href="http://www.ncbi.nlm.nih.gov/pubmed/12824425">(Link)</a></p><p>Predicting the effects of coding non-synonymous variants on protein function using the SIFT algorithm, Nature Protocols 4, - 1073 - 1081 (2009) <a href="/www/nprot.2009.86.pdf">(Link)</a></p>'});

my $main_content = '<p><strong>SIFT</strong> predicts whether an <strong>amino acid substitution affects protein function</strong> based on sequence homology and the physical properties of amino acids. SIFT can be applied to naturally occurring <strong>nonsynonymous polymorphisms and laboratory-induced missense mutations.</strong></p>'."\n";
$main_content .= '<table class="contenttable">'."\n";
$main_content .= '  <tbody>'."\n";
$main_content .= '    <tr class="tableHeader">'."\n";
$main_content .= '      <td><p>Genome Tools</p></td>'."\n";
$main_content .= '      <td><p>Tool Description/Input</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowOdd">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_intersect_coding_submit.html">1. Restrict to coding variants</a></p></td>'."\n";
$main_content .= '      <td><p>Extract coding variants from a large list of genomic variants</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowEven">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_chr_coords_submit.html">2a. SIFT nonsynonymous single nucleotide variants (genome-scale)</a></p></td>'."\n";
$main_content .= '      <td><p>Provides SIFT predictions for variants with genome coordinates</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowOdd">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_chr_coords_indels_submit.html">2b. Classify coding indels (Insertion / Deletion polymorphisms)</a></p></td>'."\n";
$main_content .= '      <td><p>Indel variants in genome coordinates</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableHeader">'."\n";
$main_content .= '      <td><p>Batch Tools</p></td>'."\n";
$main_content .= '      <td><p>Tool Description/Input</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowOdd">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_pid_subst_all_submit.html">SIFT Batch Protein</a></p></td>'."\n";
$main_content .= '      <td><p>Multiple proteins identifiers</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowEven">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_dbSNP.html">SIFT dbSNP</a></p></td>'."\n";
$main_content .= '      <td><p>Multiple rsIDs (NCBI dbSNP)</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableHeader">'."\n";
$main_content .= '      <td><p>Single Protein Tools</p></td>'."\n";
$main_content .= '      <td><p>Tool Description/Input</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowOdd">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_BLink_submit.html">SIFT BLink</a></p></td>'."\n";
$main_content .= '      <td><p>Single protein (RefSeq ID /gi number)</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowEven">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_seq_submit2.html">SIFT Sequence</a></p></td>'."\n";
$main_content .= '      <td><p>Single protein (Query sequence)</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowOdd">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_related_seqs_submit.html">SIFT Related Sequences</a></p></td>'."\n";
$main_content .= '      <td><p>Single protein (Query with related sequences)</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '    <tr class="tableRowEven">'."\n";
$main_content .= '      <td><p><a href="www/SIFT_aligned_seqs_submit.html">SIFT Aligned Sequences</a></p></td>'."\n";
$main_content .= '      <td><p>Single protein (Query with aligned sequences)</p></td>'."\n";
$main_content .= '    </tr>'."\n";
$main_content .= '  </tbody>'."\n";
$main_content .= '</table>'."\n";
$main_content .= '<hr />'."\n";
$main_content .= '<p>Other prediction tools for amino acid substitutions: <a href="http://www.bork.embl-heidelberg.de/PolyPhen">PolyPhen</a>, <a href="http://mendel.stanford.edu/SidowLab/downloads/MAPP/index.html">MAPP</a>, <a href="http://www.snps3d.org">SNPs3D</a></p>'."\n";

my $tt = Template->new({
  ABSOLUTE => 1,
});
my $template_file = '/opt/www/common/perl_templates/3_column_fluid_width_full_screen.tpl';
my $vars = {
  title => $title,
  top_menu => \@top_menu,
  side_menu => \@side_menu,
  left_content => $left_content,
  right_content => \@right_content,
  main_content => $main_content,
  site => $site,
  project_name => $project_name,
  search_site => 'www',
};

$tt->process($template_file, $vars) || die $tt->error();
