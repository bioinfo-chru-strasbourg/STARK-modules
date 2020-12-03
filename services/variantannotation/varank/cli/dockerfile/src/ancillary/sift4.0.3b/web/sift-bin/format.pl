#!/usr/bin/perl
#  format.pl <pid>
#	Formats SIFT output in ../tmp/<pid>.*

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.


$| = 1;
require 'SIFT_subroutines.pm';

#	Set file permissions to rw-rw----
system("umask 006");
$bin = ".";
$tmp = "/opt/www/sift/tmp";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
#print "<TITLE>SIFTING Results</TITLE>\n";
print "<body bgcolor=white>\n";

if ( @ARGV < 1 )
{
   print "USAGE: format.pl <pid>.\n";
   exit(-1);
}
else { $pid = $ARGV[0]; }

$queryseqfile = $tmp . "/$pid.seq";
$out = $tmp . "/$pid.siftresults"; # where sifting results are printed
$comments = $tmp . "/$pid.comments";

### substitution file #######
$subst_file = $tmp . "/$pid.substfile";

######  printing results #########

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font color=blue>T</font> results</center></H1></A><BR>\n";

#>>>>Generates message in Apache error_log if none of these
@error_files = `ls $tmp/$pid.*.error`;
if (@error_files > 0) {
	# there's an error
	print "<b>ERROR!</b>  Please <A HREF=\"/sift-bin/contact.pl\">email us</A> with your sequence if you do not understand the error.<BR><BR>";
	foreach $errorfile (@error_files) {
#print "$errorfile\n";
		open (OUT, "<$errorfile");
		while ($_ = <OUT>) {print;}
		close (OUT);
	}
	if (open (OUT, "$tmp/$pid.seq.query.globalX")) {
		close (OUT);
		print "<BR>";
		system ("sed s/X/-/g $tmp/$pid.seq.query.globalX > $tmp/$pid.aln");
		print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$pid.aln+Alignment+PRE\">";
		print "PSIBLAST search results</A><BR>\n";
	}
#>>>	system ("rm -f $tmp/$pid.*");   <<<<NO!
	exit (-1);
}
 
print "</PRE>";

# SIFT successful -- OUTPUT RESULTS
$commentscsh = $tmp . "/$pid.commentscsh"; 
$psiblastout = $tmp . "/$pid" . ".alignedfasta";
$seqno = `grep ">" $psiblastout | wc -l`;
$seqno--; # subtract the QUERY sequence

print "<b>$seqno</b> sequences were selected to be closely related to your query sequence.<BR>\n";
 
print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$pid.msf+Alignment+PRE\" TARGET=_blank>";
print "PSIBLAST alignment of submitted sequences</A><BR>\n";

print "<A HREF=\"/sift-bin/catfile_notitle.csh?$tmp/$pid.alignedfasta+PRE\" TARGET=_blank>Alignment in FASTA format </A> (for modification)<BR>"; 

print_psiblast_alignment_desc();
print "<i>Please check the sequences that have been chosen.  If the sequences are too diverged from your query or the alignment is questionable, we suggest you modify the fasta-formatted file above and <A HREF=\"/sift/SIFT_aligned_seqs_submit.html\" TARGET=_blank>resubmit</A>.</i><BR>"; 
print "<BR>";


$no_table_files = `ls $tmp/$pid.aatable* |grep -v html | wc -l`;

tolerated_nottolerated_aminoacidtable($tmp,$no_table_files, $pid);

$commentscsh = $tmp . "/$pid.commentscsh";
$loc = rindex ($out, "/");
$file = substr ($out, $loc+1, 100);
print "<BR>\n";


$loc = rindex ($out, "/");
$file = substr ($out, $loc+1, 100);

print_score_desc ($tmp,$file);
print "<BR>\n";

print_predictions ($tmp,$file);
print "<BR>\n";

#print "<HR>";
#print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$pid.subblock+Blocks+PRE\">"
;
#print "Block(s)</A> extracted from the alignment of the query sequence and";
#print " the $seqno closely related sequences<BR>";
#print "<i>For wider blocks, submit the sequences to ";
#print "<A HREF=\"/blocks/make_blocks.html\">";
#print "Block Maker</A></i><BR>";
#print "<BR>\n";

print "<HR>";

#open (OUT, "<$tmp/$pid.alignedfasta");
#print "<PRE>";
#while ($_= <OUT>) {
#print;
#}
#close (OUT);

#open (OUT, "<$tmp/$pid.seq.query.selectedclumped.log");
#while ($_ = <OUT>) {
#print;
#}
#close (OUT);
#print "</PRE>";
#print "<HR>";

#errors or warnings that don't disrupt program
open (OUT, "<$tmp/$pid.err");
while ($_ = <OUT>) {
# don't want to print Jorja's messages, only my own.  so my error
# messages all have ERROR in the line, but don't print blastpgp error for
#reading checkpoing file
      if (/ERROR/ || /WARNING/ || /\*\*\*/ ) {
            unless (/blastpgp/) {
		print;
      		}
	}
}
close(OUT);

#comments
#open (OUT, "<$comments");
#while ($_ = <OUT>) {print;}
#close (OUT);

system ("rm -f $tmp/$pid.comments*");
system ("rm -f $tmp/$pid.tmp");

#system ("rm -f $tmp/$pid*");
#system ("rm -f $tmp/*$pid");

#print "</PRE>";
#print "</PRE><A HREF=\"#top\">[return to top]</A><P>";
#-------------------------------------------------------------------------
exit (0);
