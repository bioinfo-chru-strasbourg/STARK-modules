#!/usr/bin/perl

# getblock.pl
# 6/7/06 Added cat of $pros
# 2/18/07 blockmap from here instead of proweb
# 4/ 9/07 treeviewer from here instead of proweb

#	Set file permissions to rw-rw----
system("umask 006");

select(STDOUT); $| = 1;

print "Content-type: text/html\n\n";
print "<HTML><TITLE>Blocks Information</TITLE>\n";

if (@ARGV > 0 ) { $group = $ARGV[0]; }
else
{
   read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
   %names = &parse_query($QUERY_STRING);
   $group = $names{group};
}
#	Strip off any alpha chars at end
$group =~ s/[A-Z]$//;
#print "group=$group\n";
if ($group eq "")
{
   print "<P><B>ERROR:</B> Please enter a Blocks family (e.g. IPB001525).<P>\n";
   exit(-1);
}
$blocks_flag = 0;
$prints_flag = 0;
if ($group =~ m/^IPB/) { $blocks_flag = 1; }
elsif ($group =~ m/^IPR/) { $group =~ s/^IPR/IPB/; $blocks_flag = 1; }
elsif ($group =~ m/^PR/) { $prints_flag = 1; }
if ($blocks_flag == 0 && $prints_flag == 0)
{
   print "<P><B>ERROR:</B> $group is not a recognized Blocks family.<P>\n";
   exit(-1);
}

#print "<A NAME=top><H1>Get Block $group</H1></A>\n";
#print "<A NAME=top></A>\n";
print "<A NAME=top><H1 ALIGN=CENTER><IMG ALIGN=MIDDLE SRC=\"/blocks/icons/hutch_logo.gif\">Blocks Information for $group<IMG ALIGN=MIDDLE SRC=\"/blocks/icons/blocks.gif\"></H1></A>\n";


#=========================================================================
   
print "<P>\n";
if ( $blocks_flag )
{
   $ipr = $group; $ipr =~ s/^IPB/IPR/;
   $bdir = "../data-blocks";
   $maptype = "IPB";
}
elsif ( $prints_flag )
{
   $bdir = "../data-prints";
   $maptype = "PR";
}
$head = "$bdir/bhead";
$blks = "$bdir/blks/$group.blks";
$map = "$bdir/maps/$group.map";
$tree = "$bdir/trees/$group.tree";
$cob = "$bdir/cobs/$group.cob";
$pdbs = "$bdir/pdbs/$group.pdb";
$pros = "$bdir/pros/$group.pros";

#	Read the blocks into memory to get the DE and list of ACs
open (IN, "<$blks");
@blocks = <IN>;
close(IN);
$id = $de = "";
@acs = ();
foreach $in_rec (@blocks)
{
   if ($in_rec =~ m/^AC   /)
   {
     ($ac) = $in_rec =~ m/^AC   (\S+)/;
     $ac =~ s/;//g;
     @acs = (@acs, $ac);
   }
   elsif ($id eq "" && $in_rec =~ m/^ID   /)
   {
     ($id) = $in_rec =~ m/^ID   (\S+)/;
     $id =~ s/;//g;
   }
   elsif ($de eq "" && $in_rec =~ m/^DE   /)
   {
     ($de) = $in_rec =~ m/^DE   (.+)/;
   }
}
#----------------------------------------------------------------------
#  Start printing
print "<TITLE>Blocks Information for $group</TITLE>\n";
if ($#acs < 0)
{
   print "<P><B>ERROR:</B> No blocks found for $group.<P>\n";
   exit(-1);
}

print "<A NAME=top><H1>$group: $id</H1></A>\n";
print "<H2>$de</H2><P>\n";

print "<MENU>\n";
print "<UL>\n";
print "<LI> <A HREF=\"/blocks/help/blocks_format.html\">Introduction</A>\n";

foreach $ac (@acs)
{
   print "<LI> <A HREF=\"#", $ac, "\">Block number ", $ac, "</A>\n";
}

print "<P>\n";
if ($blocks_flag)
{
   print "<LI>";
   print "InterPro entry <A HREF=\"http://www.ebi.ac.uk/interpro/IEntry?ac=$ipr\"TARGET=\"getblock\"> $ipr</A>" ;
   print " (source of sequences used to make blocks)\n";
}
elsif ($prints_flag)
{
   print "<LI>";
   print "PRINTS Entry <A HREF=\"http://www.bioinf.man.ac.uk/cgi-bin/dbbrowser/PRINTS/DoPRINTS.pl?cmd_a=Display&qua_a=/Brief&fun_a=Text&qst_a=$group\"TARGET=\"getblock\"> $group</A>";
   print " (source of blocks)\n";
}

print "<P>\n";
print "<LI>Protein Sequences Used to Make Blocks.";
print "[<A HREF=\"/blocks-bin/catfile.sh?$pros\">Sequences in fasta format</A>]";

print "<P>\n";
print "<LI>Block Maps.";
#print "[<A HREF=\"http://www.proweb.org/proweb-bin/blockmap.cgi?name=$group&type=FAM&condensed=YES&dbtype=$maptype\"TARGET=\"getblock\">Graphical Map</A>]";
print "[<A HREF=\"/blocks-bin/blockmap.pl?name=$group&type=FAM&condensed=YES&dbtype=$maptype\"TARGET=\"getblock\">Graphical Map</A>]";
print "     [<A HREF=\"/blocks-bin/map.csh?$map\"TARGET=\"getblock\">Text Map</A>]";
print "     [<A HREF=\"/blocks-bin/catfile.sh?$map\"TARGET=\"getblock\">Map Positions</A>]";
print "     [<A HREF=\"/blocks/help/about_maps.html\"TARGET=\"getblock\">About Maps</A>]";

print "<P>\n";
print "<LI>Logos.";
print "[<A HREF=\"/blocks/help/about_logos.html\">About Logos</A>]";
print "<BR>Select display format: ";
if (-e "../logos/$group.gif")
{  
   print "<A HREF=\"/blocks/logos/$group.gif\">[GIF]</A> \n"; 
   print "<A HREF=\"/blocks-bin/convert.csh?../logos/$group.gif+pdf\">[PDF]</A> ";
}
else
{  
   print "<A HREF=\"/blocks-bin/logo.csh?$blks+gif\">[GIF]</A> ";
   print "<A HREF=\"/blocks-bin/logo.csh?$blks+pdf\">[PDF]</A> ";
}
print "<A HREF=\"/blocks-bin/logo.csh?$blks+ps\">[Postscript]</A> ";

#	Some groups with < 4 seqs don't have trees
if ( -s "$tree" )	# $ntree > 0 ???
{
print "<P>\n";
print "<LI>Tree from blocks alignment. ";
print "[<A HREF=\"/blocks/help/about_trees.html\">About Trees</A>]  ";
print "[<A HREF=\"http://www.proweb.org/treeinfo.html\">About ProWeb TreeViewer</A>]  ";
print "<BR>[<A HREF=\"/blocks-bin/drawgram.csh?$tree+dat\">Data</A>]  ";
#print " [<A HREF=\"http://www.proweb.org/proweb-bin/trees.cgi?family=$group\" TARGET=\"_blank\">ProWeb TreeViewer</A>] ";
print " [<A HREF=\"/blocks-bin/treeviewer.pl?family=$group\" TARGET=\"_blank\">TreeViewer</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?$tree+xbm\">XBitmap</A>]  ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?$tree+gif\">GIF</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?$tree+pdf\">PDF</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?$tree+ps\">Postscript</A>]  ";
#print "[<A HREF=\"/blocks-bin/drawgram.csh?$tree+new\">Newick</A>] ";
}

if (!(-z $pdbs))
{
   print "<P><LI><A HREF=\"/blocks-bin/pdbmast.pl?$pdbs+$group\">Structures</A>\n";
}

print "<P>";
print "<LI>Search blocks vs other databases:<BR><UL>\n";
print "<LI><A HREF=\"#cobbler\">COBBLER sequence</A>\n";
print " and BLAST searches";
print "     [<A HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]";

print "<LI>";
print "<A HREF=\"/blocks-bin/mast.sh?$blks+$group\" TARGET=\"getblock\">MAST Search</A>";
print " of all blocks vs a sequence database ";
print "     [<A HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">About MAST</A>]";

print "<LI>";
print "<A HREF=\"/blocks-bin/LAMA_search.sh?$blks\" TARGET=\"_blank\">LAMA search</A>";
print " of all blocks vs a blocks database  ";
print "     [<A HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]";
print "</UL>";

print "<P>";
print "<LI>";
print "<A HREF=\"/blocks-bin/codehop.sh?$blks\">CODEHOP</A>";
print " to design PCR primers from blocks";
print "     [<A HREF=\"/blocks/help/CODEHOP/tips.html\">About CODEHOP</A>]";

print "<P>";
print "<LI>";
print "<A HREF=\"/blocks-bin/sift.sh?$blks\" TARGET=\"_blank\">SIFT</A>";
print " to predict amino acid substitutions in blocks";
print "     [<A HREF=\"/sift/SIFT_on_blocks.html\">About SIFT</A>]";
print "<P>\n";

print "<P>";
print "<LI>";
print "<A HREF=\"/blocks-bin/blalign.sh?$blks\">Re-format</A>";
print " blocks as a multiple alignment ";
print "<P>\n";

#
print "</UL>\n";
print "</MENU>\n";
print "<P><HR>\n";

#------------------------------------------------------------------------
# Header
print "<PRE>";
open (IN, "<$head");
foreach $in_rec (<IN>)  { print "$in_rec"; }
close (IN);

print "</PRE><A HREF=\"#top\">[Return to top]</A><P>\n";

#-----------------------------------------------------------------------
# Output Blocks
foreach $in_rec (@blocks)
{
  if ($in_rec =~ m/^ID   /)
  {  $id_rec = $in_rec;   }
  elsif ($in_rec =~ m/^AC   /)
  {
     $ctemp = $in_rec;
     ($ac) = $ctemp =~ m/^AC   (\S+)/;
     $ac =~ s/;//g;
     print "</PRE><H3><A NAME=\"", $ac,"\">Block ", $ac, "</A>";
     print "</H3><P><PRE>\n";
     print "$id_rec";
     print "$in_rec";
  }
  elsif ($in_rec =~ m/^DE   / || $in_rec =~ m/^BL   /)
  {  print "$in_rec";  }
  else
  {
    $ctemp = $in_rec;
    $ctemp =~ s/\s+//g;
    if ($ctemp eq "") { print "$in_rec"; }	#  blank line
    elsif ($in_rec =~ /\s*\S+\ +\(/) 		#  sequence?
    {
      ($in_rec =~ /^\s*(\S+)/) && ($sid = $1) if $in_rec =~ /\s*\S+\ +\(/;
      #	Assume we have a sequence name == $sid; could be "id|ac", 
      #	if so, use just	the ac part in the link
      @bar = split(/\|/, $sid);
      if (@bar[1] ne "") {$sid = @bar[1]; }
      $in_rec =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$sid">\2</A>\3| if  $in_rec =~ /\s*\S+\ +\(/;
      print $in_rec;
    }  # end of sequence line
    else { print $in_rec; }	#  unknown line
  }
}  # end of blocks file
close(IN);
print "</PRE><A HREF=\"#top\">[Return to top]</A><P>\n";

#-------------------------------------------------------------
# Output the COBBLER sequence
#  Read it into memory
open (IN, "<$cob");
@cobseq = <IN>;
close(IN);

print "<A NAME=\"cobbler\"><H3>COBBLER sequence (region containing Blocks only)</H3></A>\n";

print "To do a BLAST search, copy the cobbler sequence";
print " below then click on a BLAST link<BR>\n";

#	Print sequence for blast search
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=0&PROGRAM=blastp&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
foreach $in_rec (@cobseq)
{ 
#  $in_rec =~ s/\s+//g;		# remove whitespace
   print "$in_rec"; 
}
print "\" TARGET=\"getblock\">Blast Search</A>]</B>\n";

#	Print sequence again for gap-blast search
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=1&PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
foreach $in_rec (@cobseq)
{ 
#  $in_rec =~ s/\s+//g;		# remove whitespace
   print "$in_rec"; 
}
print "\" TARGET=\"getblock\">Gap-Blast Search</A>]</B>\n";

#	Print sequence again for psi-blast search
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?AUTO_FORMAT=Semiauto&PROGRAM=blastp&RUN_PSIBLAST=on&CDD_SEARCH=on&COMPOSITION_BASED_STATISTICS=on&DESCRIPTIONS=250&ALIGNMENTS=100&ALIGNMENT_VIEW=Pairwise&DATABASE=nr&END_OF_HTTPGET=YES&QUERY=";
foreach $in_rec (@cobseq)
{ 
#  $in_rec =~ s/\s+//g;		# remove whitespace
   print "$in_rec"; 
}
print "\" TARGET=\"getblock\">PSI-Blast Search</A>]</B>\n";

#	Print sequence to the screen
print "<PRE>";
foreach $in_rec (@cobseq)
{ print "$in_rec"; }
print "</PRE><A HREF=\"#top\">[Return to top]</A><P>\n";
print "</HTML>";

#-------------------------------------------------------------------------
exit (0);

#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment variable
# returns: an associative array of name/value pairs.  The name is the key.
sub parse_query {
  local($query_string) = @_;
  local(%ans, @q, $pair);
#print $query_string;
  # break up into individual name/value lines
  @q = split(/&/, $query_string);

  foreach $pair (@q) {
    # break the name/value pairs up
    # use split rather than regular expressions because the value may have
    #  newlines in it
    split(/=/, $pair, 2);

    # change '+' to ' '
    $_[1] =~ s/\+/ /g;

    # change the escaped characters (has to be after the split on '&' and '=')
    $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

    $ans{$_[0]} = $_[1];
  }

  return %ans;
}

# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}
