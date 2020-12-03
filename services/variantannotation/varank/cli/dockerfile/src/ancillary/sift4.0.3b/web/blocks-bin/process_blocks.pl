#!/usr/bin/perl
#		process_blocks.pl (executed by www/process_blocks.html)
#
#	NOTE: Would have to input the sequences that went into
#	making the blocks in order to make a cobbler sequence.
#	Also have to input sequences in order to add sequences to blocks.
#--------------------------------------------------------------------------
#  7/19/02 Fix 3dblocks link; add tree viewer link
#  8/14/02 Put logos, etc. in separate window
#  8/29/02 Use hidden fields for 3dblocks because
#          sometimes the input is too big for Apache on www.proweb.org
#  1/ 3/03 Use hidden fields for treeviewer
#  1/ 8/03 Pass treeviewer blocks & seqs
#  4/14/03 Proweb block map
#  8/27/03 Added DIY search
#  8/30/03 Added minwidth & maxwidth
#  3/26/04 Run LAMA & SIFT in a blank window
#  2/18/07 Use blockmap.pl here instead of proweb
#  4/ 9/07 Use treeviewer.pl here instead of proweb
#  7/10/07 Limit input size to $maxchar characters
# 12/11/07 Capture HTTP_USER_AGENT just for info
#--------------------------------------------------------------------------
#
$ID = $$;
$maxchar = 100000;
#$maxchar = 10000;		#>>>testing

$bin = ".";
$tmp = "../tmp";
#		These files are removed after use
$blin = "$tmp/$ID.blin";
$blks = "$tmp/$ID.blks";
$cblks = "$tmp/$ID.cblks";
$seqs = "$tmp/$ID.seqs";
$err = "$tmp/$ID.err";
$out = "$tmp/$ID.out";
#		These files are left
$wblks = "$tmp/$ID.wblks";
$tree = "$tmp/$ID.treefile";
#	htmlize-map.pl expects this format!
$map = "$wblks.map";
$pssm = "$tmp/$ID.pssm";
$mast = "$tmp/$ID.mast";

#	Be sure files and directories are created with rw-rw---- permissions
#system("umask 006");

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>Block Information</TITLE>\n";

open(TMP, ">>../tmp/agents");
print TMP "agent=$ENV{\"HTTP_USER_AGENT\"} referer=$ENV{\"HTTP_REFERER\"}\n";
close(TMP);

# check that this is a POST submission

if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
    print "This script should be referenced with a METHOD of POST. \n";
    exit;
}

#if ($ENV{"CONTENT_TYPE"} ne "application/x-www-form-urlencoded") {
#    print "This script can only be used to decode form results. \n";
#    exit;
#}


# get the QUERY_STRING environment variable and the name/value pairs
read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

# check that there are blocks
if ($names{seqfile} eq "" && $names{sequences} eq "") {
    print "<H1>Error</H1>Please enter some blocks.<P>\n";
    exit;
}

#	Get rid of everything but numbers
#print "minwidth=$names{minwidth} maxwidth=$names{maxwidth}\n";
$names{minwidth} =~ s/\D//g;
$names{maxwidth} =~ s/\D//g;
$minwidth = 10;
if ($names{minwidth} > 0) { $minwidth = $names{minwidth}; }
$maxwidth = 55;
if ($names{maxwidth} >= $minwidth) { $maxwidth = $names{maxwidth}; }


#  create a file of the input multiple alignment
#  Limit the size
if (length($names{seqfile}) > $maxchar)
{ 
   $x = substr($names{seqfile}, 0, $maxchar);
   $names{seqfile} = $x;
   print "<B>Warning:</B> Your input has been truncated to $maxchar characters.<P>\n";
} 
if (length($names{sequences}) > $maxchar)
{ 
   $x = substr($names{sequences}, 0, $maxchar);
   $names{sequences} = $x;
   print "<B>Warning:</B> Your input has been truncated to $maxchar characters.<P>\n";
} 
#  Notice input from both sources will be concatenated
open(BLK, ">$blin");
if ($names{seqfile} ne "")
{  print BLK $names{seqfile}; }
if ($names{sequences} ne "")
{  print BLK $names{sequences};   }
print BLK "\n";
close(BLK);

#  Be sure blocks are in the correct format: puts output in
#	$tmp/$ID.blks = $blks and $tmp/$ID.seqs 
#>>>> empty .seqs file when blocks input?
#  Can handle blocks, fasta and clustal format multiple alignments
#  Adds sequence weights as well, but in floating point format
#print "minwidth=$minwidth maxwidth=$maxwidth\n";
system("$bin/blimps-bin/mablock $blin $tmp/$ID B $minwidth $maxwidth > $err 2>$err");

#	Check that $blks has non-zero size
if (-s "$blks") 
{
   #  calibrate the blocks for diy searching
   system("$bin/calibrate.csh $blks $cblks > /dev/null");
   #  create a file of PB-weighted blocks
   # NOTE: Should give them the option of using existing weights (if any)
   #       or other weighting schemes
   system("$bin/blimps-bin/blweight $cblks $wblks P M > /dev/null");

   #  create the PSSM files == $pssm & $mast
   system("$bin/blimps-bin/blk2pssm $wblks $pssm B > /dev/null");
   system("$bin/blimps-bin/blk2pssm $wblks $mast M > /dev/null");

   #  create the treefile == $tree
   system("$bin/maketree.csh $wblks > /dev/null");

   #  create the map file == $map
   system("$bin/blimps-bin/makeblockmap $wblks $map > /dev/null 2>&1");

   #  remove unnecessary files
#  system("rm $blin $blks $tmp/$ID.seqs");
}

#=========================================================================
#	Assuming here that all the other files got made if $wblks has
#	non-zero size
if (-s "$wblks") 
{
  print "<H1>Block Information</H1>\n";
# print "<STRONG>";
  print "<PRE>";
  open(OUT, "cat $out |");
  while ($_ = <OUT>) { print; }
  close(OUT);	
  print "</PRE>";

  print "<P><I>Sequence-Weighted Blocks:</I><BR>";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/catfile.sh?$wblks\">Blocks Format</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/catfile.sh?$pssm\">Blimps PSSM</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/catfile.sh?$mast\">MAST PSSM</A>] ";
  print "[<A  TARGET=\"blocks\"HREF=\"/blocks/help/PSSM_def.html\">About PSSMs</A>] ";
  print "<BR>";
  # The blockmap name must match what's in the mapfile
#<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/blockmap.cgi\">
  open(GREP, "grep \"^>\" $map |");
  $map_rec = <GREP>;
  ($mapname) = $map_rec =~ m/^>(\S+)/;
  print "<P><I>Map:     </I><BR>";
  print "
<FORM METHOD=POST ACTION=\"/blocks-bin/blockmap.pl\">
<INPUT TYPE=hidden NAME=\"name\" VALUE=\"$mapname\">
<INPUT TYPE=hidden NAME=\"type\" VALUE=\"FAM\">
<INPUT TYPE=hidden NAME=\"condensed\" VALUE=\"YES\">
<INPUT TYPE=hidden NAME=\"dbtype\" VALUE=\"USER\">
<INPUT TYPE=hidden NAME=\"map\" VALUE=\"";
  open(MAP, "<$map");
  while ($maprec = <MAP>)
  { print $maprec; }
  close(MAP);
  print "\"><INPUT TYPE=submit VALUE=\"Graphical Map\"></FORM>  \n";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/map.csh?$map\">Text Map</A>]  \n";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/tmp/$wblks.map\">Map Positions</A>]  \n";
  print "[<A HREF=\"/blocks/help/about_maps.html\"TARGET=\"blocks\">About Maps</A>]";

  print "<P><I>Logos:     </I><BR>";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/logo.csh?$wblks+ps\">Postscript</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/logo.csh?$wblks+pdf\">PDF</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/logo.csh?$wblks+gif\">GIF</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/about_logos.html\">About logos</A>]\n ";

  print "<P><I>Tree:     </I><BR>";
# print "<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/trees.cgi\"
  print "<FORM METHOD=POST ACTION=\"/blocks-bin/treeviewer.pl\"
 TARGET=\"blocks\">
<INPUT TYPE=hidden NAME=\"format\" VALUE=\"png\">
<INPUT TYPE=hidden NAME=\"length\" VALUE=\"yes\">
<INPUT TYPE=hidden NAME=\"tree\" VALUE=\"";
  open(TRE, "<$tree");
  while ($trerec = <TRE>)
  { print $trerec; }
  close(TRE);
  print "\"><INPUT TYPE=hidden NAME=\"blocks\" VALUE=\"";
  open(BLK, "<$wblks");
  while ($blkrec = <BLK>)
  { print $blkrec; }
  close(BLK);
  print "\"><INPUT TYPE=hidden NAME=\"seqs\" VALUE=\"";
  open(SEQ, "<$seqs");
  while ($seqrec = <SEQ>)
  { print $seqrec; }
  close(SEQ);
  print "\"><INPUT TYPE=submit VALUE=\"ProWeb TreeViewer\"></FORM>\n";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/maketree.csh?$wblks+xbm\">XBitmap</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/maketree.csh?$wblks+ps\">Postscript</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/maketree.csh?$wblks+pdf\">PDF</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/maketree.csh?$wblks+gif\">GIF</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/maketree.csh?$wblks+new\">NEW</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/maketree.csh?$wblks+dat\">Data</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/about_trees.html\">About trees</A>]\n";

#	Don't have a cobbler sequence to send to 3dblocks.cgi
# print "&SEQ=";
  print "<P><I>Structures (takes several minutes):</I><BR>";
  print "
<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/3dblocks.cgi\"
 TARGET=\"3dblocks\">
<INPUT TYPE=hidden NAME=\"METHOD\" VALUE=\"Mast\">
<INPUT TYPE=hidden NAME=\"CUTOFF\" VALUE=\"0.0001\">
<INPUT TYPE=hidden NAME=\"SAVE\" VALUE=\"NO\">
<INPUT TYPE=hidden NAME=\"blocks_data\" VALUE=\"";
  open(BLK, "<$wblks");
  while ($blkrec = <BLK>)
  { print $blkrec; }
  close(BLK);
  print "\"><INPUT TYPE=submit VALUE=\"3D Blocks\"></FORM>\n";
  print "  [<A TARGET=\"blocks\" HREF=\"http://www.proweb.org/3dblocks_intro.html\">About 3D Blocks</A>]</FORM>\n";

  print "<P><I>Multiple alignment search:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/LAMA_search.sh?$wblks\">LAMA</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]   ";
  print "<SPACER type=horizontal size=10>";
  $group = "UserBlocks";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/mast.sh?$wblks+$group\">MAST</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">About MAST</A>]\n";

  print "<P><I>COBBLER sequence and single sequence search:</I><BR>";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/cobbler.sh?$wblks\">COBBLER</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]\n";

  print "<P><I>Search these blocks:</I><BR>";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/diy.sh?$wblks\">DIY Search</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/about_diy.html\">About DIY</A>]\n";

  print "<P><I>Primers:</I><BR>";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/codehop.sh?$wblks\">CODEHOP</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/CODEHOP/CODEHOP_help.html\">About CODEHOP</A>]\n";

  print "<P><I>Substitutions:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/sift.sh?$wblks\">SIFT</A>] ";
  print "[<A TARGET=\"blocks\" HREF=\"/sift/SIFT_on_blocks.html\">About SIFT on Blocks</A>]\n";

}
else 
{ 
  print "<PRE>Problem processing input:<BR>";
  open(ERR, "cat $err |");
  while ($_ = <ERR>) { print; }
  close(ERR);	
  print "\n<P><\PRE>Please refer to "; 
  print "[<A TARGET=\"blocks\" HREF=\"/blocks/blocks_format.html\">acceptable formats</A>].";
}

print "<P>\n";

exit(0);

#=========================================================================
# parameter: a string that is the html QUERY_STRING environment variable
# returns: an associative array of name/value pairs.  The name is the key.
sub parse_query {
  local($query_string) = @_;
  local(%ans, @q, $pair);

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

#-------------------------------------------------------------------------
#  NOTE: Not really general, does special stuff for names{seqfile}
#	 and for names{sequences}
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:	CONTENT_TYPE
#		QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the key.

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#		Content-Disposition: form-data; name="seqfile"; filename="file.pros"
#		<blank line>
#		<content of file.pros>
#		<boundary>
#		Content-Disposition: form-data; name="sequences"
#		<blank line>
#		<data from input form>
#		<boundary>

sub general_parse {
  local($content_type, $query_string) = @_;
  local(%ans, @q, $pair, $loc, $boundary, $temp, $temp1);

#print "<PRE>$content_type\n\n";
#print "query_string:\n$query_string\n\n";

  if ($content_type eq "application/x-www-form-urlencoded")
  {
     # break up into individual name/value lines
     @q = split(/&/, $query_string);

     foreach $pair (@q) {
       # break the name/value pairs up
       # use split rather than regular expressions because the value may have
       #  newlines in it
       split(/=/, $pair, 2);

       # change '+' to ' '
       $_[1] =~ s/\+/ /g;
   
       # change the escaped characters (must be after the split on '&' and '=')
       $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

       $ans{$_[0]} = $_[1];
     }

  }
  else
  {
     $loc = index($content_type, "boundary=");
     if ($loc > 0)
     { 
        $temp = substr($content_type, $loc+9);
#	Why is this necessary? (boundary= sometimes doesn't match actual?)
#	$boundary = "--".$temp;
        $boundary = $temp;
#print "boundary=$boundary\n\n";
        # break up into individual name/value lines
        @q = split(/$boundary/, $query_string);
        foreach $pair (@q) {
          # break the name/value pairs up
#print "pair=$pair\n\n";
          $loc = index($pair, "name=");
	  $temp = substr($pair, $loc+5);
#	  $loc = index($temp, "\n\n");
 	  $loc = index($temp, "\n");
	  $temp1 = substr($temp, $loc+2);
#print "1 temp=$temp\n";
#		Get rid of stuff after the name
#	Need to look for a ; before the first \n ... <<< find another way!
	  $loc1 = index($temp, ";");
          $loc = index($temp, "\n");
          if ($loc1 > 0 && $loc1 < $loc) { $loc = $loc1; }
 	  if ($loc > 0) { $temp = substr($temp, 0, $loc); }
#print "2 loc=$loc temp=$temp\n";
#		Get rid of quotes around the name
          $temp =~ s/\"//g;
#print "3 temp=$temp\n";
#		Still has a trailing whitespace character ...
          $temp =~ s/\s//g;

#		Need to strip leading/ending whitespace off of $temp1,
#		but be careful not to strip off internal CRs in "seqfile:
#		and "sequences"
#print "4 temp=$temp\ntemp1=$temp1\n";
	  if ($temp ne "seqfile" && $temp ne "sequences")
 	  { $temp1 =~ s/\s//g; }

#		MAC file lines end in just \r, no \n;
#		DOS file lines end in \r\n; UNIX in \n.
#		mablock uses fgets() which requires \n.
#		Change \r\n to \n, then change \r to \n
 	  if ($temp eq "seqfile")
  	  { $temp1 =~ s/\r\n/\n/g;  $temp1 =~ s/\r/\n/g; }
#print "temp=$temp\ntemp1=$temp1\n</PRE>";
	  if ($temp ne "") { $ans{$temp} = $temp1; }
        }
     }
     else
     {  print "Cannot parse $content_type\n"; }
  }
  return %ans;
}

# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

