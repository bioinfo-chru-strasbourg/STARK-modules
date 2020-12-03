#!/usr/bin/perl
#	htmlize-getblock
#  getblock output is assumed to consist of:
#	* Introduction from blksort.stp or prints.stp
#	* The blocks, followed by "%d blocks processed"
#	* The cobbler.pros sequence, followed by "end cobbler"
#	* The block map lines, followed by "end map"
#	* The tree lines, followed by "end tree"
# 	* The proweb.dat entries, followed by "end proweb"
# 	* The blinks.dat entries, followed by "end blinks"
# 	* The pdb.dat entries, followed by "end pdb"
#	* prosite entry (prosite.dat, prosite.doc)
#----------------------------------------------------------------
#  2/25/96  Modified for getblock with COBBLER sequence  JGH
#           Got rid of while (1) statements
#  5/29/96  Modified for getblock with PROWEB data.  JGH
#  7/11/96  Added limit on cputime, fixed error in while loop at
#           line 286 which caused looping.  JGH
#  7/12/96  Put eof checks in all while statements JGH
#  8/ 7/96  Revised internal prosite links.  JGH
#  9/ 3/96  Added tree link.  JGH
#  9/11/96  Added frames. JGH
# 11/ 6/96  Added MAST link.  JGH
#  3/ 6/97  Moved sequence links from SwissProt to NCBI to generalize
#	    for Prints sequence IDs. Check for PR or BL block_ACs.
#		logo.sh & LAMA_search still need work ...
#  3/ 9/97  logo.sh just once
#  7/31/97  Added maps (map.sh)
#  8/11/97  Moved logo, lama & mast links to the top
#  8/13/97  Removed frames, use target window for documentation instead.
#	    Use target window for Blast search.
#  9/ 4/97  Added about_trees, about_maps
#  9/ 8/97  Added "Mime" type for trees, newick.sh
#  9/18/97  Added link to Psi-Blast
# 10/ 1/97  Changed Blastp & Psi-Blastp for DESCRIPTIONS=500
# 10/ 9/97  Added link to Emotif
# 10/12/97  Renamed logo.sh logo.ps.sh; moved trees to top
#  5/ 8/98  Added pdb link. 
#  7/ 7/98  Added ProDom.
# 10/ 8/98  Added Domo.
# 12/ 1/98  Pre-made gif logos
# 12/ 6/98  Revised tree display shells (drawgram.csh)
# 12/ 8/98  logo.csh replaces logo.*.csh|sh
# 12/28/98  Don't edit ">", "<", etc. Causes problems with some browsers
#		when <PRE>
#  1/ 7/99  Expect blk2pdb.dat for other than BL; format changes
#  3/22/99  Use sprot links for all but prints
#  3/27/99  Add internal link for group number at top
#  6/11/99  New sequence names: id|ac
#  7/23/99  Changed pdbfile to call pdbfile.pl
#  8/21/99  Prints links to swiss now (23.1); may not be a tree
# 11/ 6/99  Blast links are broken ...
#  4/29/00  Blocks 12 (IPB...)
#  2/ 2/01 Added link to ProWeb tree display
# 11/29/01 Fix PSI-BLAST link
#  1/ 6/02 Rearranged tree options
#  1/18/02 Rearranged output
#  2/18/03 pdbfile changes
#  3/ 8/03 Replace pdbfile.pl with pdbmast.pl
#  4/13/03 Link to proweb graphical block maps
#  8/25/03 Link to blalign
#  3/26/04 Run LAMA & SIFT in a blank window
#  3/14/06 Removed Metafam link (no longer works)
#  4/ 8/07 Get pdb info from pdbs/ instead of pdbmast.dat
#----------------------------------------------------------------
#
#	Prevent this script from running forever - 30 seconds
system("limit cputime 30");

if (@ARGV == 0) {
  $infile = STDIN;
}
elsif (@ARGV > 1) {
  print "Usage:\nhtmlize-getblock.pl [getblock output]\n";
}
else {
  unless (open($infile, @ARGV[0])) {
    print STDERR "Can't open @ARGV[0]: $!\n";
    return;
  }
}

@input = <$infile>;

#     If a <FILEHANDLE> is used in a context that is  looking  for
#     an  array,  an  array  consisting  of all the input lines is
#     returned, one line per array element.  It's easy to  make  a
#     LARGE data space this way, so use with care.


#----------------------------------------------------------------
#	Read entire input file, saving some lines in arrays

$past_blocks = 0;
$past_intro = 0;
$past_cobb = 0;
$past_map = 0;
$past_tree = 0;
$past_web = 0;
$nweb = 0;
$past_lnk = 0;
$nlnk = 0;
$past_pdb = 0;
$npdb = 0;
$ntree = 0;

#	Write the blocks for blk2pssm to use later for MAST
open(FBLK, ">../tmp/$$.blocks");
open(FMAP, ">../tmp/$$.blocks.mapfile");
open(FTREE, ">../tmp/$$.treefile");
open(FPDB, ">../tmp/$$.pdbfile");
system("chmod -f 660 ../tmp/$$.*");
#
for ($j = 0; $j < @input; $j++)
{
  # beyond the blocks yet?
  if (@input[$j] =~ /^\d+ blocks processed/ ) 
  {   $past_blocks = 1; close(FBLK);   }

  # beyond the cobbler sequence?
  if (@input[$j] =~ /end cobbler/ ) 
  {   $past_cobb = 1;   }

  # beyond the block map?
  if (@input[$j] =~ /end map/ ) 
  {   $past_map = 1; close(FMAP);  }

  # beyond the tree?
  if (@input[$j] =~ /end tree/ ) 
  {   $past_tree = 1; close(FTREE); }

  # beyond proweb yet?
  $past_web = 1 if @input[$j] =~ /end proweb/ ;

  # beyond blinks yet?
  $past_lnk = 1 if @input[$j] =~ /end blinks/ ;

  # beyond the pdb?
  if (@input[$j] =~ /end pdb/ ) 
  {   $past_pdb = 1; close(FPDB); }

  # read the line numbers that the blocks start at
  if (@input[$j] =~ /^ID   / && !$past_blocks)
  {  @block_lines = (@block_lines, $j); 
     @block_ID = (@block_ID, (@input[$j] =~ /^ID\s+(.*);/)) ;
     $past_intro = 1; 
  }

  # read the block numbers for use in the TOC
  @block_AC = (@block_AC, (@input[$j] =~ /^AC\s+(.*);/))  
    if @input[$j] =~ /^AC   / && !$past_blocks ;

  # 	Description
  @block_DE = (@block_DE, (@input[$j] =~ /^DE\s+(.*)$/)) 
    if @input[$j] =~ /^DE   / && !$past_blocks ;

  #	Prepare the blocks
  if ($past_intro && !$past_blocks)
  {  print FBLK @input[$j];  }

  #	Prepare the block map file, don't write "end cobbler" line
  if ($past_cobb && !$past_map && !(@input[$j] =~ /end cobbler/))
  {  print FMAP @input[$j];  }

  #	Prepare the tree file, don't write "end map" or ">" title line
  if ($past_map && !$past_tree && !(@input[$j] =~ /end map/) &&
      !(@input[$j] =~ /^>/) )
  {  print FTREE @input[$j]; $ntree +=1;  }

  #	Keep track of proweb lines
  if ($past_tree && ! $past_web && !(@input[$j] =~ /end tree/))
  { $nweb +=1; }
  #	Keep track of prodom lines
  if ($past_web && ! $past_lnk && !(@input[$j] =~ /end proweb/))
  { $nlnk +=1; }
  #	Keep track of pdb lines
  if ($past_lnk && ! $past_pdb && !(@input[$j] =~ /end blinks/))
  { $npdb +=1; }

# Double subscript blink? Save line numbers instead of array?
# if ($past_proweb && !past_lnk && !(@input[$j] =~ /end blinks/) )
#    @blink_AC = (@blink_AC, (@input[$j] =~ /^AC\s+(.*);/))  
#	blink[0] = AC, blink[1] = link, blink[2] = description
#         @blink = split(/!/, @input[$j]);

  #	Prepare the pdb file, don't write "end blinks" or ">" title line
  if ($past_lnk && !$past_pdb && !(@input[$j] =~ /end blinks/) )
  {
     #	Just pass the info to pdbmast.pl now
     print FPDB @input[$j];
  }
#>>>>>>> PDB stuff is redone later, look for pdbmast.dat

#>>>>>>> prosite stuff is currently missing; see getblock.sh
  # read the prosite ID for use in the title and TOC
  @prosite_ID = (@prosite_ID, (@input[$j] =~ /^ID\s+(.*);/)) 
    if @input[$j] =~ /^ID   / && $past_blocks && $past_web ;

  # read the prosite AC for use in the TOC
  @prosite_AC = (@prosite_AC, (@input[$j] =~ /^AC\s+(.*);/)) 
    if @input[$j] =~ /^AC   / && $past_blocks && $past_web ;

  # read the prosite DE for use in the TOC and title
  @prosite_DE = (@prosite_DE, (@input[$j] =~ /^DE\s+(.*)$/)) 
    if @input[$j] =~ /^DE   / && $past_blocks && $past_web ;

  # read the prosite DO for use in the TOC and title
  @prosite_DO = (@prosite_DO, (@input[$j] =~ /^DO\s+(.*);/)) 
    if @input[$j] =~ /^DO   / && $past_blocks && $past_web ;

  # read in the prosite sequence numbers and names
  # this reads a string like this (minus the '#'): 
#DR   P02820, OSTC_BOVIN, T; P02822, OSTC_CHICK, T; P15504, OSTC_DRONO, T;
#  $names{OSTC_BOVIN} = P02820 and $prosite{P02820} = OSTC_BOVIN

  if (@input[$j] =~ /^DR   /) {
    @bar = split(/; */, @input[$j]);

    foreach $bar (@bar) {
      ($bar =~ /(\S+)\s*,\s+(\S+)\s*,/) 
       && ($names{$2} = $1)
       && ($prosite{$1} = $2);
    }
  }
}     #  end of input[j]
#=========================================================================
#   Now entire input file has been read & stored
#   Check to see if any blocks were found
if (@block_lines == 0) {
  print "</PRE>

<TITLE>Getblock Error</TITLE>
<H1>Getblock Error</H1>
Invalid block number or nonexistent block.<P>
 
The correct format for the block number should be IPB??????.<P>
 
For example:  IPB001525, where '0' is the number zero.<P>
  
";
  exit;
}

#--------------------------------------------------------------------------

#	Get the number part of the block name
$interpro_flag = 1;
$blocks_flag = 0;
$prints_flag = 0;
$prodom_flag = 0;
$domo_flag = 0;
$pfam_flag = 0;

(@block_AC[0] =~ /^IPB(\d+)/) && ($block_group_num = $1);
if (! $block_group_num)
{
#		Try for old blocks
   (@block_AC[0] =~ /^BL(\d+)/) && ($block_group_num = $1);
   if ($block_group_num) 
   { $blocks_flag = 1; $interpro_flag = 0; }
#		Try for Prints
   else
   {
      (@block_AC[0] =~ /^PR(\d+)/) && ($block_group_num = $1);
      if ($block_group_num) 
      { $prints_flag = 1; $interpro_flag = 0; }
#		Try for Prodom
   else
   {
      (@block_AC[0] =~ /^PD(\d+)/) && ($block_group_num = $1);
      if ($block_group_num) { $prodom_flag = 1;  $interpro_flag = 0;}
      else
      {
         (@block_AC[0] =~ /^BP(\d+)/) && ($block_group_num = $1);
         if ($block_group_num) { $prodom_flag = 1;  $interpro_flag = 0;}
#		Try for Domo
         else
         {
            (@block_AC[0] =~ /^DM(\d+)/) && ($block_group_num = $1);
            if ($block_group_num) { $domo_flag = 1;  $interpro_flag = 0;}
#		Try for Pfam
            else
            {
               (@block_AC[0] =~ /^PF(\d+)/) && ($block_group_num = $1);
               if ($block_group_num) { $pfam_flag = 1;  $interpro_flag = 0;}
            }
         }
      }
   }
   }
}     #  end if not interpro

if (! $interpro_flag)
{   $group = substr(@block_AC[0], 0, 7); }
else
{   
      $group = substr(@block_AC[0], 0, 9); 
      $ipr =  "IPR".$block_group_num;
}
   
#  Output title 
print "</PRE>\n";
print "<A NAME=\"",$group,"\"></A>\n";
print "<TITLE>Blocks for $group</TITLE>\n";
print "<A NAME=top><H1>$group: @block_ID[0]</H1></A>\n";
print "@block_DE[0]<P>\n";

print "<MENU>\n";
print "<UL>\n";
#print "<LI> <A HREF=\"#bintro\">Introduction</A>\n";
print "<LI> <A HREF=\"/blocks/help/blocks_format.html\">Introduction</A>\n";
$nblock = 0; $last_ac = "";
foreach $blk (@block_AC) {
  print "<LI> <A HREF=\"#", $blk, "\">Block number ", $blk, "</A>\n";
  $nblock++; $last_ac = $blk;
}

print "<P>";
if ( $prints_flag )
{
   print "<LI>";
   print "PRINTS Entry <A HREF=\"http://www.bioinf.man.ac.uk/cgi-bin/dbbrowser/PRINTS/DoPRINTS.pl?cmd_a=Display&qua_a=/Brief&fun_a=Text&qst_a=$group\"TARGET=\"getblock\"> $group</A>";
   print " (source of blocks)\n";
}
if ( $blocks_flag )
{
#>>>>>>  prosite_DO isn't picked up unless getblock reads prosite.dat ...
   $pdat = "PS".$block_group_num;
   $pdoc = "PS".$block_group_num;
   print "<LI>";
   print "PROSITE entry <A HREF=\"http://www.expasy.ch/cgi-bin/get-prosite-entry?$pdat\"TARGET=\"getblock\"> $pdat</A>" ;
   print " (source of sequences used to make blocks)\n";
}
if ( $interpro_flag )
{
   print "<LI>";
   print "InterPro entry <A HREF=\"http://www.ebi.ac.uk/interpro/IEntry?ac=$ipr\"TARGET=\"getblock\"> $ipr</A>" ;
   print " (source of sequences used to make blocks)\n";
}

if ($prints_flag) { $maptype = "PR"; }
else { $maptype = "IPB"; }
print "<P>\n";
print "<LI>Block Maps.";
print "[<A HREF=\"http://www.proweb.org/proweb-bin/blockmap.cgi?name=$group&type=F
AM&condensed=YES&dbtype=$maptype\"TARGET=\"getblock\">Graphical Map</A>]";
print "     [<A HREF=\"/blocks-bin/map.csh?../tmp/$$.blocks.mapfile\"TARGET=\"getblock\">Text Map</
A>]";
print "     [<A HREF=\"/blocks/tmp/$$.blocks.mapfile\"TARGET=\"getblock\">Map Positions</A>]";
print "     [<A HREF=\"/blocks/help/about_maps.html\"TARGET=\"getblock\">About Maps</A>]";

print "<P>";
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
   print "<A HREF=\"/blocks-bin/logo.csh?../tmp/$$.blocks+gif\">[GIF]</A> ";
   print "<A HREF=\"/blocks-bin/logo.csh?../tmp/$$.blocks+pdf\">[PDF]</A> ";
}
print "<A HREF=\"/blocks-bin/logo.csh?../tmp/$$.blocks+ps\">[Postscript]</A> ";

#	Some groups with < 4 seqs don't have trees
if ( -s "../tmp/$$.treefile" )	# $ntree > 0 ???
{
print "<P>\n";
print "<LI>Tree from blocks alignment. ";
print "[<A HREF=\"/blocks/help/about_trees.html\">About Trees</A>]  ";
print "[<A HREF=\"http://www.proweb.org/treeinfo.html\">About ProWeb TreeViewer</A>]  ";
print "<BR>[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+dat\">Data</A>]  ";
print " [<A HREF=\"http://www.proweb.org/proweb-bin/trees.cgi?family=$group\" TARGET=\"proweb\">ProWeb TreeViewer</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+xbm\">XBitmap</A>]  ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+gif\">GIF</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+pdf\">PDF</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+ps\">Postscript</A>]  ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+new\">Newick</A>] ";
}

#>>>> Re-do the pdb stuff from pdbmast.dat
#if ($prints_flag > 0) { $pdbfile = "../data-prints/pdbmast.dat"; }
#else                  { $pdbfile = "../data-blocks/pdbmast.dat"; }
#$pdbtemp = "../tmp/$$.pdbmast";
#system("grep \"^$group\" $pdbfile > $pdbtemp");
if ($prints_flag > 0) { $pdbtemp = "../data-prints/pdbs/$group.pdb"; }
else                  { $pdbtemp = "../data-blocks/pdbs/$group.pdb"; }
#if ($npdb > 0)
if (!(-z $pdbtemp))
{
   print "<P><LI><A HREF=\"/blocks-bin/pdbmast.pl?$pdbtemp+$last_ac\">Structures</A>\n";
}

print "<P>";
print "<LI>Search blocks vs other databases:<BR><UL>\n";
print "<LI><A HREF=\"#cobbler\">COBBLER sequence</A>\n";
print " and BLAST searches";
print "     [<A HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]";

print "<LI>";
print "<A HREF=\"/blocks-bin/mast.sh?../tmp/$$.blocks+$group\" TARGET=\"getblock\">MAST Search</A>";
print " of all blocks vs a sequence database ";
print "     [<A HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">About MAST</A>]";

print "<LI>";
print "<A HREF=\"/blocks-bin/LAMA_search.sh?../tmp/$$.blocks\" TARGET=\"_blank\">LAMA search</A>";
print " of all blocks vs a blocks database  ";
print "     [<A HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]";
print "</UL>";

print "<P>";
print "<LI>";
print "<A HREF=\"/blocks-bin/codehop.sh?../tmp/$$.blocks\">CODEHOP</A>";
print " to design PCR primers from blocks";
print "     [<A HREF=\"/blocks/help/CODEHOP/tips.html\">About CODEHOP</A>]";

print "<P>";
print "<LI>";
print "<A HREF=\"/blocks-bin/sift.sh?../tmp/$$.blocks\" TARGET=\"_blank\">SIFT</A>";
print " to predict amino acid substitutions in blocks";
print "     [<A HREF=\"/sift/SIFT_on_blocks.html\">About SIFT</A>]";
print "<P>\n";

print "<P>";
print "<LI>";
print "<A HREF=\"/blocks-bin/blalign.sh?../tmp/$$.blocks\">Re-format</A>";
print " blocks as a multiple alignment ";
print "<P>\n";

#
if ($nweb > 0 || $nlnk >0)
{print "<LI> <A HREF=\"#proweb\">Additional Links</A>\n"};
print "</UL>\n";
print "</MENU>\n";
print "<P><HR>\n";

print "<PRE>\n";

#------------------------------------------------------------------------
# Output Intro (header & body) (to first block)
#print "<A NAME=\"bintro\"><H3>Introduction</H3></A><P>\n";

for ($i=0; $i<@block_lines[0]; ++$i) {
  print @input[$i];
}

print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";


#-----------------------------------------------------------------------
# Output Blocks headers and bodies with links inserted
#>>> attach blinks to appropriate block somehow

for ($j=0; $j<@block_lines ; ++$j) {
  print "<H3><A NAME=\"", @block_AC[$j],"\">Block ", @block_AC[$j], "</A>";
  print "</H3>";
  print "<P>\n";
  $i = @block_lines[$j];
#	could just use $id here instead of @names{$id} ...
  while ($i < @input) {
    if (@input[$i] =~ /\s*\S+\ +\(/) 
    {
      (@input[$i] =~ /^\s*(\S+)/) && ($id = $1) if @input[$i] =~ /\s*\S+\ +\(/;
#	Assume we have a sequence name == $id; could be "id|ac", if so, use just
#	the ac part in the link
      @bar = split(/\|/, $id);
      if (@bar[1] ne "") {$id = @bar[1]; }
      @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$id">\2</A>\3| if  @input[$i] =~ /\s*\S+\ +\(/;
#     if (! $prints_flag)		# link to Swiss-Prot
#     {
#        @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$id">\2</A>\3| if  @input[$i] =~ /\s*\S+\ +\(/;
#     }
#     else			# link to Entrez
#     {
#        @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www3.ncbi.nlm.nih.gov:80/htbin-post/Entrez/query?uid=$id&form=6&db=p&Dopt=g">\2</A>\3| if  @input[$i] =~ /\s*\S+\ +\(/;
#     }
      print @input[$i];
    }
    else {
      print @input[$i];
    }
    last if @input[$i] =~ /^\/\//;
    ++$i;
  }
  print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";
}

# Output lines following blocks ?

++$i;
while ($i < @input) {
  last if @input[$i] =~ /^\d* blocks processed/;
  print @input[$i];
  ++$i;
}
#		Don't print the "blocks processed" line
++$i;
 
#-------------------------------------------------------------
# Output the COBBLER sequence
print "<A NAME=\"cobbler\"><H3>COBBLER sequence (region containing Blocks only)</H3></A>";

print "</PRE>To do a BLAST search, copy the cobbler sequence";
print " below then click on a BLAST link<PRE><BR>\n";

#	Print sequence for blast search
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/cgi-bin/BLAST/nph-blast?PROGRAM=blastp&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=0&PROGRAM=blastp&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
#   Save the position, put the sequence in the blast link
$j = $i;
$i = $j + 3;				#  First three lines are headings
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
print "\" TARGET=\"getblock\">Blast Search</A>]</B>";
#print " (copy the COBBLER sequence first)\n";

#	Print sequence again for gap-blast search
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=1&PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
$i = $j + 3;
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
print "\" TARGET=\"getblock\">Gap-Blast Search</A>]</B>";
#print " (copy the COBBLER sequence first)\n";

#	Print sequence again for psi-blast search
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/psiblast.cgi?PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?AUTO_FORMAT=Semiauto&PROGRAM=blastp&RUN_PSIBLAST=on&CDD_SEARCH=on&COMPOSITION_BASED_STATISTICS=on&DESCRIPTIONS=250&ALIGNMENTS=100&ALIGNMENT_VIEW=Pairwise&DATABASE=nr&END_OF_HTTPGET=YES&QUERY=";
$i = $j + 3;
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
print "\" TARGET=\"getblock\">PSI-Blast Search</A>]</B>";
print "\n";

#	Print sequence to the screen
$i = $j;
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
#	Don't print the "end cobbler" line
++$i;
print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";

#-------------------------------------------------------------------------
#	Skip to the proweb entries

while ($i < @input) {
  last if @input[$i] =~ /end tree/;
  ++$i;
}
#		Skip the "end tree" line
++$i;

#-------------------------------------------------------------
if ($nweb > 0 || $nlnk > 0)
{
print "<A NAME=\"proweb\"><H3>Additional Links (separate browser window)</H3></A><P>\n";
}

# Output the PROWEB data
while ($i < @input) {
  last if @input[$i] =~ /end proweb/;
  if (!(@input[$i] =~ /^ID   /) && !(@input[$i] =~ /^AC   /) &&
      !(@input[$i] =~ /^\/\//) ) {
          @bar = split(/!/, @input[$i]);
          print "<A HREF=\"", @bar[0], "\"TARGET=\"getblock\">",  @bar[1], "</A>";
          print "<BR>";
  }
  ++$i;
}
#	Don't print the "end proweb" line
++$i;

#	Metafam link
#if (! $blocks_flag)
#{
#   print "<A HREF=\"http://metafam.ahc.umn.edu/superSet.jsp?SearchType=ProteinFamily&SearchString=$group&ShowApplet=TRUE&ShowKeys=TRUE\"TARGET=\"getblock\">MetaFam $group</A>";
#}


#	Output a blank line before starting more links
print "<P>\n";

# Output the ProDom links in blinks.dat - note @bar[] subscripts!
while ($i < @input) {
  last if @input[$i] =~ /end blinks/;
  if (!(@input[$i] =~ /^ID   /) && !(@input[$i] =~ /^AC   /) &&
      !(@input[$i] =~ /^\/\//) ) {
          @bar = split(/!/, @input[$i]);
          print "@bar[0]: <A HREF=\"", @bar[1], "\"TARGET=\"getblock\">",  @bar[2], "</A>";
          print "<BR>";
  }
  ++$i;
}
#	Don't print the "end blinks" line
++$i;

#--------------------------------------------------------------

exit;
