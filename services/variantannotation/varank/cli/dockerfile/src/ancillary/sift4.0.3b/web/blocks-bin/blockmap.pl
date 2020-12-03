#!/usr/bin/perl
#>>>wais doesn't work for DCM_ECOLI|P11876
#>>>seq lengths in prints maps are truncated to last block

# blockmap.html: name=,type=FAM|SEQ,condensed=YES|NO
# getblock.pl: name=IPB...|PR...,dbtype=IPB|PR,type=FAM,condensed=YES
# process_blocks.pl: name=,dbtype=USER,type=FAM,condensed=YES,map=
# bm_format.pl: name=,dbtype=USER,type=FAM,condensed=YES,map=

# Expected CGI input to display a block family
#  name      : BlockFamilyName|SeqID [,SeqID]* (not needed for user defined input)
#  condensed : YES, NO
#  type      : FAM, SEQ
#  dbtype    : IPB, PR, USER
#  map       : <filename>

#--------------------------------------------------------------------------
# 2/15/07 Moved from bateson to morgan, dropped mysql stuff,
#		dropped chromdb stuff, simplified considerably

#-----------------------------MODULES--------------------------------------
use strict;
use CGI::Carp;
use CGI qw/:standard/;

my $tmpdir = "../tmp";
my $waisd = "../index";

#database files
my $BLOCKS_DATABASE   = "../data-blocks/maps";
my $PRINTS_DATABASE   = "../data-prints/maps";

#my $htmlUrl                  = "http://blocks.fhcrc.org/blockmap";
my $htmlUrl                  = "/blocks/blockmap";
my $linePicture              = "$htmlUrl/graphics/line2.jpg";
my $blockPicture             = "$htmlUrl/graphics/blocks.jpg";
my $colorsFile           = "../www/blockmap/colors.txt"; 
my $databaseFile;


# database type
my $ALL            = 0;		# never used
my $BLOCK          = 1;
my $PRINT          = 2;
my $CHROM          = 3;		# never used
my $USER           = 4;

# input types
my $FAM            = 0;
my $SEQ            = 1;
my $HITS           = 2;		# never used

# graphical display variables
my $DISPLAYHEIGHT  = 7;
my $SCALINGFACTOR  = 1/2;

#is Condensed type
my $TRUE            = 1;
my $FALSE           = 0;

#--------------------MAIN-----------------------------------------
# declare variables
my $arraySize;
my $key;
my $cgi = new CGI;
my %input;
my %inputVarHash;
#	Defaults
#my @name;			why list?
my $name = "";
my $type = $FAM;
my $dbType = "";
my $isCondensed = 1;
my $mapFile = "";
my ($in_rec, @fams, $x, $done, $in_name, $in_len, $in_nblk, $ib);

# Hashes for input variables from form
%inputVarHash      = ('ALL',$ALL, 'CH',$CHROM, 'IPB', $BLOCK, 'PR', $PRINT, 'USER', $USER, 'FAM', $FAM, 'SEQ', $SEQ,'HITS', $HITS, 'YES', $TRUE, 'NO', $FALSE);

print header();


#read input from the forms
for $key ( $cgi->param() ) 
{
    $input{$key} = $cgi->param($key);
    
    if ($key eq "type")		# SEQ or FAM
    {
	$type = $inputVarHash{$input{$key}};
        if ($type ne $FAM && $type ne $SEQ)
        {
           print "<P><B>WARNING:</B> Ignoring unknown type $type, assuming FAM.<P>\n";
           $type = $FAM;
        }
    }
    # dbtype=blank from form, IPB or PR from getblock, USER from process_block
    elsif($key eq "dbtype")
    {
	$dbType = $inputVarHash{$input{$key}};
        if ($dbType ne "$BLOCK" && $dbType ne "$PRINT" && $dbType ne "$USER")
        {  
           print "<P><B>WARNING:</B> Ignoring unknown dbtype $dbType<P>\n";
           $dbType = ""; 
        }
    }
    elsif ($key eq "condensed")
    {
	   $isCondensed =$inputVarHash{$input{$key}};
    }
#   elsif ($key =~ /name*/)		???multiple names ever used???
    # name = protein seq id or IPB_id or PR_id
    elsif ($key eq "name")
    {
	#make sure it's not a blank entry
	if (($input{$key}=~/^\D.*/))
        {
	    $name=$input{$key};
            if ($type eq $FAM && $dbType ne $USER)
            {
               if ($name =~ m/^IPB\d\d\d\d\d\d/) { $dbType=$BLOCK; }
               elsif ($name =~ m/^PR\d\d\d\d\d/) { $dbType=$PRINT; }
            }
	}
    }
    elsif ($key eq "map")	# user-provided map
    {
	my $mapUser = $input{$key};

        #write contents of map to a tempfile
	if ( $mapUser  ne "")
        {
           $mapFile = "$tmpdir/$$.map";
	   open (OUT, ">$mapFile");
	   print OUT $mapUser;  
	   close (OUT);
	}
    }
}  # end of input fields

#-----------CHECK INPUT-----------------------------------
my $nerr = 0;
if ( $name eq "")
{
   print "<P><B>ERROR</B>: Missing sequence or family identifier.<P>\n";
   $nerr++;
} 

if ( $dbType eq $USER )
{
    if ($type ne $FAM)
    {
	print "<P><B>Warning:</B> Only family identifier is valid for user-provided map.<P>\n";
	$type = $FAM;
    }
    if ($mapFile eq "")
    {
	print "<P><B>ERROR</B>: User-provided map expected.<P>\n";
	$nerr++;
    }
} 
else		# $dbType ne $USER 
{
    if ($type eq $FAM)
    {
       if ($dbType eq $PRINT) { $mapFile = "$PRINTS_DATABASE/$name.map"; }
       else                   { $mapFile = "$BLOCKS_DATABASE/$name.map"; }
    }
    else	# $type eq $SEQ, use the wais index to find the families
    {
      $x = "'$name'";
      open(WAIS, "-|") || exec ("./waisq", "-c", $waisd, "-s", $waisd, "-m", "10", "-f", "-", "-S", "blocks.src", "-g", $x);

      #         :headline "IPB001438   /home/blocks/wais_files/"
      #         :headline "PR00722   /home/blocks/wais_files/"
      while ($in_rec = <WAIS>)
      {
         if ( ($in_rec =~ m/headline "IPB0/) ||
              ($in_rec =~ m/headline "PR0/) )
         {
            ($x) = $in_rec =~ m/headline "(\S+)/;
            push (@fams, $x);
         }
      }
      close(WAIS);
      $mapFile = "$tmpdir/$$.map";
      open(OUT, ">$mapFile");
      foreach $x (@fams)
      {
         #  am assuming these files exist and are formatted perfectly
         if ($x =~ m/^IPB/)
         { open(MAP, "< $BLOCKS_DATABASE/$x.map"); }
         elsif ($x =~ m/^PR/)
         { open(MAP, "< $PRINTS_DATABASE/$x.map"); }
         $in_rec = <MAP>;
         print OUT "$in_rec";	# should be >$x ...
         $done = 0;
         while (($done==0) && ($in_rec = <MAP>))
         {
            if ($in_rec =~ m/$name/)
            {
               $done = 1;
               print OUT "$in_rec";
               ($in_name, $in_len, $in_nblk) = split(/\s+/, $in_rec);
               for ($ib=0; $ib < $in_nblk; $ib++)
               {
                  $in_rec = <MAP>;
                  print OUT "$in_rec";
               }
            }

         }
         close(MAP);
      }
      close(OUT);
    } # end of SEQ
} # end of not USER

if (!(-e $mapFile))
{
   print "<P><B>ERROR</B>: Map file does not exist $mapFile<P>\n";
   $nerr++;
}

# Read $mapFile into memory 
my @maprecs;
open (INPUT, "<$mapFile");
@maprecs = <INPUT>;
close(INPUT);

if ($#maprecs < 1)
{
   print "<P><B>ERROR</B>: Map file is empty $mapFile<P>\n";
   $nerr++;
}
if ( !($maprecs[0] =~ m/^>/) )
{
   print "<P><B>ERROR</B>: Invalid map file $mapFile<P>\n";
   $nerr++;
}

if ($nerr > 0) { exit(0); }
#print "<PRE>";
#foreach $x (@maprecs) { print "$x\n"; }
#print "</PRE>";

#--------------------------------------------------------------------
#	Load the colors into a list
open(COLORS, $colorsFile);
my @colors = <COLORS>;
close (COLORS);
#print "<PRE>";
#my $x;
#foreach $x (@colors) { print "$x\n"; }
#print "</PRE>";

#--------------------------------------------------------------------
#	These routines use global vars $isCondensed and @maprecs
printHeader();
if ($type eq $SEQ)
{ process_seq(); }
else
{ process_fam(); }
printFooter();

exit(0);

#==================================================================
#------------------MAIN PROCESSING SUBROUTINES-------------------------

# reads in a Block Family and will display a block map 
# of the Sequences that it belongs to
#  Assumes global vars @colors, $isCondensed and @maprecs
#>IPB000104 3 9
#Type I antifreeze protein signature
#ANP4_PSEAM|P02734 85 3
#A 45 59
#B 60 70
#C 71 79
#ANPX_PSEAM|P07835 91 3
#A 51 65
#B 66 76
#C 77 85
#...

sub process_fam()
{
    my $input;
    my $nseqs;
    my $nblocks;
    my $outfile;
    my $famID;
    my $description;
    my ($minrecs,$is,$ib,$im, $seqID, $seqlen, @seqstarts, @seqends, @blocknames);

    ($famID,$nblocks,$nseqs) = split(/\s+/, $maprecs[0]);
    $description = $maprecs[1];
    chomp ($description);
    printFileDescription($description,$nseqs,$nblocks);

    $minrecs = 1 + $nseqs * ($nblocks + 1);
    if ($#maprecs < $minrecs)
    {
       print "<P><B>WARNING:</B> Map file appears to be incomplete.<P>\n";
    }

#   for ($is=0; $is < $nseqs; $is++)
#   {
#      $im = 2 + $is * ($nblocks + 1);
#      #seq info is in $maprec[$im}
#      for ($ib=0; $ib < $nblocks; $ib++)
#      {
#         $im++;
#         #block info for seq is in $maprec[$im}
#      }
#   }
    my @tempArray;
    my $sizeOfArray;
    my $isDoubleLinks;
    my $link;
    my $link2;
    my $seqIDUncut;
    
    #print table header
    printTableHeader();
    
    $im = 2;				# @maprecs index
    for ($is = 0; $is < $nseqs; $is++) 
    {
        #  Read a sequence line
	$input = $maprecs[$im++];
	($seqID,$seqlen,$nblocks) = split (/\s+/, $input);
	$seqIDUncut 	  = $seqID;
	
	#parse ID correctly for ChromDB and Swiss Prot/Trembl
	#parse second link
	@tempArray      = split(/\|/,$seqID);
	$sizeOfArray    = @tempArray;
	$isDoubleLinks  = 0;
	
        if ($dbType ne $USER)
        {
	    # Parsing Links for Sequence IDs
	    if ( $sizeOfArray > 1)
            {
		$link =  makeSwissProtLink((split(/\|/,$seqID))[1]);
	    }
	    else
            {
                $link = makeSwissProtLink($seqID);
            }
	} # for the $USER check;

	#  Get the block starts/ends in this sequence
	for ($ib = 0; $ib < $nblocks; $ib++)
        {
	    $input = $maprecs[$im++];
	    ($blocknames[$ib], $seqstarts[$ib], $seqends[$ib]) = split (/\s+/, $input);
	}
	# Display the blocks for this sequence

        print  "\t\t<tr>\n";
	if ($dbType eq $USER){
                print  "\t\t<td>$seqID\n";
	}
	elsif ($isDoubleLinks){
		print  "\t\t<td><font><a target= blank href=\"$link\">$tempArray[0]</a></font>/";
		print  "<font><a target= blank href=\"$link2\">$tempArray[1]</a></font>\n";
	}
	else{
       		print  "\t\t<td><font><a target= blank href=\"$link\">$seqID</a></font>\n";
	}
       	print  "\t\t<td>$seqlen\n";
	print  "\t\t<td NOWRAP>";

	if ($isCondensed){
		displaySequenceCondensed($seqID, $nblocks, $seqlen, \@blocknames, \@seqstarts,\@seqends, \@colors);	
	    }
	else {
		displaySequenceLong($seqID, $nblocks, $seqlen, \@blocknames, \@seqstarts,\@seqends, \@colors);
	    }
    }  # end of sequence is
    printTableFooter();
    
#   $outfile 	           = $famID . ".html";
#   printFastaOutput($mapFile);
} # end of process_fam


# Display all blocks from multiple block families for a single sequence
#	Expects global vars $name, isCondensed, @maprecs, @colors
#	Assumes @maprecs has the records for one sequence, e.g.
#>IPB008146 5 2
#Q9K134 472 5
#A 49 91
#B 129 138
#C 174 224
#D 266 278
#E 340 383
#>IPB008147 7 2
#Q9K134 472 7
#A 52 94
#B 128 138
#C 168 187
#D 207 254
#E 266 278
#F 317 351
#G 358 372

sub process_seq()
{
    my $currentIndex;
    my $input;
    my $seqID;
    my $seqlen;
    my $nfams;
    my $famID;
    my $famIDUncut;
    my $nblocks;
    my @blocknames;
    my @seqstarts;
    my @seqends;
    my ($if, $im, $ib, $temp);
    
    # print the description table at the top
    printFileDescriptionBySequence($name);

    # read maprecs once to see how many families
    $nfams = 0;
    foreach $input (@maprecs)
    {  if ($input =~ m/^>/) { $nfams++; }   }
    
    $im = 0;
    for ($if =0; $if< $nfams; $if++)
    {
        $input = $maprecs[$im++];
        chomp($input);
	printTableHeaderBySequence($name);
        ($famID, $nblocks) = split(/\s+/, $input);
        $famID =~ s/^>//;
	
        $input = $maprecs[$im++];
        ($seqID, $seqlen, $nblocks) = split(/\s+/, $input);
	    
        for ($ib=0; $ib < $nblocks; $ib++)
        {
	    $input = $maprecs[$im++];
            ($blocknames[$ib], $seqstarts[$ib], $seqends[$ib]) = split(/\s+/, $input);
	} # end of blocks for this family
        
        #  Display the blocks for this family
	    
	    print  "\t\t<tr>\n";
	    printSeqLink($famID);
	    print  "\t\t<td>$seqlen\n";
	    print  "\t\t<td>";
	    
	    if ($isCondensed){
		displaySequenceCondensed($famID, $nblocks, $seqlen, \@blocknames, \@seqstarts,\@seqends, \@colors);	
	    }
	    else {
		displaySequenceLong($famID, $nblocks, $seqlen, \@blocknames, \@seqstarts,\@seqends, \@colors);
	    }
	    
	printTableFooter();
	print "<BR>";
#	printFastaOutput($inputFile[$i],$type); 
	
    }  # end of a block family
    exit;

} # end of process_seq
#----------------------------------------------------------------------

#----------------HTML SUBROUTINES---------------------------------------

sub printSeqLink()
{
    my $ID = $_[0];
	    
    if(($ID =~/IPB.*/) || ($ID =~/PR.*/)){
	print  "\t\t<td><a target=blank href= \"http://blocks.fhcrc.org/blocks-bin/getblock.pl?$ID\">$ID</a></font>\n";
    }
} # end of printSeqLink

sub printErrorUserDefinedInput{
    # Print Error In formPage
    print "<HTML><TITLE>Error In Submission</TITLE><BODY BGCOLOR = WHITE>";
    printHeader();
    print "<b>Your input format was incorrect.  Be sure that all fields are labeled correctly.<br>";
    print "If you use dbtype = USER, please make sure there is type = FAM and the map field is filled in correctly.<br>";
    print "</BODY></HTML>";

}
sub printNoHitsFound{
     # Print Error Page
    print "<HTML><TITLE>Error In Submission</TITLE><BODY BGCOLOR = WHITE>";
    printHeader();
#   print "<b><font color= blue> @name</font> yielded no hits.<br>";
    print " Please go to the <a href= 'http://www.proweb.org/blockmap'> previous page</a> and resubmit";
    print "</BODY></HTML>";
}

sub printErrorPage{
    # Print Error Page
    print "<HTML><TITLE>No Results Found</TITLE><BODY BGCOLOR = WHITE>";
    printHeader();
#   print "<b>Your request for <font color= blue> @name</font> did not yield any results.<br>";
    print " Please go to the <a href= 'http://www.proweb.org/blockmap'> previous page</a> and resubmit another entry.";
    print "</BODY></HTML>";
}
sub printErrorInFormPage{
    # Print Error In formPage
    print "<HTML><TITLE>Error In Submission</TITLE><BODY BGCOLOR = WHITE>";
    printHeader();
    print "<b>Your input format was incorrect.  Be sure that all fields are labeled correctly.<br>";
    print "The correct labels are: name, type, dbtype, and condensed.<br>";
    print " Please go to the <a href= 'http://www.proweb.org/blockmap'> previous page</a> to resubmit or type directly as a URL";
    print "</BODY></HTML>";
}
sub printErrorInputPage{
    # Print Error Page
    print "<HTML><TITLE>Error In Submission</TITLE><BODY BGCOLOR = WHITE>";
    printHeader();
    print "<b>Not all fields are filled in correctly.  Make sure that all fields have a valid value.<br>";
#   print "Also, be aware that the HITS option for type is only valid with the Chrom Database <br>";
    print " Please go to the <a href= 'http://www.proweb.org/blockmap'> previous page</a> to resubmit or type directly as a URL";
    print "</BODY></HTML>";
}


sub printFastaOutput
{
    my $inputFile = $_[0];
    my $input;

    print "\t<!-- The following is the Fasta input used to make this display:\n";

    open(INPUT, "<$inputFile");

    while ($input = <INPUT>){
	print "\t\t$input";
    }
    close (INPUT);
    print "\t------------------------------------------------------------------>";
}  # end of printFastaOutput


sub printFileDescriptionBySequence {

    my $name    = $_[0];

    #print the description table at the top
    print  "\t<table border = 1>\n";
    print  "\t\t<tr>\n\t\t<td><font face=\"Verdana\" size=\"2\"><b>SequenceID(s):</b></font> </td>\n";
    print  "\t\t<td><font face=\"Verdana\" size=\"2\"> $name </font></td>\n\n";

    printLegend();
    print  "\t</table>\n\n\n";
    print  "\t<br><br>\n\n";
}


sub printTableHeaderBySequence(){
	my $sequenceName = $_[0];
    print  "\t<table border = 1>\n";
    print  "\t\t<tr>\n\t<td align = \"center\"><font face=\"Verdana\" size=\"2\"><b>Block Family</b></font></td>\n";	
    print  "\t\t<td align = \"center\"><font face=\"Verdana\" size=\"2\"><b>Length</b></font></td>\n";	
    print  "\t\t<td align = \"center\"><font face=\"Verdana\" size=\"2\"><b>Sequence [$sequenceName]</b></font></td>\n\n";	
}

sub printFileDescription {

    my @inputVariables 		= $_;
    my $description 		= $_[0];
    my $nseqs 	= $_[1];
    my $blocks 			= $_[2];

    #print the description table at the top
    print  "\t<table border = 1>\n";
    print  "\t\t<tr>\n\t\t<td><font face=\"Verdana\" size=\"2\"><b>Description:</b></font> </td>\n";
    print  "\t\t<td><font face=\"Verdana\" size=\"2\"> $description </font></td>\n\n";
    print  "\t\t<tr><td><font face=\"Verdana\" size=\"2\"><b>Sequences:</b></font> </td>\n";
    print  "\t\t<td><font face=\"Verdana\" size=2>$nseqs</font></td>\n\n";
    print  "\t\t<tr><td><font face=\"Verdana\" size=\"2\"><b>Distinct blocks:</b></font></td>\n";
    print  "\t\t<td><font face=\"Verdana\" size=\"2\">$blocks</font>\n";
    printLegend();
    print  "\t</table>\n\n\n";
    print  "\t<br><br>\n\n";
}

sub printLegend{
    my $scalingWidth = 100;
    my $displayHeight = $DISPLAYHEIGHT;

    print  "\t\t<tr><td><font face=\"Verdana\" size=\"2\"><b>Map Scaling:</b></font> </td>\n";
    print  "\t\t<td><font face=\"Verdana\" size=2>\|";
    displayLine($scalingWidth, $displayHeight);
    print  "\| [100 amino acids]</font></td>\n\n";
    print  "\t\t<tr><td><font face=\"Verdana\" size=\"2\"><b>Notes:</b></font></td>\n";
    print  "\t\t<td><font face=\"Verdana\" size=\"2\"> Mouse over to show start and end positions </font>\n";
    
}

sub printHeader{
    print  "<HTML>\n<TITLE>Block map</TITLE>\n\t<BODY BGCOLOR = white>\n\t<table>\n\t\t<tr>";
    print  "<td><p><img border=0 src= $blockPicture width=91 height=113></p>\n";
    print  "\t\t<td><font size=6 face=\"Times New Roman\"> BLOCK MAP</font></tr>\n";
    print  "\t</table>\n\n\n\t<hr>\n\n\n";
}

sub printFooter{
    print  "\n\n\n\t</BODY>\n";
    print  "</HTML>";
}

sub printTableHeader(){

    print  "\t<table border = 1 bordercolor = gray>\n";
    print  "\t\t<tr>\n\t<td align = \"center\"><font face=\"Verdana\" size=\"2\"><b>Sequence ID</b></font></td>\n";	
    print  "\t\t<td align = \"center\"><font face=\"Verdana\" size=\"2\"><b>Length</b></font></td>\n";	
    print  "\t\t<td align = \"center\"><font face=\"Verdana\" size=\"2\"><b>Blocks</b></font></td>\n\n";	
}

sub printTableFooter(){
	print "\n\t</table>\n\n";
}

sub displayLine {
    my $scaling_factor  = $SCALINGFACTOR;
    my $displayWidth    = $scaling_factor * $_[0];
    my $displayHeight   = $DISPLAYHEIGHT;

   # if the width is less than one, the line will display incorrectly
   # just ignore.  The size is negilible 
   if (!(($displayWidth > 0) && ($displayWidth < 1))){
	print  "<img border=0 src=$linePicture width=$displayWidth height=$displayHeight>";
   }
}

sub displayBlock {
    my $scaling_factor = 1/2; 
    my $seqstarts  = $_[0];
    my  $seqends   = $_[1];
    my  $colors    = $htmlUrl. "/". $_[2];
    my  $displayWidth  = $scaling_factor *$_[3];
    my  $motifName     = $_[4];
    my  $ID            = $_[5];
    my  $displayAlt    = $ID. $motifName ." ". "[$seqstarts-$seqends]";
    my  $displayHeight = $DISPLAYHEIGHT;

    # For now we just anchor the page to itself.  
    print  "<a href=# ";
    displayMouseOver($displayAlt);
    print ">";
    print "<img alt='$displayAlt' border=1 src=\"$colors\" width=$displayWidth height=$displayHeight></a>";

}
    
sub displayMouseOver{
        my $displayString = $_[0];

        print "onMouseOver=\"window.status='$displayString'; return true\"";
        print " onMouseOut=\"window.status=''; return true\"";
}

#----------------------------------------------------------------------




#----------------OTHER SUBROUTINES--------------------------------------

#Creates a swiss/trembl link given a swiss/trembl ID
sub makeSwissProtLink 
{
    my $displayID = $_[0];
    my $link;
    my $swissUrl = "http://www.expasy.ch/cgi-bin/get-sprot-entry?";

    $link = $swissUrl . $displayID;
    return $link;
}  # end of make SwissProtLink

#-----------------------------------------------------------------
# function displays series of blocks on multiple lines
sub displaySequenceLong 
{
    my $ID 				= $_[0];
    my $nblocks	= $_[1];
    my $seqlen           = $_[2];
    
    my $blocknamesRef	 	        = $_[3];
    my $seqstartsRef 		= $_[4];
    my $seqendsRef		        = $_[5];
    my $colorsRef		        = $_[6];		
    
    my @blocknames 			= @$blocknamesRef;
    my @seqstarts 		        = @$seqstartsRef;
    my @seqends			= @$seqendsRef;
    my @colors			= @$colorsRef;	

    my $displayWidth;
    my $j;

    for ($j =0; $j < $nblocks; $j++){
	
	#begin graphics for sequence display
	print  "<font size=1>$ID$blocknames[$j]</font>     ";
	displayLine ($seqstarts[$j]);
	
	$displayWidth = $seqends[$j] - $seqstarts[$j];
	displayBlock ($seqstarts[$j], $seqends[$j], $colors[$j], $displayWidth, $blocknames[$j], $ID);
	
	while (($nblocks >= $j+1) && ($blocknames[$j] eq $blocknames[$j+1])){
	    
	    $j++;
	    $displayWidth = $seqstarts[$j+1] - $seqends[$j];
	    displayLine ($displayWidth);
	    
	    $displayWidth = $seqends[$j] - $seqstarts[$j];
	    displayBlock ($seqstarts[$j], $seqends[$j], $colors[$j], $displayWidth, $blocknames[$j],$ID );
	}
	
	$displayWidth = $seqlen - $seqends[$j];
	if ($displayWidth){
	    displayLine ($displayWidth);
	    print  "<br>";
	}
    }
} # end of displaySequenceLong

#-----------------------------------------------------------------------
# function to display a series of blocks on one line
sub displaySequenceCondensed
{
    my $ID 				= $_[0];
    my $nblocks	= $_[1];
    my $seqlen           = $_[2];

    my $blocknamesRef	 	        = $_[3];
    my $seqstartsRef 		= $_[4];
    my $seqendsRef		        = $_[5];
    my $colorsRef		        = $_[6];		
    
    my @blocknames 			= @$blocknamesRef;
    my @seqstarts 		        = @$seqstartsRef;
    my @seqends			= @$seqendsRef;
    my @colors			= @$colorsRef;

    my @newSequenceStart;
    my @newSequenceEnd;
    my @newColorImage;
    my @newMotifNames;
    
    my $displayWidth;
    my @temp2;
    my @maps;
    my $j;
    my $k;
    my $l;
    my $m;
    my $isOverlap;

    # sort images/colors by ascending starting sequences
    @temp2 = sort {$a<=>$b} (@seqstarts);
    for ($k=0; $k < $nblocks; $k++){
	for ($l=0; $l < $nblocks; $l++){
	    #print " test: $temp2[$k] eq $seqstarts[$l] \n";
	    if ($temp2[$k] eq $seqstarts[$l]){
		$maps[$k] = $l;
	    }
	}
    } 

    # new arrays contain *actual* positions
    # non "new" ones will be manipulated according to adjustments in overlap 
    for ($m=0; $m < $nblocks; $m++){
	$newSequenceStart[$m] = $seqstarts[$maps[$m]];
	$newSequenceEnd [$m]  = $seqends [$maps[$m]];
	$newColorImage [$m]   = $colors [$maps[$m]];
	$newMotifNames [$m]   = $blocknames [$maps[$m]];
    }
    
    for ($m=0; $m < $nblocks; $m++){
	$seqstarts[$m]  = $newSequenceStart[$m];
	$seqends[$m]    = $newSequenceEnd[$m];
	$colors[$m]     = $newColorImage[$m];
	$blocknames[$m]     = $newMotifNames[$m];
    }
    
    #being graphics for sequence display
    displayLine ($seqstarts[0]);
    
        
    for ($j =0; $j < $nblocks; $j++){
	#check for overlap 
	$isOverlap = 0;
	if (($nblocks > $j+1) && ($seqstarts[$j+1] < $seqends[$j])){
	    $isOverlap = 1;	
	    $seqends[$j] = $seqstarts[$j+1];
	}
	
	$displayWidth = $seqends[$j] - $seqstarts[$j];
	displayBlock ($newSequenceStart[$j], $newSequenceEnd[$j], $colors[$j], $displayWidth, $blocknames[$j], $ID);
	
	if (!$isOverlap){
	    #if last sequence display to last sequence; else display line width
	    if ($nblocks > $j+1) {$displayWidth = $seqstarts[$j+1] - $seqends[$j];}
	    else {$displayWidth = $seqlen   - $seqends[$j];}
	    
	    if ($displayWidth){
		displayLine ($displayWidth);
	    }
	}
    }
} # end of displaySequenceCondensed
