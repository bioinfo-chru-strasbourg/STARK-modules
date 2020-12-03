#!/usr/bin/perl

#	rpsblast.pl
# Execute rpsblast & display the output

# 4/13/07 Execute getblock.pl instead of getblock.sh

#	Set file permissions to rw-rw----
system("umask 006");
$bin = "./";
$tmpdir = "../tmp";
$bdir = "../";
$ipbdir = "$bdir/data-blocks";
$prdir = "$bdir/data-prints";
$pfdir = "$bdir/data-pfam";
$bpdir = "$bdir/data-prodom";
$dmdir = "$bdir/data-domo";


# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>RPS-BLAST Results</TITLE>\n";
if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

if ($names{database} =~ /minus/) { $db = "blminus"; }
elsif ($names{database} =~ /plus/) { $db = "blplus"; }
elsif ($names{database} =~ /prints/) { $db = "prints"; }
else { $db = "blplus";  }

#print "$names{database} $db\n";
if ($db eq "") { $db = "blplus"; }
#$db = $names{database};
#$db =~ /\s+//g;

if ($names{ex} eq "") { $names{ex} = "5"; }
$ex = $names{ex};

#	Get the sequence & write it to a file
if ($names{sequence} eq "") {
   print "<H1>Error</H1> Please enter a query sequence.<P>\n";
   exit;
}
$seq = "../tmp/$$.seq";
open(SEQ, ">$seq");
print SEQ $names{sequence};
print SEQ "\n";
close(SEQ);

#	Filter flag, default is not to filter
$filt = "F";
if ($names{filter} eq "T") { $filt = "T"; }

#=========================================================================
#	Run shell now

$out = "../tmp/$$.out";
$arguments = $seq.' '.$db.' '.$ex.' '.$filt;
#print $arguments;
#	Only the first 2 arguments are getting passed this way
system("$bin/rpsblast.csh $seq $db $ex $filt > $out 2>&1");

#	why doesn't umask take care of this?
system("chmod 660 ../tmp/$$.*");

#=========================================================================
#	Format the output

print "<HTML><BASE HREF=\"http://blocks.fhcrc.org\">\n";
print "<TITLE>Block Searcher RPS-BLAST Results</TITLE>\n";
print "<H1>Block Searcher RPS-BLAST Results</H1>\n";

print"<PRE>\n";

$ac = "";
@blist = (); $nblock = 0; $bline = ""; 
open(OUT, "<$out");
while ($_ = <OUT>)
{
  #   Put links on the list of hits
  if (/^IPB\d\d\d\d\d\d\S*/) {
    s|(^IPB\d\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^BL\d\d\d\d\d\S*/) {
    s|(^BL\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^DM\d\d\d\d\d\S*/) {
    s|(^DM\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^PD\d\d\d\d\d\S*/) {
    s|(^PD\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^BP\d\d\d\d\d\S*/) {
    s|(^BP\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^PR\d\d\d\d\d\S*/) {
    s|(^PR\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^PF\d\d\d\d\d\S*/) {
    s|(^PF\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  #   Mark the blocks in the alignment
  elsif ($_ =~ m/^>/) 
  {
     @words = split(/\s+/, $_);
     $ac = $words[0]; $ac =~ s/^>//;
     if ($ac =~ m/^IPB/) { $dir = $ipbdir; }
     elsif ($ac =~ m/^PR/) { $dir = $prdir; }
     elsif ($ac =~ m/^PF/) { $dir = $pfdir; }
     elsif ($ac =~ m/^BP/) { $dir = $bpdir; }
     elsif ($ac =~ m/^DM/) { $dir = $dmdir; }
#print "ac=$ac dir=$dir\n";

      # Find >$ac in $dir/cobbler.pros and read cobbler sequence name
      #     & starting aa ">$ac seq_name from xxx to yyy"
      @blist = (); $nblock = 0; $bline = ""; 
      $cob_ac = "";
      open(GREP, "grep \"^>$ac\" $dir/cobbler.pros |");
      while($cob = <GREP>)
      {
#print "$cob\n";
         ($cob_ac, $cob_seq, $junk1, $cob_start, $junk2, $cob_end)
		= split(/\s+/, $cob);
         $cob_ac =~ s/^>//;
         if ($cob_ac eq $ac)
         {
         # Find >$ac in $dir/maps.dat and read starting & ending locations of
         #     each block for cobbler sequence seq_name
         # NOTE: too slow to read $dir/maps.dat directly, thus getseq
            $mapfile = "$tmpdir/$ac.map";
            system("$bin/blimps-bin/getseq $ac $dir/maps.dat $mapfile > /dev/null");
#           open(MAP, "<$dir/maps.dat");
            open(MAP, "<$mapfile");
            while ( ($map = <MAP>) && !($map =~ m/^>$ac/) )
            { ;  }
            if ($map =~ m/^>$ac/)
            {
               @words = split(/\s+/, $map);
               $nblock = $words[1];
               while ( ($map = <MAP>) && !($map =~ m/^>/) &&
                 !($map =~ m/^$cob_seq/) )
               { ;  }
               if ($map =~ m/^$cob_seq/)
               {
                  $i = 0;
                  while( ($i < $nblock) && ($map = <MAP>) )
                  {
                     ($block, $start, $end) = split(/\s+/, $map);
                     @blist = (@blist, $block, $start, $end);
                     $i++;
#print "$block $start $end\n";
                  }
               }
            }
            close(MAP);
            system("rm -f $mapfile");
         }
      }  # end of cobbler file
#print "cobac=$cob_ac $cob_seq $cob_start $cob_end\n";
  } # end of ^>

  elsif ($_ =~ m/^Sbjct: /) 
  {
  # Read each Sbjct line until next ^> line "Sbjct: <start>"
  #     Subtract (xxx-1) from <start>, compare it with block start/end
  #     If Sbjct line is within a block, insert line afterwards
     $bline = ""; $bpos = 0;
     if ($nblock > 0)
     {
        ($junk1, $sstart, $sseq, $send) = split(/\s+/, $_);
        $bsstart = $sstart + $cob_start - 1;	# Sbjct location in block
        $bsend = $send + $cob_start - 1;	# Sbjct location in block
        $gaps = 59 - ($bsend - $bsstart);	# Number of gaps in Sbjct
        $bline = "$ac";

        # distance between ^Sbjct: and start of sequence varies, us. 11 or 12
        # try to figure out what it is
        $subseq = substr($_, 6, 10);	# 10 chars after ^Sbjct:
        $sublen1 = length($subseq);
#print "subseq=$subseq\n";
        $subseq =~ s/[0-9]//g;		# remove numbers & spaces
#print "subseq=$subseq\n";
        $subseq =~ s/\s+//g;
#print "subseq=$subseq\n";
        $sublen2 = length($subseq);
        $skip = 6 + ($sublen1 - $sublen2); # Sbjct: + numbers & spaces
#print "sublen1=$sublen1 sublen2=$sublen2 skip=$skip\n";

        for ($i=length($ac); $i < $skip; $i++)	# Assumes "^Sbjct: nnn  "
        { $bline = "$bline" . " "; }
        $bpos = $bsstart; $sbjctpos = 0;
        for ($i=0; $i<$nblock; $i++)
        {
           $k = $i * 3;			# 3 elements in @blist
           ($block, $start, $end) = ($blist[$k], $blist[$k+1], $blist[$k+2]);
#print "bsstart=$bsstart bsend=$bsend block=$block start=$start end=$end\n";
           if ($end >= $bpos && $start <= $bsend)	# block touches
           {
              if ($end > $bsend) { $last = $bsend; }
              else { $last = $end; }
              if ($start > $bpos)
              {   
                 # If there are gaps in Sbjct before $start, have to skip over
		 #>>> wrong if gaps within the block <<<
                 $nextgap = index($sseq, "-", $sbjctpos); # -1 if none, else pos
                 if ($nextgap >= 0)
                 {
                    $sublen = $start - $bpos;
                    $subseq = substr($sseq, $sbjctpos, $sublen);
#print "subseq=$subseq\n";
                    $sublen1 = length($subseq);
                    $subseq =~ s/-//g;		# remove the gap characters
                    $sublen2 = length($subseq);
                    $ngaps = $sublen1 - $sublen2; # number of gaps in subseq
#print "start=$start bpos=$bpos sbjctpos=$sbjctpos sublen=$sublen ngaps=$ngaps\n";
                 }
                 for ($pos = $bpos; $pos < $start + $ngaps; $pos++)
                 { $bline = "$bline" . " "; $bpos++; $sbjctpos++; }
                 $bpos -= $ngaps;  # should be true block position
              }
#print "bpos=$bpos end=$end bsend=$bsend last=$last\n";
              for ($pos = $bpos; $pos <= $last; $pos++)
              { $bline = "$bline" . "$block"; $bpos++; $sbjctpos++; }
           }
        }  # end of blocks
        $bline = "$bline" . "\n";
     } # end of if blocks
  }
  print $_;
  if ($bline ne "") { print $bline; $bline = ""; }
}
close(OUT);

print "</PRE>";
print "</PRE><A HREF=\"#top\">[return to top]</A><P>";
#-------------------------------------------------------------------------
exit (0);


#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:	CONTENT_TYPE
#		QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the key.

# WARNING:  Some of this routine is program-dependent!!!

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#		Content-Disposition: form-data; name="key1"
#		<blank line>
#		val1
#		<boundary>
#		Content-Disposition: form-data; name="key2"
#		<blank line>
#		val2
#		<boundary>

sub general_parse {
  local($content_type, $query_string) = @_;
  local(%ans, @q, $pair, $loc, $boundary, $temp, $temp1);

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
#		Why is this necessary? (boundary= doesn't match actual)
	$boundary = "--".$temp;
#print "<PRE>";
#print "$query_string\n";
#print "boundary=$boundary\n";
        # break up into individual name/value lines
        @q = split(/$boundary/, $query_string);
        foreach $pair (@q) {
          # break the name/value pairs up
#print "pair=$pair\n";
          $loc = index($pair, "name=");
	  $temp = substr($pair, $loc+5);
#	  $loc = index($temp, "\n\n");
 	  $loc = index($temp, "\n");
	  $temp1 = substr($temp, $loc+2);
#print "1temp=$temp\ntemp1=$temp1\n";
#		Get rid of stuff after the name; including semicolon if any
	  $loc_semi = index($temp, ";");
	  $loc_eol = index($temp, "\n");
	  $loc = $loc_eol;
          if ($loc_semi > 0 && $loc_semi < $loc) {$loc = $loc_semi; }
	  if ($loc > 0) { $temp = substr($temp, 0, $loc); }
#		Get rid of quotes around the name
          $temp =~ s/\"//g;
#		Still has a trailing whitespace character ...
          $temp =~ s/\s//g;
#		Need to strip leading/ending whitespace off of $temp1,
#		but be careful not to strip off internal CRs
#>>>>>       Following line is program-dependent !!!
	  if ($temp ne "dbfile" && $temp ne "sequence" && $temp ne "database")
 	  { $temp1 =~ s/\s//g; }
#		MAC file lines end in just \r, no \n, so makelis won't find all
#		of the sequences; DOS file lines end in \r\n, UNIX in \n.
	  if ($temp eq "seqfile")
#	  { $temp1 =~ s/\r/\n/g; $temp1 =~ s/\n\n/\n/g; }
#		Change \r\n to \n and then \r to \n
 	  { $temp1 =~ s/\r\n/\n/g; $temp1 =~ s/\r/\n/g; }
#print "2temp=$temp\ntemp1=$temp1\n";
	  if ($temp ne "") { $ans{$temp} = $temp1; }
        }
     }
     else
     {  print "Cannot parse\n";  
        print "content_type=$content_type\n";
        print "query_string=$query_string\n";
     }
  }
  return %ans;
#print "</PRE>";
}   # end of general_parse
#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

