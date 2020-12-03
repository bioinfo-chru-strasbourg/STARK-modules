#!/usr/bin/perl
#       makeblocks.pl
#	executes protomat.csh
#
#  2/ 9/98 Added pdf versions of trees
#  2/12/98 Moved codehop link to blocks/bin/codehop.sh; removed flush of
#	   stdout, caused output problems...
#  3/ 8/98 Possibly get sequences from a file. general_parse()
#  5/27/98 Gapped blast link
#  11/20/98 Fixed bug with semicolons in general_parse()
#           Changed blocks-bin to blockmkr-bin
#	    Relative file name for $tmp
#  12/ 7/98 drawgram.csh instead of drawgram.*.sh|csh & newick.sh
#           logo.csh instead of logo.*.sh|csh
#   7/20/99 Link to proweb's 3dmotifs
#   6/ 6/00 Fixed Gibbs Blast links
#   8/27/01 3dblocks changes
#  11/29/01 Fixed Blast links
#   5/13/02 Be sure input sequences have different titles
#   7/19/02 Fix proweb 3dblocks link (runs on separate system now)
#           Add proweb tree viewer link
#           Fix blastp links
#   8/14/02 Output in separate target windows
#           Link to file with hints on saving output
#   9/17/02 3dblocks link again (are METHOD & SAVE parameters used?)
#   1/ 3/03 Make proweb treeviewer link a button
#   1/ 6/03 Pass proweb treeviewer blocks & pros
#   4/16/03 Added proweb block mapper
#   4/23/03 User reports Gibbs format problem (?), major format changes
#   4/28/03 Don't use $r.out file
#   5/ 9/03 Execute protomat.csh instead of processmail.sh
#   5/ 9/03 Put htmlize stuff in bm_format.pl
#=========================================================================
#
# "constants" and variables
#	Assume are executing in ~/bin subdirectory
$tmp = "../tmp/bm";
$www_address = "WWW";


#	Be sure files and directories are created with rw-rw---- permissions
system("umask 006");

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>Block Maker results</TITLE>\n";

# check that this is a POST submission
if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
    print "This script should be referenced with a METHOD of POST.\n";
    exit;
}
#	Check nixlist; IP addresses with a record of bad behavior
$ipaddr = $ENV{"REMOTE_ADDR"};
$nixstatus = system("./nixcheck.pl $ipaddr");
if ($nixstatus != 0)  {  exit(-1);  }

#print "<PRE>";
#print "<P>";
#print "HOST:$ENV{\"REMOTE_HOST\"}<BR>";
#print "ADDR:$ENV{\"REMOTE_ADDR\"}<BR>";
#print "USER:$ENV{\"REMOTE_USER\"}<BR>";
#print "IDENT:$ENV{\"REMOTE_IDENT\"}<BR>";
#print "$ENV{\"CONTENT_TYPE\"}";
#print "</P>";

#	Now multipart/form-data
#if ($ENV{"CONTENT_TYPE"} ne "application/x-www-form-urlencoded") {
#    print "This script can only be used to decode form results. \n";
#    exit;
#}


# get the QUERY_STRING environment variable and the name/value pairs
read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});


# debug
#print "<PRE>\n";
#print "<P>\n";
#print $QUERY_STRING ;
#print "<P>\n";
#exit;

#%names = &parse_query($QUERY_STRING);
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

#debug
#print "<PRE>\n";
#while (($key, $val) = each %names) {
#   print "key:$key\nval:$val\n";
#   $n = length($key); print " $n\n";
#   $loc = index($key, "address");
#   if ($loc >= 0) {$address = $val; }
#   $loc = index($key, "desc");
#   if ($loc >= 0) {$desc = $val; }
#}

#print "<PRE>\n";
#print "address:$names{address}\n";
#print "address:$address\n";
#print "desc:$names{desc}\n";
#print "desc:$desc\n";
#print "seqfile:$names{seqfile}\n";
#exit;


#
if ($names{desc} eq "")
{ $names{desc} = "unknown"; }

# check that there is sequence data
$temp = $names{sequences};
$temp =~ s/\s//g;
$n1 = length($temp);
$temp = $names{seqfile};
$temp =~ s/\s//g;
$n2 = length($temp);
if ($n1 <= 0 && $n2 <= 0)
{
    print "<H1>Error</H1>
	Please provide at least two protein sequences.<P>\n";
    exit;
}


# get $id
$id = $$;
system("mkdir -p $tmp/$id");
system("chmod -f 770 $tmp/$id");
$prefix = "$tmp/$id/$id";

$mail_flag = 0;
if ($names{address} ne "") { $return = $names{address}; $mail_flag = 1; }
else                       { $return = $www_address;    }
#	This file still looks like an email message
$seq_file = "$prefix.in";
open(SEQ, ">$seq_file");
print SEQ "From $return\n";
print SEQ "From: $return\n";
print SEQ "Reply-To: $return\n";
print SEQ "Subject: $names{desc}\n";
print SEQ "\n";

#	Notice that the sequences can come from both sources
#	Neither will ever eq "" because of whitespace ...
if  ($names{seqfile} ne "" )
{ print SEQ $names{seqfile}; }
if  ($names{sequences} ne "" )
{ print SEQ $names{sequences}; }
print SEQ "\n";
close SEQ;



#=========================================================================
# run protomat.csh, creates several output files in $prefix
#	This flush causes everything printed until OUT is referenced
#	again to be printed ... that is, stuff after the wait is printed
#	before it finishes waiting. It only seems to wait for stuff from
#	OUT and continues to process other statements. 
#	If email, queue it and exit
print "<PRE>";
if ($mail_flag)
{
   system("./add_queue_entry.pl MAKER_queue ./protomat.csh $prefix $return");
   print "Your results will be emailed to $return\n";
}
else
{
   print "Making your blocks, please wait...\n";
   select(STDOUT); $| = 1;
   system("./protomat.csh $prefix $return 2>&1 /dev/null");
}

#===========================================================================
$nerr = 0;
open(ERR, "<$prefix.err");
while ($err = <ERR>) { print "$err"; $nerr++; }
close(ERR);
if ($nerr > 0) { exit(-1); }

print "Your results will be available for approximately 24 hours ";
print "<A HREF=\"/blocks-bin/bm_format.pl?$id\">here</A>.<P>\n";

exit(0);

#=========================================================================
#-------------------------------------------------------------------------
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

# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:	CONTENT_TYPE
#		QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the key.

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
#print "temp=$temp\ntemp1=$temp1\n";
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
	  if ($temp ne "seqfile" && $temp ne "sequences")
 	  { $temp1 =~ s/\s//g; }
#		MAC file lines end in just \r, no \n, so makelis won't find all
#		of the sequences; DOS file lines end in \r\n, UNIX in \n.
	  if ($temp eq "seqfile")
#	  { $temp1 =~ s/\r/\n/g; $temp1 =~ s/\n\n/\n/g; }
#		Change \r\n to \n and then \r to \n
 	  { $temp1 =~ s/\r\n/\n/g; $temp1 =~ s/\r/\n/g; }
#         print "temp=$temp\ntemp1=$temp1\n";
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

