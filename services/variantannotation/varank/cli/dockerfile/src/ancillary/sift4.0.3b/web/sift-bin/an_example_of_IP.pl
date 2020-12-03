#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use Data::Dumper;
#
# Nov 28, 2000 added option to email results
# 7/4/01 gap option is turned off permanently
#11/11/03 Added SIFT_queue stuff  JGH
# 3/8/04 Don't wait for search to finish, send them format.pl URL  JGH
#-------------------------------------------------------------------------

$| = 1;
require 'SIFT_subroutines.pm';
#	Set file permissions to rw-rw----
system("umask 006");
$bin = ".";
$tmp = "../tmp";
$pid = $$;
$program_call = $bin . "/SIFTING2.csh";
$return_address = "sift\@fhcrc.org";
my $ip_logfile = "IP_log";
# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<body bgcolor=white>\n";

if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

if ($names{address} ne "") {
	$address = $names{address};
}

# check that there's a query sequence
if ($names{gi_number} eq "") {
	print "<H1>Error</H1>Please enter a GI number.<P>\n";
	exit;
}

### substitution file #######
$subst_file = $tmp . "/$pid.substfile";
if ($names{substitution_file} eq "" && $names{substitutions} eq "") {
        $subst_file = "-";
}
if ($subst_file ne "-") { # have some substitutions
        open (SUBSTITUTION, ">$subst_file");
        if ($names{substitution_file} ne "") {
                $names{substitution_file} =~ tr/a-z/A-Z/;
                print SUBSTITUTION $names{substitution_file};
        }
        if ($names{substitutions} ne "") {
                $names{substitutions} =~ tr/a-z/A-Z/;
                print SUBSTITUTION $names{substitutions};
        }
        close (SUBSTITUTION);
}
## Check that this IP address hasn't been used too much
 	my $IP_address = $ENV{REMOTE_ADDR};
	my %IP_log_hash = make_hash ($ip_logfile);
	if ($IP_log_hash{$IP_address} > 10) {
		print "<H1>Your computer has exceeded its daily limit.</H1><BR>"; 
		print "Please download <A HREF=\"http://blocks.fhcrc.org/sift/SIFT.html\">SIFT software</A HREF> directly to your computer or <A HREF=\"http://blocks.fhcrc.org/blocks-bin/contact.pl\">contact</A HREF> us so that we can help you.  Thank you for using SIFT. <BR>";
		exit;
	} else {
		update_IP_logfile ($ip_logfile, $IP_address);
		# i need to put a lock on this file
	}

 	 
## SIFT prediction operations

$exp_option = 1;
#$info = $names{info};
$comments = "$tmp/$pid.comments"; 
$out = $tmp . "/$$.siftresults";
$seq_identity_filter = $names{seq_identity_filter};
$seq_identity_filter = trim ($seq_identity_filter);
chomp ($seq_identity_filter);


######  Calling the program #########

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font 
color=blue>T</font> results</center></H1></A><BR>\n";

$sequences_to_select = "BEST";
if ($names{sequences_to_select} =~ /ALL/) {
	$sequences_to_select = "ALL";
}
if ($address ne "") {
print "Results will also be mailed to $address.<BR><BR>\n";	
}
select (STDOUT); $| = 1;
	system("$program_call $pid NCBI_blink $names{gi_number} $out 0 $subst_file 0 $seq_identity_filter $sequences_to_select $address > $comments 2>&1");
	system ("echo $address >> /home/sift/email_log");


print "</PRE>";
# SIFT successful -- OUTPUT RESULTS
$commentscsh = $tmp . "/$pid.commentscsh"; 
$psiblastout = $tmp . "/$pid" . ".alignedfasta";
$seqno = `grep ">" $psiblastout | wc -l`;
$seqno--; # subtract the QUERY sequence

if ($seqno >= 0) {
print "<b>$seqno</b> sequences were selected to be closely related to your query sequence.<BR>\n";
 
print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$pid.msf+Alignment+PRE\" TARGET=_blank>";
print "PSIBLAST alignment of submitted sequences</A><BR>\n";

print "<A HREF=\"/sift-bin/catfile_notitle.csh?$tmp/$pid.alignedfasta+PRE\" TARGET=_blank>Alignment in FASTA format </A> (for modification)<BR>";


print_psiblast_alignment_desc();
print "<i>Please check the sequences that have been chosen.  If the sequences are too diverged from your query or the alignment is questionable, we suggest you modify the fasta-formatted file above and <A HREF=\"http://blocks.fhcrc.org/sift/SIFT_aligned_seqs_submit.html\" TARGET=_blank>resubmit</A>.</i><BR>"; 
  
 

print "<BR>";


$no_table_files = `ls $tmp/$$.aatable* | wc -l`;

tolerated_nottolerated_aminoacidtable($tmp,$no_table_files, $$);

$commentscsh = $tmp . "/$$.commentscsh";
$loc = rindex ($out, "/");
$file = substr ($out, $loc+1, 100);
print "<BR>\n";


$loc = rindex ($out, "/");
$file = substr ($out, $loc+1, 100);

print_score_desc ($file);
print "<BR>\n";

print_score_desc ($tmp,$file);
print "<i>If you received a warning that the sequences were not diverse enough, you can have SIFT choose more diverse sequences <A HREF=\"http://blocks.fhcrc.org/sift/SIFT_seq_submit2.html\">here.</A></i><BR><BR>"; 
}

system ("cat $tmp/$pid.*.error >> $tmp/$pid.err");
#errors or warnings that don't disrupt program
open (OUT, "<$tmp/$pid.err");
while ($_ = <OUT>) {
# don't want to print Jorja's messages, only my own.  so my error
# messages all have ERROR in the line, but don't print blastpgp error 
# for reading checkpoing file
      if (/ERROR/ || /WARNING/ || /\*\*\*/ ) {
            unless (/blastpgp/) {
		print;
      		}
	}
}
close(OUT);
# Pauline added printint error files
@error_files = `ls $tmp/$pid.*.error`;
if (@error_files > 0) {
	# there's an error
#	print "<b>ERROR!</b> Please <A HREF=\"/blocks-bin/contact.pl\" > email us </A> with your sequence if you do not understand the error. <BR><BR>";
	if (open (OUT, "$tmp/$pid.seq.query.globalX")) {
		close (OUT);
		print "<BR>";
		system ("sed s/X/-/g $tmp/$pid.seq.query.globalX > $tmp/$pid.aln"); 
		print "<A HREF=\"/pauline-bin/catfile.csh?btest-tmp/$pid.aln+Alignment+PRE\">";
		print "PSIBLAST search resullts</A><BR>\n";
	}
#	exit (-1);
#}

print "<HR>";


}

#-------------------------------------------------------------------------
exit (0);

#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment 
#variable
# returns: an associative array of name/value pairs.  The name is the 
#key.
sub parse_query {
  local($query_string) = @_;
  local(%ans, @q, $pair);
#print $query_string;
  # break up into individual name/value lines
  @q = split(/&/, $query_string);

  foreach $pair (@q) {
    # break the name/value pairs up
    # use split rather than regular expressions because the value may 
   # have
    #  newlines in it
    split(/=/, $pair, 2);                        

    # change '+' to ' '
    $_[1] =~ s/\+/ /g;

    # change the escaped characters (has to be after the split on '&' 
	# and '=')
    $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

    $ans{$_[0]} = $_[1];
  }

  return %ans;
}                                     
#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a 
# string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:   CONTENT_TYPE
#               QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the 
# key.

# WARNING:  Some of this routine is program-dependent!!!

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#               Content-Disposition: form-data; name="key1"
#               <blank line>
#               val1
#               <boundary>
#               Content-Disposition: form-data; name="key2"
#               <blank line>
#               val2
#               <boundary>

sub general_parse {
  local($content_type, $query_string) = @_;
  local(%ans, @q, $pair, $loc, $boundary, $temp, $temp1);

  if ($content_type eq "application/x-www-form-urlencoded")
  {
     # break up into individual name/value lines
     @q = split(/&/, $query_string);

     foreach $pair (@q) {
       # break the name/value pairs up
       # use split rather than regular expressions because the value 
	# may have
       #  newlines in it
       split(/=/, $pair, 2);

       # change '+' to ' '
       $_[1] =~ s/\+/ /g;

       # change the escaped characters (must be after the split on '&' 
# and '=')
       $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

       $ans{$_[0]} = $_[1];
     } #end of foreach $pair 

  } #end of if ($content_type) 
  else
  {
     $loc = index($content_type, "boundary=");
     if ($loc > 0)
     {
        $temp = substr($content_type, $loc+9);
#               Why is this necessary? (boundary= doesn't match actual)
        $boundary = "--".$temp;
        # break up into individual name/value lines
        @q = split(/$boundary/, $query_string);


        foreach $pair (@q) {
	 # break the name/value pairs up
          $loc = index($pair, "name=");
          $temp = substr($pair, $loc+5);
#         $loc = index($temp, "\n\n");
          $loc = index($temp, "\n");
          $temp1 = substr($temp, $loc+2);
#   Get rid of stuff after the name; including semicolon if any
          $loc_semi = index($temp, ";");
          $loc_eol = index($temp, "\n");
          $loc = $loc_eol;
          if ($loc_semi > 0 && $loc_semi < $loc) {$loc = $loc_semi; }
          if ($loc > 0) { $temp = substr($temp, 0, $loc); }
#               Get rid of quotes around the name
          $temp =~ s/\"//g;
#               Still has a trailing whitespace character ...
          $temp =~ s/\s//g;
#               Substitute spaces with nothing
#		Need to strip leading/ending whitespace off of $temp1,
#               but be careful not to strip off internal CRs
#               MAC file lines end in just \r, no \n, so makelis won't 
# find all
#               of the sequences; DOS file lines end in \r\n, UNIX in 
#\n.
#               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
#In my case, I want to keep the newlines in fields which have "file" or 
# 'seq"
# and remove newlines everywhere else.
	   if  ($temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/) 
 	   { $temp1 =~ s/\r\n/\n/g; $temp1 =~ s/\r/\n/g; }
	   # for other variables that are NOT files or file-like, remove extra
	   #whitespace
	    else { $temp1 =~ s/\s//g;} 
	if ($temp ne "") { $ans{$temp} = $temp1; }
     }# end of foreach
     } #end of if loc > 0
     else
     {  print "Cannot parse\n";
        print "content_type=$content_type\n";
        print "query_string=$query_string\n";
     }
  }
  return %ans;
#print "</PRE>";
}   # end of general_parse

# returns hash for a file, 2nd field is the key and the 3rd field
# is the value 4th field, is the delimiter
sub
make_hash
{
        my ($file) = @_;
        my %hash;
        open (HASH, $file) || die "can't open $file";
        my $line;
        while ($line = <HASH>) {
                chomp ($line);;
                if (exists ($hash{$line})) {
			$hash{$line}++;
		} else {
			$hash{$line} = 1;
		}
        }
        close (HASH);
        return (%hash);
}

sub
update_IP_logfile
{
	my ($queuefile, $IP_address) = @_;
	
	
	$lockqueuefile = "$queuefile.lock";
	
	# lockfile will wait until it can lock the file
	`./lockfile $lockqueuefile`;
	
	
	# append the address and command to the queue file
	open(FILE, ">>$queuefile");
	print FILE "$IP_address\n"; 
	close(FILE);

	chmod(0664, $queuefile);

	# remove the lock file
	unlink($lockqueuefile);

}

