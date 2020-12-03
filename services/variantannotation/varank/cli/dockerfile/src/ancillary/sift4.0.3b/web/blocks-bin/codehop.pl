#!/usr/bin/perl
#			codehop.pl
#9/25/08 Check for codon usage table here

# Execute codehop & process the output

#	Set file permissions to rw-rw----
system("umask 006");

#	Location of codon usage tables
$codon_dir = '../docs';

$program = './codehop.csh';
$htmlize = './htmlize-codehop';
$bdir = '../tmp/';
$blocks = $bdir.$$.'.blks';
$out = $bdir.$$.'.pride';
$err = $bdir.$$.'.err';

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>CODEHOP Results</TITLE>\n";
#print "<PRE>";


if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

#	Get the parameter info
#	codonuse	-Ccodonuse
#	gcode		-Ggenetic_code
#	core_degen	-Dcore_degen
#	core_strict	-Score_strict
#	clamp_strict	-Lclamp_strict
#	clamp_temp	-Tclamp_temp
#	polyx		-Apolyx
#	clamp_conc	-Nclamp_conc
#	verbose		-V 
#	outoligo	-O
#	rose		-R
#	most		-M
#       begin	        -B
#
read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
#print $QUERY_STRING;

%names = &parse_query($QUERY_STRING);

#	Get the blocks & write them to a file for codehop
if ($names{blocks} eq "") {
   print "<H1>Error</H1> Please enter a block.<P>\n";
   exit;
}
open(BLK, ">$blocks");
print BLK $names{blocks};
print BLK "\n";
close(BLK);


#	Make a string of the other parameters
if ($names{codonuse} ne "")
{  $codon_file = "$codon_dir/$names{codonuse}"; }
else
{  $codon_file = "$codon_dir/default.codon.use"; }
if (!(-e $codon_file) || (-z $codon_file))
{  $codon_file = "$codon_dir/default.codon.use"; }
print "<PRE>names{codonuse}=$names{codonuse}\n";
print "codon_file=$codon_file\n";
$parameters = '-C'.$codon_file;

$rose = ''; $most = ''; $verbose = ''; $outoligo = ''; $begin = '';
if ($names{gcode} ne "")
{  $parameters = $parameters.' -G'.$names{gcode}; }

#if ($names{codonuse} ne "")
#{  $parameters = $parameters.' -C'.$codon_dir.$names{codonuse};
#   $codon = '-C'.$codon_dir.$names{codonuse};  }
#else
#{  $parameters = $parameters.' -C'.$codon_dir.'default.codon.use';
#   $codon = '-C'.$codon_dir.'default.codon.use';  }
#
if ($names{core_degen} ne "")
{  $parameters = $parameters.' -D'.$names{core_degen}; }
if ($names{core_strict} ne "")
{  $parameters = $parameters.' -S'.$names{core_strict}; }
if ($names{clamp_strict} ne "")
{  $parameters = $parameters.' -L'.$names{clamp_strict}; }
if ($names{clamp_temp} ne "")
{  $parameters = $parameters.' -T'.$names{clamp_temp}; }
if ($names{polyx} ne "")
{  $parameters = $parameters.' -A'.$names{polyx}; }
if ($names{clamp_conc} ne "")
{  $parameters = $parameters.' -N'.$names{clamp_conc}; }
#		Checkboxes
if ($names{rose} eq "TRUE")
{  $parameters = $parameters.' -R';  $rose = '-R'; }
if ($names{most} eq "TRUE")
{  $parameters = $parameters.' -M'; $most = '-M';  }
if ($names{verbose} eq "TRUE")
{  $parameters = $parameters.' -V'; $verbose = '-V'; }
if ($names{outoligo} ne "")
{  $parameters = $parameters.' -O'.$names{outoligo}; }
if ($names{begin} eq "TRUE")
{  $parameters = $parameters.' -B'; $begin = '-B'; }


print "<PRE>$parameters\n<\PRE>";

#=========================================================================
#	Run codehop now
#  NOTE: Following are unreliable, perl executes sh -c, sometimes hangs
#open(OUT, "$program $blocks $parameters  2>&1 | $htmlize |");
#open(ERR, "$program $blocks $out '$parameters' 2>&1 |");
#wait;
#while ($_ = <ERR>) {
#   print;
#}
#close(ERR);

print "<PRE>";
system("$program $blocks $out '$parameters' > $err 2>&1");
open (ERR, "<$err");
while ($_ = <ERR>) { print; }
close(ERR);
print "</PRE>";

#	Make sure group can remove these files
system("chmod -f 660 ../tmp/$$.*");
#system("rm $err");
#
open (OUT, "<$out");
while ($_ = <OUT>) { print; }
close(OUT);


#print "</PRE>";
exit;

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

