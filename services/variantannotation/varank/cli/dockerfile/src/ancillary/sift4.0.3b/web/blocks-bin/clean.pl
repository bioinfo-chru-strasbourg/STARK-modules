#!/usr/bin/perl
#
#	Remove all ../tmp text files > <days> old

#Default is 0.5 days = 12 hours
if ( $#ARGV >= 0 ) { $days = $ARGV[0]; }
else               { $days = 0.5; }
if ($days < 0.0 || $days > 1.0) { $days = 0.5; }
#print STDERR "$#ARGV days=$days\n";

$tmp = "$$.clean.tmp";

system("ls -t ../tmp > $tmp");

open(IN, "< $tmp");
foreach $x (<IN>)
{
   chomp($x);
   $file = "../tmp/$x";

   # Look for text files > $days old
   if ((-T $file) && ((-M $file) > $days))
   {  system("rm -f $file"); }
#  {  print "$file\n"; }
}
close(IN);
system("rm -f $tmp");
exit(0);
