#!/usr/bin/perl

@spam_list = ("ku.name");
$email = "ISH\@ku.name";
($e1,$e2) = split(/\@/, $email);
print "email=$email e1=$e1 e2=$e2\n";
foreach $spam (@spam_list)
{
    if ($e2 eq $spam) { print "e2=$spam\n"; }
    if ($email =~ m/($spam)/) { print "$spam\n"; } 
#   if ($email =~ m/ku\.name/) { print "$spam\n"; } 
#   if ($email =~ m/ku/) { print "ku\n"; } 
#   if ($email =~ m/ISH/) { print "ISH\n"; } 
#   if ($email =~ m/name/) { print "name\n"; } 
}
exit;

for ($i=0;$i<10;$i++)
{
   $randnum = int(rand(5));
   print "$randnum\n";
}
exit;

$de_rec = "DE   one two three four\n";
($de) = $de_rec =~ m/^DE   (.+)/;
print "de=$de\n";
exit;

$mailprog = "/usr/bin/mailx";
$webmaster = "jorja\@fhcrc.org";
$siftmaster = "sift\@fhcrc.org";
$tmp = "../tmp";

$names{topic} = "ch";
$names{email} = "henikoff\@yahoo.com";
$com = "$tmp/test";


if ($names{topic} eq "si")
{
   system("$mailprog -s \"comment $names{topic}\" -r $names{email} $siftmaster < $com");
}
else
{
   print "$names{topic} $webmaster\n";
   system("$mailprog -s \"comment $names{topic}\" -r $names{email} $webmaster < $com");
}
print "<P>Thank you for your comment.<P>\n";


#-------------------------------------------------------------------------
exit (0);
