#!/usr/bin/perl
#	cyrcaLama.pl
#	Executed by LAMA_search.c
#	Args are 1=lama-output 2=zscore 3=cyrca-output
# Output codes:
#       0 - success, Cyrca output not empty
#    Errors:
#       1 - wrong input
#       2 - No Connected graphs found
#       5 - Could not open a file or execute a program

#-==-=-=-=-=-=-=-=-=-=- Modules   -=-=-==-=-=-==-==-=-=-=-=-=-=

use strict;

#-=-=-=-==-=-=-=-=-=-=- Constants  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

my $CyrcaEnvFile = "./cyrca/perlib/CyrcaEnv.pm";

require $CyrcaEnvFile;

my $depth = 1;
my $tmpDir     = "../tmp"; 
my $tmpCyrca   =  $tmpDir . $$ . "cyrca_tmp";
my $tmpAnnotated= $tmpDir . $$ . "cyrca_anotated";
my $CyrcaOutputNotEmpty;

#-=-=-=-=-=--=--==-=-=-=-=- Parameters -=-=-=-=--=-=-=-=-=--=-=-=-

my ($lamaFile, $zscore, $outputFile) = @ARGV;

unless ($lamaFile && $outputFile) { exit (1);}

my $title            = "CYRCA output\n";
my $cyr;
my $MY_ADD;

#-=-=-=-=-=-=-=-=-=-=-= The flow of the program -=-=-=-=-=-=-==-=-=-

unless (open (OUTPUTFILE, ">$outputFile"))
{	print STDERR "ERROR: can't open $outputFile for writing\n"; exit (5);}


$cyr = "$cyrcaEnv::SingleBlockCyrca $cyrcaEnv::lamaOutput $lamaFile $depth $zscore";

unless (open (CYR, "$cyr |") ) 
{ print STDERR "open $cyr failed\n";	exit (5);    }
my @CyrArray = <CYR>;
close CYR;

my $line;
foreach $line (@CyrArray)
{
    if ($line =~ /Connected Graph:/)
    {
	$CyrcaOutputNotEmpty = 1;
	last;
    }
}

unless ($CyrcaOutputNotEmpty){ exit (2);}

unless (open (TMP, ">$tmpCyrca"))
{ print STDERR "ERROR: can't open $tmpCyrca\n"; exit (5); }
print TMP @CyrArray;

close TMP;
undef @CyrArray; #Just to realease the memory.

$MY_ADD = "$cyrcaEnv::AddInfo $CyrcaEnvFile $tmpCyrca html"; 

unless (open (OUT, "$MY_ADD |"))
{	print STDERR "ERROR: execution of $MY_ADD failed: $?\n";	exit(5); }

open (TMP, ">$tmpAnnotated") || exit (5);
print TMP <OUT>;
close TMP;
close OUT;

my $multipleGraphs = "$cyrcaEnv::findBlocksWithMultipleGraphs $tmpAnnotated";

unless (open (OUT, "$multipleGraphs | "))
{	print STDERR "ERROR: execution of $multipleGraphs failed: $?\n";	exit(5); }
my $line;

# print the output, sciping the info lines:
while ($line = <OUT>)
{
    if ($line =~/Add_info version/)
    {next;}
    if ($line =~/Used blocks database/)
    {next;}
    if ($line =~/blocks2pdb links/)
    {next;}
    if ($line =~/cyrcash version/)
    {next;}
    if ($line =~/Input file:/)
    {next;}
    if ($line =~/Depth of the search:/)
    {next;}
    if ($line =~/Z-score:/)
    {next;}
    print OUTPUTFILE $line;
    if ($line =~ /Blocks participating in multiple graphs:/)
    {
	last;
    }
}

print OUTPUTFILE <OUT>;

close OUT;
close OUTPUTFILE;

# delete temporary file
unlink  $tmpCyrca;
unlink  $tmpAnnotated;
exit (0);
 
 
 
