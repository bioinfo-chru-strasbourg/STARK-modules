#! /usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];

$num_split = @ARGV[1];
$suffix = @ARGV[2];
open( FASTA, $file ) || die("can't open");

$num_files = -1;
$num_dirs  = -1;

while (<FASTA>) {
	
	if ( $_ =~ /Sear/ ) {
		close (homologs);
		$num_files = $num_files + 1;
		if ( $num_files % $num_split == 0 ) {
			$num_dirs = $num_dirs + 1;
			mkdir "partition\_$suffix\_$num_dirs";			
			$outfile = "homologs\_$suffix\_$num_dirs\_$num_files.fasta";
			open( homologs, ">>partition\_$suffix\_$num_dirs/$outfile" );
		}
		else{
			
			$outfile = "homologs\_$suffix\_$num_dirs\_$num_files.fasta";
			open( homologs, ">>partition\_$suffix\_$num_dirs/$outfile" );	
		}
			
	}
	if (!($_ =~ /Sear/ || /^NP_/ || /^XP_/)){
		print homologs $_;
	}
	
	
}

