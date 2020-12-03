#!/usr/local/bin/perl -w

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use strict;
my $chrfile = $ARGV[0];
my $chrfile_name;
my $pid;
if ($chrfile =~ /.+?(\d+\.chrfile)/){
	$chrfile_name = $1;
}
if ($chrfile_name =~ /(\d+)\.chrfile/){
	$pid = $1;
}
my $tmp = "/opt/www/sift/tmp/$pid";
my $tempfile = "$tmp/$chrfile_name.temp";

open (TEMPFILE, ">$tempfile") || die ("cannot open tmp file");
while (<>){
	chomp;
	$_ =~ s/^\s+|\s+$//g;
	my $chr = "";
	my $start = "";
	my $stop = "";
	my $orn = "";
	my $allele = "";
	my $comment = "";
	my @elts = split /,/, $_;
	$chr = $elts[0];
	$start = $elts[1];
	$stop = $elts[2];
	$orn = $elts[3];
	$allele = $elts[4];
	$comment = $elts[5];
	next if ($chr eq "");
	#rare bad case if stop < start
	if ($start > $stop){
		next;
	}	

	#if insertion
	if ($start == $stop){
		if ($allele =~ /(\w+)/){	#if already in snp_classifier format (no slash): ATGGC
			$allele = uc ($1);
		}
		elsif($allele =~ /^\-*\/(\w+)/){  #if of the form -/ATTGCA -> ATTGCA
			$allele = uc ($1)
		}
		elsif ($allele =~ /^(\w+)\/\-*/){  #if of the form ATTG/- -> ATTG
			$allele = uc ($1);
		}
	}

	#if deletion
	if ($stop - $start >= 1){
		$allele = "-/";
	}
	if ($comment !~ /\w|\d/){
		print TEMPFILE "$chr,$start,$stop,$orn,$allele\n";
	}
	else{
		print TEMPFILE "$chr,$start,$stop,$orn,$allele,$comment\n";
	}
}
system ("mv $tempfile $chrfile");
