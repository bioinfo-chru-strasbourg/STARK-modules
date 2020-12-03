#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

#use strict;
use HTTP::Status;
use HTTP::Response;
use HTTP::Cookies;
use LWP::UserAgent;
use URI::URL;
($DAY, $MONTH, $YEAR) = (localtime)[3,4,5];
$y = $YEAR + 1900;
$m = $MONTH + 1;
$d = $DAY;
$datestamp = "$y-$m-$d";
my $agent = LWP::UserAgent->new;

# build our user agent
$agent->timeout(60);
$agent->cookie_jar(
	HTTP::Cookies->new( file => "$ENV{HOME}/.cookies.txt", autosave => 1 ) )
  ;    # ~/.cookies.txt must be RW by the owner of this script!
push @{ $agent->requests_redirectable }, 'POST';

# we have to get redirected by the server due to a strange and invalid character contained in the uidbadge
# plus the need to extract a couple of state values to pass to the login script
my $pre_response =
  $agent->get('http://sift.jcvi.org/sift-bin/SIFT_dbSNP_submit.pl');

if ( !$pre_response->is_success ) {
	die $pre_response->status_line
	  ;    # if it fails at this point, it's likely a network or server issue.
}

my $req = $pre_response->request; # needed for getting the proper URL to post to


# build the form to send. the server is quite anal about what it will take.
my $script_response = $agent->post(
	$req->uri,
	[
		'SNP' => "rs1098"

	]
);
$rsid_found= 0;
$subst_found = 0;
$protein_found = 0;
$prediction_found = 0;
$content = $script_response->content;
chomp($content);
$content =~ s/\r//g;
@lines = split( "\n", $content );
for ( $i = 0 ; $i <= scalar @lines ; $i++ ) {
	$line      = @lines[$i];
	$prev_line = @lines[ $i - 1 ];
	$next_line = @lines[ $i + 1 ];
	$line      =~ s/\s+$//;
	$line      =~ s/^\s+//;
	$prev_line =~ s/\s+$//;
	$prev_line =~ s/^\s+//;
	$next_line =~ s/\s+$//;
	$next_line =~ s/^\s+//;

	if ( $line eq "" ) {
		next;
	}
	#print "$line\n";
	if ($line =~ /rs1098/){
		$rsid_found = 1;
	}
	if ($line =~ /L8Q/){
		$subst_found = 1;
	}
	if ($line =~ /NP_078974/){
		$protein_found = 1;
	}	
}
if ($rsid_found == 0 || $subst_found ==0 || $protein_found == 0){
	#print "SIFT_dbSNP Failed!";	
	open ( OUTFILE, ">/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_dbSNP_Failed.html");
	print OUTFILE $content;
}
else{
	print "SIFT_dbSNP Passed!\n";
	print "/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_dbSNP_Passed.html\n";
	open ( OUTFILE, ">/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_dbSNP_Passed.html") || die ("cannot open file");
	print OUTFILE $content;
}




