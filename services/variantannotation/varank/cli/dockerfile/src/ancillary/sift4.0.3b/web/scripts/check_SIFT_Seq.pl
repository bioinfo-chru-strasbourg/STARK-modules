#!/usr/local/bin/perl
#use strict;

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
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
$agent->timeout(300);
$agent->cookie_jar(
	HTTP::Cookies->new( file => "$ENV{HOME}/.cookies.txt", autosave => 1 ) )
  ;    # ~/.cookies.txt must be RW by the owner of this script!
push @{ $agent->requests_redirectable }, 'POST';

# we have to get redirected by the server due to a strange and invalid character contained in the uidbadge
# plus the need to extract a couple of state values to pass to the login script
my $pre_response =
  $agent->get('http://sift.jcvi.org/sift-bin/SIFT_seq_submit2.pl');

if ( !$pre_response->is_success ) {
	die $pre_response->status_line
	  ;    # if it fails at this point, it's likely a network or server issue.
}

my $req = $pre_response->request; # needed for getting the proper URL to post to


# build the form to send. the server is quite anal about what it will take.
my $script_response = $agent->post(
	$req->uri,
	[
		'queryseq' => ">gi|532319|pir|TVFV2E|TVFV2E envelope protein\nELRLRYCAPAGFALLKCNDADYDGFKTNCSNVSVVHCTNLMNTTVTTGLLLNGSYSENRTQIWQKHRTSNDSALILLNKHYNLTVTCKRPGNKTVLPVTIMAGLVFHSQKYNLRLRQAWC",
		'database' => "sp_tr",
		'info' => "3.00",
		'seq_identity_filter' => "90"

	]
);
$positions_found= 0;
$matrix_found = 0;
$pred_found = 0;

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
	
	if ($line =~ /ERROR! Unable to contact NCBI/i){
		print "SIFT_Seq Failed!";
		last;
	}
	if ($line =~ /timeout/i){
		print "SIFT_Seq Failed!";
		last;
	}
	
	if ($line =~ /Positions 1 to 100/){
		$positions_found = 1;
	}
	if ($line =~ /Scaled Probabilities for Entire Protein/){
		$matrix_found = 1;
	}
	if ($line =~ /Predictions of substitutions entered/){
		$pred_found = 1;
	}	
}
if ($positions_found == 0 || $matrix_found ==0 || $pred_found == 0){
	#print "SIFT_Seq Failed!";	
	open ( OUTFILE, ">/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_Seq_Failed.html");
	print OUTFILE $content;
}
else{
	#print "SIFT_Seq Passed!";
	open ( OUTFILE, ">/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_Seq_Passed.html");
	print OUTFILE $content;
}




