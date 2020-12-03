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
$tmp = "/opt/www/sift/tmp";
($DAY, $MONTH, $YEAR) = (localtime)[3,4,5];
$y = $YEAR + 1900;
$m = $MONTH + 1;
$d = $DAY;
$datestamp = "$y-$m-$d";
my $agent = LWP::UserAgent->new;
$agent->timeout(60);
# build our user agent

$agent->cookie_jar(
	HTTP::Cookies->new( file => "$ENV{HOME}/.cookies.txt", autosave => 1 ) )
  ;    # ~/.cookies.txt must be RW by the owner of this script!
push @{ $agent->requests_redirectable }, 'POST';

# we have to get redirected by the server due to a strange and invalid character contained in the uidbadge
# plus the need to extract a couple of state values to pass to the login script
my $pre_response =
  $agent->get('http://sift-dev.jcvi.org/sift-bin/SIFT_BLink_submit.pl');

if ( !$pre_response->is_success ) {
	die $pre_response->status_line
	  ;    # if it fails at this point, it's likely a network or server issue.
}

my $req = $pre_response->request; # needed for getting the proper URL to post to


# build the form to send. the server is quite anal about what it will take.
my $script_response = $agent->post(
	$req->uri,
	[
		'gi_number' => "NP_009059",
		'sequences_to_select' => "BEST",
		'seq_identity_filter' => "90"

	]
);
$positions_found= 0;
$matrix_found = 0;
$pred_found = 0;
open (OUT,">$tmp/$$.out") || die ("Cannot open file");

$content = $script_response->content;
chomp($content);

print OUT $content;


