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
use LWP::Simple;
use URI::URL;
($DAY, $MONTH, $YEAR) = (localtime)[3,4,5];
$y = $YEAR + 1900;
$m = $MONTH + 1;
$d = $DAY;
$datestamp = "$y-$m-$d";
$tmp = "/opt/www/sift/tmp/";
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
  $agent->get('http://sift.jcvi.org/sift-bin/SIFT_feed_to_chr_coords.pl');

if ( !$pre_response->is_success ) {
	die $pre_response->status_line
	  ;    # if it fails at this point, it's likely a network or server issue.
}

my $req = $pre_response->request; # needed for getting the proper URL to post to


# build the form to send. the server is quite anal about what it will take.
my $script_response = $agent->post(
	$req->uri,
	[
		'organism' => "Human_db",
		'CHR' => "10,115912482,-1,C\/T\n10,115900918,-1,G\/T\n16,69875502,-1,G\/T",
		'oo1' => "1",
		'oo2' => "1",
		'oo3' => "1",
		'oo4' => "1",
                'oo5' => "1",
                'oo6' => "1",
		'oo7' => "1",
                'oo8' => "1",
                'oo9' => "1",
		'oo10' => "1",
                'oo11' => "1",
                'oo12' => "1"


	]
);
$job_id = "";
$output_exists = 1;
$num_output_fields=0;
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
	#print "$line\n";
	if ( $line eq "" ) {
		next;
	}
	
	if ($line =~ /Your job id is (\d+) and/){
		$job_id = $1;
	}
}
$status = getstore("http://sift.jcvi.org/tmp/$job_id\_predictions.html","/opt/www/sift/tmp/$job_id\_predictions.html");
$output = "$tmp/$job_id\_predictions.html";
$output_content = "";
if (-e "$output"){
	$output_exists = 1;
	open (OUTPUT, "$output") || die ("Cannot open output");
	while (<OUTPUT>){
		chomp;
		@elts = split /<td>/,$_;
		$num_output_fields+= scalar @elts;
		$output_content.="$_\n";
	}
}
else{
	$output_exists = 0;
}	

print "$job_id\t$output_exists\t$num_output_fields\n";
if ($job_id eq "" || $output_exists ==0 || $num_output_fields != 84){
	print "SIFT_Genome Failed!\n";	
	open ( OUTFILE, ">/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_Genome_Failed.html");
	print OUTFILE $output_content;
	close (OUTFILE);
}
else{
	print "SIFT_Genome Passed!\n";
	open ( OUTFILE, ">/opt/www/sift/htdocs/daily_log/$datestamp/SIFT_Genome_Passed.html");
	print OUTFILE $output_content;
	close(OUTFILE);
}




