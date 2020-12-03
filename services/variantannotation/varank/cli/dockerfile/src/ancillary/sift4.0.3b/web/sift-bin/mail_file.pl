#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use MIME::Lite;
use Net::SMTP;

my $from_address = "sift";
my $to_address = $ARGV[0];
my $file = $ARGV[1];
my $subject_header = $ARGV[2];

my $msg = MIME::Lite->new (
	From => $from_address,
	To => $to_address,
	Subject => $subject_header,
	Type => "multipart/mixed"
); 

$msg->attach (	
	Type => 'TEXT/HTML',
	Path => "$file",
	Filename => "$file",
	Disposition => 'attachment'
)
;

$msg->send;

