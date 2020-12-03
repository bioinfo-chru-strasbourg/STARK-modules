#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$tmp = "/opt/www/sift/tmp/";
$pid = @ARGV[0];

$address = @ARGV[1];

$num_tables = `ls $tmp/$pid.aatable* | grep -v html | wc -l`;

if ( $address ne "" ) {
	open( MESSAGE, ">$tmp/$pid.email_message.txt" );
	print MESSAGE
"Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent query attached with this message.\nRemember this job id \"$pid\" for any future correspondance.\nHere is the description of the attached files:\n1. $pid.alignedfasta: Multiple alignment of your protein sequence with its homologs.\n2. $pid.siftresults.matrix.html: Conditional Probability Matrix.\n3. $pid.siftresults.predictions.html: Predictions on substitutions entered by you. If you did not submit any, this file should be empty.\n4. $pid.aatable(n).html: Predictions for all positions in the sequence submitted by you.\n\nPlease do not hesitate to contact us if you have any questions about SIFT.\n\nThanks\nSIFT Team\nJ Craig Venter Institute (West Coast Campus)\n10355 Science Center Drive\nSan Diego, CA 92121\nUSA";
	close(MESSAGE);

	my $mutt_command =
"mutt -F /opt/www/sift/htdocs/.muttrc -a $tmp/$pid.alignedfasta -a $tmp/$pid.siftresults.matrix.html -a $tmp/$pid.siftresults.predictions.html";
	system(
		"zip $tmp/$pid.aatable.html.zip $tmp/$pid.aatable*.html"
	);
	$mutt_command .=
" -a $tmp/$pid.aatable.html.zip -s \"SIFT Results for Job ID $pid\" $address <$tmp/$pid.email_message.txt";
	system("$mutt_command");

}
