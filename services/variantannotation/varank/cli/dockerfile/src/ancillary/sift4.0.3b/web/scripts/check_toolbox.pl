#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
($DAY, $MONTH, $YEAR) = (localtime)[3,4,5];
$y = $YEAR + 1900;
$m = $MONTH + 1;
$d = $DAY;
$datestamp = "$y-$m-$d";
$tests_completed = 0;
$failed = 0;

mkdir "/opt/www/sift/htdocs/daily_log/$datestamp";
system("perl /opt/www/sift/htdocs/scripts/check_SIFT_dbSNP.pl");
system("perl /opt/www/sift/htdocs/scripts/check_SIFT_pid_subst_all.pl");
system("perl /opt/www/sift/htdocs/scripts/check_SIFT_BLink.pl");
system("perl /opt/www/sift/htdocs/scripts/check_SIFT_Seq.pl");
system("perl /opt/www/sift/htdocs/scripts/check_SIFT_Genome.pl");
system("perl /opt/www/sift/htdocs/scripts/check_SIFT_indels.pl");
sleep 600; # Pauline added, give 10 minutes for it to finish running
@test_files = </opt/www/sift/htdocs/daily_log/$datestamp/*>;
$num_files = $#test_files;
#print "num _files $num_files\n";
if ($num_files >= 3){
	$tests_completed = 1;
} 
foreach $file (@test_files){
	if ($file =~ /Failed/i){
		$failed = 1;
		last;
	}
}
#print "failed status $failed\n";
if ($failed == 1 ||  $tests_completed == 0){
	system("mutt -s \"SIFT Failure Notice!\"  smurphy\@jcvi.org png\@jcvi.org < /opt/www/sift/htdocs/scripts/failure_message.txt");
	#system("mutt -s \"SIFT Failure Notice!\"  pkumar\@jcvi.org < /opt/www/sift/htdocs/scripts/failure_message.txt");
}
