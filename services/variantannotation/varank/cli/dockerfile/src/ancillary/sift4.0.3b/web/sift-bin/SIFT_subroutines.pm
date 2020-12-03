# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

sub mail_results {
	my ($address, $return_address, $tmp, $pid);

	# send alignment
        $subject = "SIFT alignment";
        $tempstring = "\"$subject\"";
        $subject = $tempstring;
        $file = $tmp . "/$pid" . ".alignedfasta";
        system ("mailx -s $subject -r $return_address $address < $file");
if ($status ==0) {
	print "Mailing alignment";
} else {
	print"Failed to mail alignment";
}
# send probabilities
        $subject = "SIFT conditional probabilities"; $tempstring ="\"$subject\""
;
        $subject = $tempstring;
        $file = $tmp . "/$pid" . ".siftresults.matrix";
        system ("mailx -s $subject -r $return_address $address < $file");

#send predictions
        @table_files = `ls $tmp/$pid.aatable*`;
        $index = 0;
        for $file (@table_files) {
# SUBJECT HEADING
# subject heading MUST be in the following order, don't try to simplify it else
# Perl may not assign what you would expect it to.
	        $subject = "SIFT positions ";
       		$subject .=  $index *100 + 1 . " to " . ($index +1) * 100;
	        $tempstring = "\"$subject\"";
       		$subject = $tempstring;

                $index++;
	        system("mailx -s $subject -r $return_address $address < $file");
        }
# mail any errors
	@error_files = `ls $tmp/$pid.*.error`;
	for $file (@error_files) {
		$subject = "SIFT error";
		$tempstring = "\"$subject\""; $subject = $tempstring;
		system ("mailx -s $subject -r $return_address $address < $file");
	}

}

sub print_psiblast_alignment_desc {

	print "The alignment taken from PSIBLAST is returned in msf format.<BR>";
	print "<b>Note:</b>  <b>X</b>es are placeholders at the beginning and end of";
	print " sequences.  While <b>-</b> means a gap in the alignment ";
	print "an X means a lack of information such as a partial alignment or ";
	print "incomplete sequence and do not contribute to the prediction.  <BR>";
} # end of print_psiblast_alignment_desc

sub mail_tolerated_nottolerated_aminoacidtable {

	my ($tmp,$address, @table_files) = @_;

	my $subject = "SIFT: Prediction Tables";
#	my $counter = 0;
	my $mailer=Mail::Mailer->new();
	$mailer->open ( {From => $SIFT_ADDRESS,
       	          To   => $address,
       	          Subject => $subject,
       	         }) or die "Can't open : $!\n";

	foreach $table_file (@table_files) {
	        open (TABLEFILE, "<$table_file");
      	 	 while ($_ = <TABLEFILE>) { 
			print $mailer "$_";  
		}
       		 close (TABLEFILE);
	}

	close ($mailer);
}

sub tolerated_nottolerated_aminoacidtable {

	my ($tmp,$no_of_tablefiles, $pid) = @_;

	my ($i, $table_file);

	print "<BR>";
	print ("S<font color=red>I</font>F<font color=blue>T</font> amino acid predictions for: <BR>");
	for ($i = 1; $i <= $no_of_tablefiles; $i++) {
		$table_file = $pid . ".aatable" . $i;
		print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$table_file\" TARGET=_blank>";
		printf "Positions %d to %d\n", ($i-1)*100 + 1 , $i *100;
		print "</A><BR>\n";
;
	}

}

sub print_logos_desc {
#
# took out 10/21/00 because took too long and wasn't looking very pretty
# anyways
# if want to put it back in, put this in SIFTING.csh file
####  START of SIFTING.csh
        # do logos stuff
#echo doing logo stuff
#       /home/sift/bin/logos.csh $pid $alignedseqs $polymorphism_file
# converts postscript to gif and pdf files.  can't combine this script with the
# previous because directory is changed in previous. having problems calling
# convert.  whatever.
#       /home/sift/bin/logo_nongrouped_web.csh $pid
# END of SIFTING.csh 

# and then in perlscript which writes output for web call this file by:
#
# START
#$no_pic_files = `ls $tmp/$$.*.ps | wc -l`;
#print_logos_desc($no_pic_files, $$);
# END 
	my ($no_of_files, $pid) = @_;
	my ($i);

	print "<b>LOGOS</b><BR>\n";
	print "If you have submitted amino acid substitutions, then these are ";
	print "indicated in the LOGO as stars.  <BR>";
	print " <font color=red>Red</font> stars indicate ";
	print "<font color=red>intolerant</font> substitutions; ";
	print "<b>black</font> stars";
	print " indicate substitution predicted to be <b>tolerated</b>.<BR>";
#	print "<BR>Postscript images also show what amino acids are allowed at each position";
#	print " after adding pseudocounts.  Amino acids that appear in the alignment are in solid color; those that have been allowed after adding pseudocounts are in outline.  ";
#	print "Outlined amino acids do not contribute to information.<BR>";
	for ($i = 1; $i <= $no_pic_files; $i++) {
 	       printf "Positions %d to %d\n", ($i-1)*50 +1, $i*50;
       	 	$pic_file = "$tmp" . $$ . "." . $i . ".ps" ;
		print "<A HREF=\"/sift-bin/display_img.csh?$pic_file+ps\">Postscript</A>\n";
		print " <A HREF=\"/sift-bin/display_img.csh?$pic_file.pdf+pdf\">PDF</A>\n";
		print " <A HREF=\"/sift-bin/display_img.csh?$pic_file.gif+gif\">GIF</A><BR>\n";
	}

}


sub
print_score_desc {
	my ($tmp,$file) = @_;

	print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$file.matrix+Probabilities\" TARGET=_blank>Scaled Probabilities for Entire Protein</A><BR>\n";
print "<i><b>May take some time to load!!</b> Please be patient if you do not see the table immediately.</i><BR>";
#	print "Scores are scaled probabilities.  \n";
	print "Amino acids with probabilities < .05 are predicted to be <font color=red>deleterious</font>.<BR>";

#, and are ";
#	print "highlighted in <font color=red>red</font>.<BR>";
}

sub 
print_predictions {
	my ($tmp,$file) = @_;

	print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$file.predictions+Predictions+PRE\" TARGET=_blank>Predictions of substitutions entered</A><BR>\n";
}

sub trim {
	my @out = @_;
	for (@out) {
		s/^\s+//;
		s/\s+$//;
	}
	return wantarray ? @out : $out[0];
}

1;
