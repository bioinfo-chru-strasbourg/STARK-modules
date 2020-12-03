#!/usr/local/bin/perl
use List::Util qw[min max];
use File::Copy;
use Getopt::Std;


# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

$| = 1;
system("umask 006");
$ENV{'SIFT_HOME'} = '/usr/local/projects/SIFT/sift4.0.2/';
my $SIFT_HOME = $ENV{'SIFT_HOME'};
use vars qw($opt_i $opt_d $opt_o $opt_A $opt_B $opt_C $opt_D $opt_E $opt_F $opt_G $opt_H $opt_I $opt_J $opt_K $opt_L);
getopts("i:d:o:A:B:C:D:E:F:G:H:I:J:K:L:");
my $usage = "usage: 
$0 
        -i <Query SNP filename with complete path>
        -d <Variation db directory path>

	<OPTIONAL>
        -o <output file with complete path: default=$SIFT_HOME/tmp>
	-A 1 to output Ensembl Gene ID: default: 0
	-B 1 to output Gene Name: default: 0
	-C 1 to output Gene Description: default: 0
	-D 1 to output Ensembl Protein Family ID: default: 0
	-E 1 to output Ensembl Protein Family Description: default: 0
	-F 1 to output Ensembl Transcript Status (Known / Novel): default: 0
	-G 1 to output Protein Family Size: default: 0
	-H 1 to output Ka/Ks (Human-mouse): default: 0
	-I 1 to output Ka/Ks (Human-macaque): default: 0
	-J 1 to output OMIM Disease: default: 0
	-K 1 to output Allele Frequencies (All Hapmap Populations - weighted average): default: 0
	-L 1 to output Allele Frequencies (CEU Hapmap population): default: 0
";
if(!(
        defined($opt_i) 
        )){
        print STDERR $usage;
        die "No SNV list entered";
}
if(!(
        defined($opt_d))){
        print STDERR $usage;
        die "No database file entered";
}

my $oo1 = defined($opt_A) ? 1:0;
my $oo2 = defined($opt_B) ? 1:0;
my $oo3 = defined($opt_C) ? 1:0;
my $oo4 = defined($opt_D) ? 1:0;
my $oo5 = defined($opt_E) ? 1:0;
my $oo6 = defined($opt_F) ? 1:0;
my $oo7 = defined($opt_G) ? 1:0;
my $oo8 = defined($opt_H) ? 1:0;
my $oo9 = defined($opt_I) ? 1:0;
my $oo10 = defined($opt_J) ? 1:0;
my $oo11 = defined($opt_K) ? 1:0;
my $oo12 = defined($opt_L) ? 1:0;
my $pid             = $$;
my $bin             = "$SIFT_HOME/bin";
my $tmp = defined($opt_o) ? "$opt_o/$pid" : "$SIFT_HOME/tmp/$pid";
my $Variation_db_dir = $opt_d;
mkdir $tmp;
my $num_coords_per_split = 1000;
require "$bin/SIFT_subroutines.pm";
my $input_chr_file = $opt_i;

my $output_options = "$oo1,$oo2,$oo3,$oo4,$oo5,$oo6,$oo7,$oo8,$oo9,$oo10,$oo11,$oo12";
my $all_chr_file = "$tmp/$pid.allchrfile";
copy ($input_chr_file , $all_chr_file);
open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file for validation");
while (<CHR_FILE>){
	if ($_ =~ /\d/ && $_ =~ /\,/){
		$first_line = $_;
		last;
	}
}
close(CHR_FILE);
if ($first_line =~ /[\d+,X,x,Y,y],\d+,\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i){
	$COORD_SYSTEM = "SPACE";
}
elsif($first_line =~ /[\d+,x,X,y,Y],\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i){
	$COORD_SYSTEM = "RESIDUE";
}
else{
	print "Error: Incorrect input format. Please see <A HREF=\"\/www\/chr_coords_example.html\">Sample format</A>";
	last;	
}

print
"Your input data has been recognized to use $COORD_SYSTEM based coordinate system. Your job id is $pid and is currently running.  Your job has been partitioned into datasets of $num_coords_per_split positions and the status of each job can be viewed $tmp/$pid.outpage.txt. Once the status of a job is 'Complete', you may view the results. A partitioned job with $num_coords_per_split input rows typically takes 6-7 min to complete.\n";


#prepare output status page
open (OUTPAGE,">>$tmp/$pid.outpage.txt") || die ("cannot open outpage");
print OUTPAGE "SIFT Results Status\n\n";
$heading = "Job\tJob size\tJob ID\tJob status\tResults file\n";


#if COORD SYSTEM is ABSOLUTE then convert chr file to space based.
if ($COORD_SYSTEM eq "RESIDUE"){
	open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
	open (CHR_FILE_NEW,">$all_chr_file.new") || die ("Cannot open new all chr file");
	while (<CHR_FILE>){
		chomp;
		if ($_ !~ /\d/){next}
		@elts = split /\,/, $_;
		$chr = @elts[0];
		$coord2 = @elts[1];
		$coord1 = $coord2-1;
		$orn = @elts[2];
		$alleles = @elts[3];
		$comment = @elts[4];
		print CHR_FILE_NEW "$chr,$coord1,$coord2,$orn,$alleles,$comment\n";
	}
	close(CHR_FILE);
	close (CHR_FILE_NEW);
	system("mv $all_chr_file.new $all_chr_file");
}


open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
$count = 0;
$new_pid = $pid+1;
$input_set = 1;
open (OUTFILE,">$tmp/$new_pid.chrfile");

while (<CHR_FILE>){
	chomp;
	if ($_ !~ /\d/){next}
	$count++;
	if ($count % $num_coords_per_split== 0){
		push @pid_list,$new_pid;
		print OUTFILE "$_\n";
                $job_size = $num_coords_per_split;
		$start = $count-$job_size+1;
		$end = $count ;
                print OUTPAGE "Partitioned set $input_set\t";
                print OUTPAGE "Input rows $start to $end\t";
                print OUTPAGE "$new_pid\t";
                print OUTPAGE "Not started\tNot available\n";
		close(OUTFILE);
		$input_set++;
                $new_pid++;

		open (OUTFILE,">$tmp/$new_pid.chrfile");
	}
	else{
		print OUTFILE "$_\n";
		
	}
}
close(CHR_FILE);
close(OUTFILE);
$job_size = $count % $num_coords_per_split;
$start = $end+1;
$end = $start + $job_size-1;
if ($job_size == 0){
	system("rm -f $tmp/$new_pid.chrfile");
	print OUTPAGE "Complete set\t";
        print OUTPAGE "Input rows 1 to $end\t";
        print OUTPAGE "$pid\t";
        print OUTPAGE "Not started\tNot available\n";

}
else{
	push @pid_list, $new_pid;
	print OUTPAGE "Partitioned set $input_set\t";
	print OUTPAGE "Input rows $start to $end\t";
	print OUTPAGE "$new_pid\t";
	print OUTPAGE "Not started\tNot available\n";
	print OUTPAGE "Complete set\t";
        print OUTPAGE "Input rows 1 to $end\t";
        print OUTPAGE "$pid\t";
        print OUTPAGE "Not started\tNot available\n";

}
print OUTPAGE "\n\n";
print OUTPAGE "Batch Report\n";
print OUTPAGE "Number of input (non-intronic) variants: \n";
print OUTPAGE "Coding variants: \n";
print OUTPAGE "Coding variants predicted: \n";
print OUTPAGE "Tolerated: \n";
print OUTPAGE "Damaging: \n";
print OUTPAGE "Nonsynonymous: \n";
print OUTPAGE "Synonymous: \n";
print OUTPAGE "Novel: \n";

close (OUTPAGE);
#system("cp $tmp/$pid.outpage.txt $tmp/$pid.outpage.swap.txt");
#direct to the outpage url
#print "Location: $outpage_url\n\n";
#print_outpage();

open (BATCH_FILE,">$tmp/$pid.batchfile");
for ($i = 0; $i < scalar @pid_list; $i++){
	$new_pid = @pid_list[$i];
	chomp $new_pid;
	if ($i == scalar @pid_list -1){
		print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tLAST\t$output_options\t$address\n";
	}
	else{
		print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tNOT_LAST\t$output_options\t$address\n";
	}
}
system("$bin/sift_feed_to_chr_coords_batch.pl $tmp/$pid.batchfile $tmp $Variation_db_dir");
close (BATCH_FILE);













sub print_outpage{
	open(OUTPAGE, "$tmp/$pid.outpage.txt") || die("cannot open outpage");
	while (<OUTPAGE>){
        	print;
	}
	close (OUTPAGE);
}








sub general_parse {
        local ( $content_type, $query_string ) = @_;
        local ( %ans, @q, $pair, $loc, $boundary, $temp, $temp1 );

        if ( $content_type eq "application/x-www-form-urlencoded" ) {

                # break up into individual name/value lines
                @q = split( /&/, $query_string );

                foreach $pair (@q) {

                        # break the name/value pairs up
                        # use split rather than regular expressions because the value
                        # may have
                        #  newlines in it
                        split( /=/, $pair, 2 );

                        # change '+' to ' '
                        $_[1] =~ s/\+/ /g;

                        # change the escaped characters (must be after the split on '&'
                        # and '=')
                        $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

                        $ans{ $_[0] } = $_[1];
                }    #end of foreach $pair

        }    #end of if ($content_type)
        else {
                $loc = index( $content_type, "boundary=" );
                if ( $loc > 0 ) {
                        $temp = substr( $content_type, $loc + 9 );

                 #               Why is this necessary? (boundary= doesn't match actual)
                        $boundary = "--" . $temp;

                        # break up into individual name/value lines
                        @q = split( /$boundary/, $query_string );

                        foreach $pair (@q) {

                                # break the name/value pairs up
                                $loc = index( $pair, "name=" );
                                $temp = substr( $pair, $loc + 5 );

                                #         $loc = index($temp, "\n\n");
                                $loc = index( $temp, "\n" );
                                $temp1 = substr( $temp, $loc + 2 );

                                #   Get rid of stuff after the name; including semicolon if any
                                $loc_semi = index( $temp, ";" );
                                $loc_eol  = index( $temp, "\n" );
                                $loc      = $loc_eol;
                                if ( $loc_semi > 0 && $loc_semi < $loc ) {
                                        $loc = $loc_semi;
                                }
                                if ( $loc > 0 ) { $temp = substr( $temp, 0, $loc ); }

                                #               Get rid of quotes around the name
                                $temp =~ s/\"//g;

                                #               Still has a trailing whitespace character ...
                                $temp =~ s/\s//g;

                  #               Substitute spaces with nothing
                  #               Need to strip leading/ending whitespace off of $temp1,
                  #               but be careful not to strip off internal CRs
                  #               MAC file lines end in just \r, no \n, so makelis won't
                  # find all
                  #               of the sequences; DOS file lines end in \r\n, UNIX in
                  #\n.
                  #               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
                 #In my case, I want to keep the newlines in fields which have "file" or
                 # 'seq"
                 # and remove newlines everywhere else.
                 #if ( $temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/ ) {
                                $temp1 =~ s/\r\n/\n/g;
                                $temp1 =~ s/\r/\n/g;

                                #}

                         # for other variables that are NOT files or file-like, remove extra
                         #whitespace
                         #else { $temp1 =~ s/\s//g; }
                                if ( $temp ne "" ) { $ans{$temp} = $temp1; }
                        }    # end of foreach
                }    #end of if loc > 0
                else {
                        print "Cannot parse\n";
                        print "content_type=$content_type\n";
                        print "query_string=$query_string\n";
                }
        }
        return %ans;

        #print "</PRE>";
}    # end of general_parse
sub hextodec {
        unpack( "N", pack( "H8", substr( "0" x 8 . shift, -8 ) ) );
}

