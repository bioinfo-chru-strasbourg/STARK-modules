#!/usr/local/bin/perl
$| = 1;
use Getopt::Std;
use File::Basename;

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software
#

#       Set file permissions to rw-rw----
system("umask 006");

$ENV{'SIFT_HOME'} = '/usr/local/projects/SIFT/sift4.0/';
my $SIFT_HOME = $ENV{'SIFT_HOME'};
use vars qw ($opt_i $opt_c $opt_o);
getopts ("i:c:o:");


my $usage = "usage:
$0
	-i <List of variants>
	-c <File with coding coordinates in gff format>

	<OPTIONAL>
        -o <output file with complete path: default=$SIFT_HOME/tmp>
";

if (! defined($opt_i)) { 
        print STDERR $usage;
        die "No SNV list entered";
}


if (! defined($opt_c)){
        print STDERR $usage;
        die "No coding file was entered";
}


my $bin             = "$SIFT_HOME/bin";
my $tmp             = "$SIFT_HOME/tmp";
my $all_chr_file =  $opt_i;
my ($filename, $dir, $ext) = fileparse  ($all_chr_file, '\.*'); 
my $outfilename = defined ($opt_o) ? "$opt_o" :  $tmp . "/" . $filename . "_codingOnly";

my $coding_file = $opt_c; 
#$coding_info_dir . "/Homo_sapien/ens.hum.ncbi36.ver41.cds.merge.gff";

#Read input list of chromosome coordinates and add to all_chr_string
my $all_chr_string;

open (CHR_FILE,"$all_chr_file") || die "Cannot open $all_chr_file";
while (<CHR_FILE>){
        if ($_ =~ /\d/ && $_ =~ /\,/){
                $first_line = $_;
                last;
        }
}
close(CHR_FILE);

if ($first_line =~ /\d+,\d+,\d+,\-?1,/i){
        $COORD_SYSTEM = "SPACE";
}
elsif($first_line =~ /\d+,\d+,\-?1,/i){
        $COORD_SYSTEM = "RESIDUE";
}
else{
        print "Incorrect input format. Please see http://sift.jcvi.org/www\/chr_coords_example.html";
        last;
}

print
"Your input data has been recognized to use $COORD_SYSTEM based coordinate system\n";
# if COORD SYSTEM is SPACE, then there are insertions with same coordinates

#if COORD SYSTEM is ABSOLUTE then convert chr file to space based.
my $new_chr_file = "$SIFT_HOME/tmp/" . $filename .  ".gff";
if ($COORD_SYSTEM eq "RESIDUE"){
	convert_residue_to_space_gff ($all_chr_file, $new_chr_file);
}
if ($COORD_SYSTEM eq "SPACE") {
	convert_space_with_insertions_gff ($all_chr_file, $new_chr_file);
}
	

 system ("sort -k1,1 -k4,4n -k5,5n $new_chr_file > $new_chr_file.sorted");
my @lines_to_get = `$bin/IntersectLocations.sh $new_chr_file.sorted gff $coding_file gff simple | grep -v Total | cut -f1 | uniq`;

my $all_chr_in_cds_file = $all_chr_file . ".cds.txt";

my ($cds_count, $non_cds_count) = get_lines ($new_chr_file, 
			$outfilename, 
			@lines_to_get);

print "Your filtered output is in $outfilename\n";

exit(0);

sub
get_lines
{
	my  ($all_chr_file, $all_chr_in_cds_file, @lines_to_get)  = @_;

	my %line_to_get_hash ;
	# more elegant way to do this by iterating through @lines_to_get, but oh well
	for (my $i=0; $i < @lines_to_get; $i++) {
		chomp ($lines_to_get[$i]);
		$lines_to_get_hash{$lines_to_get[$i]} = 1;
	}
	open (CHR_FILE, $all_chr_file) || die "can't open $all_chr_file";
	open (OUT_CHR_FILE, ">$all_chr_in_cds_file") || die "can't open $all_chr_in_cds_file";
	my $in_cds = 0;
	my $not_in_cds = 0;
        my $line;
	while ($line = <CHR_FILE>){
		chomp ($line);
		my @fields = split (/\t/, $line);
		my $line_num = $fields[1];
		if (exists ($lines_to_get_hash{$line_num})) {
			print OUT_CHR_FILE $fields[8] . "\n";
			$in_cds++;
		} else { 
			$not_in_cds++;
		}
	}
	close (CHR_FILE);
	close (OUT_CHR_FILE);
	return ($in_cds, $not_in_cds);
}

sub
convert_space_with_insertions_gff
{
       my ($all_chr_file, $new_chr_file) =  @_;

        open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
        open (CHR_FILE_NEW,">$new_chr_file") || die ("Cannot open $new_chr_file");
        my $line_num = 0;
        my $line;
	while ($line = <CHR_FILE>){
                chomp ($line);
                if ($line !~ /\d/){next}
                @elts = split (/\,/, $line);
                $chr = @elts[0];
                $coord1 = @elts[1];
                $coord2 = @elts[2]; 
                $orn = @elts[3];
                $alleles = @elts[4];
                $comment = @elts[5];
		# for insertions
		if ($coord1 == $coord2) {
			$coord1 = $coord1 - 1;
		}
                print CHR_FILE_NEW "$chr\t$line_num\t$line_num\t$coord1\t$coord2\t.\t+\t.\t$line\n"; 
                $line_num++;
        }
        close(CHR_FILE);
        close (CHR_FILE_NEW);

} # end convert_space_with_insertions_gff

sub
convert_residue_to_space_gff
{
	my ($all_chr_file, $new_chr_file) =  @_;

   	open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
        open (CHR_FILE_NEW,">$new_chr_file") || die ("Cannot open $new_chr_file");
	my $line_num = 0;
        my $line;
	while ($line = <CHR_FILE>){
                chomp ($line);
                if ($line !~ /\d/){next}
                @elts = split (/\,/, $line);
                $chr = @elts[0];
                $coord2 = @elts[1];
                $coord1 = $coord2-1;
                $orn = @elts[2];
                $alleles = @elts[3];
                $comment = @elts[4];
                print CHR_FILE_NEW "$chr\t$line_num\t$line_num\t$coord1\t$coord2\t.\t+\t.\t$line\n"; 
		$line_num++;
        }
        close(CHR_FILE);
        close (CHR_FILE_NEW);

}



# returns hash for a file, 2nd field is the key and the 3rd field
# is the value 4th field, is the delimiter
sub make_hash {
	my ($file) = @_;
	my %hash;
	open( HASH, $file ) || die "can't open $file";
	my $line;
	while ( $line = <HASH> ) {
		chomp($line);
		if ( exists( $hash{$line} ) ) {
			$hash{$line}++;
		}
		else {
			$hash{$line} = 1;
		}
	}
	close(HASH);
	return (%hash);
}


sub round {
    my($number) = shift;
    return int($number + .5);
}

