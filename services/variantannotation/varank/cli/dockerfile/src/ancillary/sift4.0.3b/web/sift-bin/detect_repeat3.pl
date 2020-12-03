#!/usr/local/bin/perl 
### Program to print the effect of indels in a coding repeat region on the protein domain###
### Input: File containing all the input indels with the following information:
###        Chromosome number, ENST, indel_aa_locations from the denormal file, the location of original protein sequence, location of modified protein sequence
### Output: Produces the normal repeat region amino acid sequence and the modified amino acid repeat region ###

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

use strict;
#use warnings;
#do "test_aa_count.pl";

$" = "\t";
our $infile1 = $ARGV[0]; ### Takes in the indel input file as a command line arg
our $lines1; 

open my $FH,"<$infile1" or die "File1 not found";

our $FH1; #File handel for the repeat database file
our $FH2; #File handel for the normal protein file
our $FH3; #File handel for the modified protein sequence file
our @ref = ();#stores the reference repeat region
our @modi = ();#stores the modified repeat region

### Reading each linr in the indel file ###
while(my $lines1 = <$FH>){
    my @f1; 
    @f1 = get_indel($lines1);
    #print $lines1,"\n";
#    our $infile = $f1[0];
    our $infile = "/opt/www/sift/db/Human_db/repeat_db.gff"; ## stores the repeat file
    open $FH1,"<$infile" or die "Repeats data file not found";

### Condition to check for the presence of amino acid coordinate from the denormal file 
### Absence of aa cords in column 5 and 6 of input file implies indel not found in the CDS

#    if(($f1[4] == "\t")&&($f1[5]=="\t")){
#	$f1[0] =~ s/chr//g;
#	my @exe = ($f1[1]);
#	$f1[1] = join(//,@exe);
#	print "@f1\t\n";
#    }
#    else{
compare(@f1);
#    }
}


### Sub to extract the repeat coordinates in protein domain for each indel from the Repeat database file ###
### Arg: @in - stores the elements of each row from the input file
sub compare{
    my (@in) = @_;
    my @f2;   
    while(my $lines2 = <$FH1>){
     	if($lines2 !~ /enst/) { #&&($lines1 !~ /.gff/)&&($lines1 !~ /^\t/)){
	    @f2 = get_aafile($lines2); ## Stores the elements of each row from the repeat database file
	}
### Condition to check if the indel lies within the coding repeat region
### Calls the sub: get_seq - to extract the repeat regions from the normal and modified protein seq
###                short_repeat - to print the modified repeat region
  
	if($in[0] eq $f2[0]){
	    if($in[1] eq $f2[5]){
		if(($in[4] >= $f2[9]) && ($in[5]<=$f2[10]) && ($in[2]>=$f2[7]) && ($in[3] <= $f2[8])){
		    push(@in,$f2[9]);
		    push(@in,$f2[10]);
		    get_seq(@in);
		    my $reference = join('',@ref);
		    my $modified = join('', @modi);
		    my @exe = ($in[1]);
		    $in[1] = join(//,@exe);
		    print "@in\t";
		    short_repeat($reference, $modified);
		    @ref =();
		    @modi = ();
		    return;
		}
		#else{
		#    chomp($in[-1]);
		#    $in[0] =~ s/chr//g;
		#    my @exe = ($in[1]);
		#    $in[1] = join(//,@exe);
		#    print "@in\t\n";
		#}
	    }
	}
    }	
    close($FH1);
}

### Sub to obtain the elements of each row from the input indel file and to store them in an array
 
sub get_indel{

    my ($lines) = @_;

    my @ele = split(/\t/,$lines);
 #   my @add_chr = ("chr",$ele[0]);
#    my $chr = join('',@add_chr);
#    $arr[0] = $chr; 
    my @arr;

    $arr[0] = $ele[0];  
    $arr[1] = $ele[1];
    #$arr[1] =~ s/ENST//ig;
    $arr[2] = $ele[2];
    $arr[3] = $ele[3];
    $arr[4] = $ele[4];
    $arr[5] = $ele[5];
    $arr[6] = $ele[6];
    $arr[7] = $ele[7];
    chomp($arr[7]);
    return @arr;
}

### Sub to obtain the elements of each row from the repeat database file and to store them in an array
sub get_aafile{
    my($lines)=@_;
    my @ele = split(/\t/,$lines);
    return @ele;
}


### Sub to obtain the repeat region from the normal protein sequence and the modified protein sequence
### Calls subs: delete_seq to delete the amino acids flanking the repeat region from the normal protein file
###             and new_prot_del_seq to delete teh amino acids flanking the repeat regions from the modified protein file
### Argument: array with the elements form a row in indel file and the corresponding amino acid repeat coordinates 
sub get_seq{
    my @f1 = (@_);
    my @del;
    my $infile2 = $f1[6];
    my $infile3 = $f1[7];
    open $FH2,"<$infile2" or die "Protein sequence file not found";
    @del = delete_seq(@f1);
    close($FH2);
    open $FH3,"<$infile3" or die "Modified protein sequence file not found";
    new_prot_del_seq(@del);
    close($FH3);
}

### Sub to extract the repeat region from the normal protein sequence
### Calls the sub parse to count the number of repeats in the repeat region 
### Returns the start flanking sequence and the end flanking seqence of the repeat region 
sub delete_seq{
    my (@in) =@_;
    my $holdTerminator = $/;
    undef $/;
    my $buf = <$FH2>;
    $/ = $holdTerminator;
    my @lines = split /$holdTerminator/, $buf;
    $buf = "init";
    $buf = join $holdTerminator, @lines;
    $buf =~ s/>ENST\d//g;
    $buf =~ s/\d//g;
    $buf =~ s/\n//g;

    my @a;
    $a[1] = substr($buf, $in[9]+1, 100000,"");
    $a[0] = substr($buf, 0, $in[8]-3,"");
    chomp($in[-1]);
    $in[0] =~ s/chr//g;
    my @exe = ("ENST",$in[2]);
#    $in[2] = join(//,@exe);
    @ref = parse($buf);
    return @a;
}


### Sub to count the number of mono, di and tri amino acids repeats within the repeat region
### Argument: A string with the repeat sequence
### Returns the parsed repeat sequence as an array 
sub parse{

    my($repeat) = @_;
    my $count = 1;
    my $position = 0;
    my @amino = split(//,$repeat);
    my $flag = 0;
    my @output;
    for(my $i=0 ; $i<scalar(@amino); ++$i){
	if($i+1 != scalar(@amino)){
            if($amino[$i] eq $amino[$i+1]){
                while($amino[$i+1] eq $amino[$i]){
		    $count =$count+1;
		    $i =$i+1;
                }
		my $aa = lc($amino[$i]);
		push(@output, "($aa)$count");
		#print @output;
		$count =1;
            }

            elsif(($amino[$i] eq $amino[$i+2])&&($amino[$i+1] eq $amino[$i+3])){
                while(($amino[$i] eq $amino[$i+2])&&($amino[$i+1] eq $amino[$i+3])){
                    $count = $count + 1;
                    $i = $i+2;
                }
		my $aa1 = lc($amino[$i]);
		my $aa2 = lc($amino[$i+1]);
               push(@output, "($aa1$aa2)$count");
	#	print $output;
                $count = 1;
                $i =$i+1;
            }

            elsif(($amino[$i] eq $amino[$i+3])&&($amino[$i+1] eq $amino[$i+4])&&($amino[$i+2] eq $amino[$i+5])){
                while(($amino[$i] eq $amino[$i+3])&&($amino[$i+1] eq $amino[$i+4])&&($amino[$i+2] eq $amino[$i+5])){
                    $count = $count + 1;
                    $i = $i+3;
                }
		my $aa1 = lc($amino[$i]);
		my $aa2 = lc($amino[$i+1]);
		my $aa3 = lc($amino[$i+2]);
                push(@output,"($aa1$aa2$aa3)$count");
		#print $output;
		$count = 1;
                $i = $i+2;
            }
            else{
		push(@output,"$amino[$i]");
                #print "$amino[$i]";
            }
        }
        else{
            if($amino[$i] ne $amino[$i-1]){
		push(@output,"$amino[$i]");
                #print "$amino[$i]";
            }
        }
    }

    return @output;
}


### Sub to obrain the repeat region from the modified protein sequence
### Arg: Array witht he start flankinf and end flanking sequence
### Calls sub parse to obtain the mono, di and tri aa repeats within the repeat region

sub new_prot_del_seq{

    my (@match) = @_;
    my $holdTerminator = $/;
    undef $/;
    my $buf = <$FH3>;
    $/ = $holdTerminator;
    my @lines = split /$holdTerminator/, $buf;
    $buf = "init";
    $buf = join $holdTerminator, @lines;
    $buf =~ s/>ENST\d//g;
    $buf =~ s/\d//g;
    $buf =~ s/\n//g;
    $buf =~ s/_//g;
    $buf =~ s/($match[0])//g;
    $buf =~ s/($match[1])//g;
    @modi = parse($buf);
#    print "\n";
    return;
}

### Sub to print just the variation in the repeat regions between the normal and modified protein sequence
### Arg: the parsed repeat regions from the modified and reference sequence

sub short_repeat{

    my($refer,$modify) = @_;
    my $mismatch_index = 0;
    my $count = 0;

    my @r = split(//,$refer);
    my @m = split(//,$modify);

#    print "\n @r \n @m\n";
    #print "\n $refer,  $modify \n";
### Checks if the parsed repeat sequence is longer than 10 characters
### Identifies the mismatch locations and the number of mismatches

   if(scalar(@r)>10){
	if(scalar(@r)>scalar(@m)){
	    for(my $i =0; $i <scalar(@m); ++$i){
		if($r[$i]!=$m[$i]){
		    $mismatch_index = $i;
		    $count = $count+1;
		}
	    }

	}
	else{
	    for(my $i =0; $i <scalar(@r); ++$i){
		if($r[$i]!=$m[$i]){
		    $mismatch_index = $i;
		    $count = $count+1;
		}
	    }

	}
	
### Checks if there is only one mismatch
### Prints only the region surrounding the mismatch
	
	if($count == 1){
	    if($mismatch_index >= 6){
		if($r[$mismatch_index-6] eq ')'){
		    print $r[$mismatch_index-8];
		    print $r[$mismatch_index-7];
		    print $r[$mismatch_index-6];
		}
		elsif($r[$mismatch_index-6] !~ /\d/g){
		    print $r[$mismatch_index-6];
		}
		for(my $j=$mismatch_index-5; $j<$mismatch_index+3; ++$j){
		    print $r[$j];
		}

		if($r[$mismatch_index+3] !~ /[(]/g){
		    print $r[$mismatch_index+3];
		    if($r[$mismatch_index+3] eq ')'){
			print $r[$mismatch_index+4];
		    }
		}
		print "\t";
		if($m[$mismatch_index-6] eq ')'){
		    print $m[$mismatch_index-8];
		    print $m[$mismatch_index-7];
		    print $m[$mismatch_index-6];
		}
		elsif($m[$mismatch_index-6] !~ /\d/g){
		    print $m[$mismatch_index-6];
		}
		for(my $j=$mismatch_index-5; $j<$mismatch_index+3; ++$j){
		    print $m[$j];
		}
		
		if($m[$mismatch_index+3] !~ /[(]/g){
		    print $m[$mismatch_index+3];
		    if($m[$mismatch_index+3] eq ')'){
			print $m[$mismatch_index+4];
		    }
		}
		print "\n";
	    }
	    else{
		for(my $j=0; $j<$mismatch_index+3; ++$j){
		    print $r[$j];
		}
		if($r[$mismatch_index+3] !~ /[(]/g){
		    print $r[$mismatch_index+3];
		    if($r[$mismatch_index+3] eq ')'){
			print $r[$mismatch_index+4];
		    }
		}
		print "\t";
		for(my $j=0; $j<$mismatch_index+3; ++$j){
		    print $m[$j];
		}
		if($m[$mismatch_index+3] !~ /[(]/g){
		    print $m[$mismatch_index+3];
		    if($m[$mismatch_index+3] eq ')'){
			print $m[$mismatch_index+4];
		    }
		}
		print "\n";
	    }
	}
	else{
	    print @r,"\t",@m,"\n";
	}
   }
    else{
	print @r,"\t",@m,"\n";
    }
}
