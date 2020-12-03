#!/usr/local/bin/perl -w

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use strict;
use DBI;

my $db = 
   DBI->connect( "dbi:SQLite:dbname=/opt/www/sift/db/Human_db/TRANSCRIPT_DB.sqlite",
         "",  "", { RaiseError => 1, AutoCommit => 1} );
my $sth_enst = $db->prepare("select * from ENST_REGION where ENST = ? order by rowid");

my $enst = $ARGV[0];
my $indel_start = $ARGV[1];
my $indel_stop = $ARGV[2];
my $enst_content;
my @result_row;
my @rows;
$sth_enst->execute("$enst");
while (@result_row = $sth_enst->fetchrow_array){
	push @rows, join("\t", @result_row);
	#print join("\t", @result_row),"\n";

}
my $outstr = "";
my $row;
my $num_exon = 0;
my $num_intron = 0;
my $num_3UTR = 0;
my $num_5UTR = 0;
my $CDS_length = 0;
my $tx_length = 0;
my $protein_length = 0;
my $indel_start_region = "";
my $indel_start_region_ctr = 0;
my $indel_stop_region = "";
my $indel_stop_region_ctr = 0;
my $indel_start_tx = 0;
my $indel_stop_tx = 0;
my $indel_start_CDS = 0;
my $indel_stop_CDS = 0;
my $indel_loc_pct = 0;
my $visual_tx = "---";
my $last_exon_length = 0;
my $last_exon_seen = 0;
our $orn;
foreach $row (@rows){
	chomp $row;
	my @elts = split /\t/, $row;
	my $chr = $elts[1];
	my $enst = $elts[0];
	my $region = $elts[5];
	my $begin = $elts[2];
	my $end =  $elts[3];
	my $length = $end - $begin;
	$tx_length += $length;
	my $region_ctr = $elts[6];
	$orn = $elts[4];
	my $gene_id = $elts[7];
	my $exon_od = $elts[8];
	
	if ($region =~ /CDS/i){
		$num_exon++ ;
		$CDS_length += $length;
		$visual_tx.="[]";
		if ($orn eq "1" ){
			$last_exon_length = $length;
		}
		else{
			if($last_exon_seen == 0){
				$last_exon_length = $length;
	                        $last_exon_seen = 1;
			}	
		}
		
	}
	if ($region =~ /intron/i){
		$num_intron++;
		$visual_tx.="--";
	}
	if ($region =~ /5UTR/i){
		$num_5UTR++;
		$visual_tx.="{}";
		
	}
	if ($region =~ /3UTR/i){
		$num_3UTR++;
		$visual_tx.="{}";
	}
	if ($indel_start >= $begin && $indel_start < $end){
		$indel_start_region = $region;
		$indel_start_region_ctr = $region_ctr;
		my $visual_tx_length = length($visual_tx);
		if ($orn eq "1"){
			substr ($visual_tx,$visual_tx_length-1,0) = ".";
		}
		else{
			substr ($visual_tx,$visual_tx_length-1,0) = "*";
		}
		if ($indel_start_region =~ /CDS/i){
			$indel_start_CDS = $CDS_length - ($end - $indel_start);
		}
		else{
			$indel_start_CDS = $CDS_length;
		}
		$indel_start_tx = $tx_length - ($end - $indel_start);
        }
	if ($indel_stop > $begin && $indel_stop <= $end){
                $indel_stop_region = $region;
                $indel_stop_region_ctr = $region_ctr;
		my $visual_tx_length = length($visual_tx);
		if ($orn eq "1"){
			substr ($visual_tx,$visual_tx_length-1,0) = "*";
		}
		else{
			substr ($visual_tx,$visual_tx_length-1,0) = ".";
		}

		if ($indel_stop_region =~ /CDS/i){
                        $indel_stop_CDS = $CDS_length - ($end - $indel_stop);
                }
		else{
			$indel_stop_CDS = $CDS_length;
		}
		$indel_stop_tx = $tx_length - ($end - $indel_stop);

        }


	#print "$chr\t$enst\t$region\t$begin\t$end\n";
}
if ($orn eq "1"){
	substr($visual_tx,0,0) = "|";
	my $visual_tx_length = length $visual_tx;
	substr($visual_tx,$visual_tx_length,0) = "--->";
	$indel_loc_pct = sprintf ("%.0f",$indel_start_CDS*100/$CDS_length);
	
}
else{
	substr($visual_tx ,0,0) = "<";
	my $visual_tx_length = length $visual_tx;
        substr($visual_tx,$visual_tx_length,0) = "---|";
	if ($indel_start_CDS != 0){
		$indel_start_CDS = $CDS_length - $indel_start_CDS;
	}
	if ($indel_stop_CDS != 0){
                $indel_stop_CDS = $CDS_length - $indel_stop_CDS;
        }
	if ($indel_start_tx != 0){
                $indel_start_tx = $tx_length - $indel_start_tx;
        }
        if ($indel_stop_tx != 0){
                $indel_stop_tx = $tx_length - $indel_stop_tx;
        }

	($indel_start_CDS,$indel_stop_CDS) = ($indel_stop_CDS,$indel_start_CDS);
	($indel_start_tx,$indel_stop_tx) = ($indel_stop_tx,$indel_start_tx);	
	$indel_loc_pct = sprintf("%.0f",$indel_start_CDS*100/$CDS_length);

}
$protein_length = $CDS_length/3;
#print "5UTR: $num_5UTR\nEXONS: $num_exon\nINTRONS: $num_intron\n3UTR: $num_3UTR\nTranscript Length: $tx_length\nCDS Length: $CDS_length\nProtein Length: $protein_length\n";

#print "Indel spans $indel_start_region $indel_start_region_ctr to $indel_stop_region $indel_stop_region_ctr\n";
#print "$visual_tx\n";
#print "$indel_start_tx\t$indel_stop_tx\t$indel_loc_pct\%\n";
#print "$indel_start_CDS\t$indel_stop_CDS\t$indel_loc_pct\%\n";
$outstr = "$indel_start_CDS\t$indel_stop_CDS\t$indel_loc_pct\%\t$indel_start_region $indel_start_region_ctr to $indel_stop_region $indel_stop_region_ctr\t$num_5UTR\t$num_exon\t$num_intron\t$num_3UTR\t$tx_length\t$CDS_length\t$protein_length\t$last_exon_length\t$visual_tx";
print $outstr;
