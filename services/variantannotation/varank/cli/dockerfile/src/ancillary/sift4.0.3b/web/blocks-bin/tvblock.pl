#!/usr/bin/perl

# Formerly proweb's blockmkr.cgi

use strict;
use lib '.';

# Fill in the Block Maker form with a fasta file created from the requested
# sequences, or just return the fasta file. Sequence names are separated by
# whitespace in param('sequences').

use CGI qw/:standard/;
use CGI::Carp;
use LWP::Simple qw(get);		# Only import get, so we don't get
					# conflicts w/ CGI.

use vars qw/$blockMakerUrl $blockMakerServer $expasyUrl/;


################################################################################
# End users should not need to edit anything below this line.

# URL of BlockMaker form
$blockMakerUrl		= 	'/blocks/blockmkr/make_blocks.html';

# URL of ExPASY script to get FASTA file for SwissProt entry.
$expasyUrl 		= 	'http://us.expasy.ch/cgi-bin/get-sprot-fasta?';

($blockMakerServer)	=	($blockMakerUrl =~ m{^(http://[\w\.]+)/});

sub fatal {
	# Die in a way that gets along with the web server.
	print "Content-type: text/plain\n\n";
	print "Fatal error: $_[0]\n";
	exit;
};

blockmkr();

sub blockmkr 
{

	my $tmpdir = retConfig('tmpdir');
	my $prositeDir = retConfig('prositedirs');

	my @sequences = split(/\s+/, param('sequences'));
	my $forceExpasy = param('expasy');
	my $blocks = param('blocks');
	my $seqs_data = param('seqs_data');
	$blocks =~ s/\s//g;

	my $family;

	if ($blocks) {
		# Don't try to read from user file if coming from
		# Blocks DB
		undef $seqs_data;
	}
	elsif ($seqs_data)
	{
		unless ($seqs_data =~ /^[\d_]+$/)
		{
			fatal("Invalid value for 'seqs_data'.");
		}
		
		unless (-e "$tmpdir/$seqs_data.seqs")
		{
			fatal("User-supplied sequence data not found.");
		}
		
		$forceExpasy = 0;
	}

	if (! $seqs_data) {
		if ($blocks =~ /^([A-Za-z]+)\d+$/)
		{
			$family = $1;
			$family =~ tr/a-z/A-Z/;
		}
		elsif ($blocks =~ /^([A-Z][a-z]+)[A-Z]+$/)
		{
			$family = $1;
		}
		else
		{
			fatal("malformed Blocks family, $blocks.");
		}

		if ($prositeDir && (! $prositeDir->{$family}))
		{
			fatal("invalid blocks family '$family'");
		}
		
		if ((! $prositeDir) || (! -e "$prositeDir->{$family}/$blocks.pros"))
		{
			$forceExpasy = 1;
		}
	}
	
	my %need;
	my %got;			# Keep track of ID's already retrieved.
	my $fasta;


	foreach (@sequences) {
		$need{$_} = 1;
	}

	if ($#sequences < 0) {
		fatal("You must select some sequences.");
	}

	if ($forceExpasy) {
		my $id;
		foreach $id (@sequences) {	# Get sequences from expasy.
			my $response = get($expasyUrl.$id);
			unless ($response) {
				fatal("Could not contact SwissProt.");
			}
			next if ($response =~ /<HTML>/i);
						# If the page contains HTML, it means an
						# error, most likely sequence not found.
			my ($got_id) = ($response =~ /^>[a-z]+\|([A-Z0-9]+)/);
						# Which ID's sequence did we actually
						# get?
			unless ($got{$got_id}) {
				$fasta .= $response."\n";
				$got{$got_id} = 1;
			}
			delete $need{$id};		# Note that we don't need it any more.
		}
	} else {
		if ($seqs_data) {
			open(DATA, "$tmpdir/$seqs_data.seqs");
		} else {
			open(DATA, "$prositeDir->{$family}/$blocks.pros");
		}
		my $current;				# Are we saving the current entry?
		while (<DATA>) {
			if (/^>([\w\|]+)/) {
				$current = 0;
				my @ids = split(/\|/,$1);
				my $id;
				foreach $id (@ids) {
					if ($need{$id}) {
						delete $need{$id};
						$current = 1; 
					}
				}
			}
			if ($current) {
				$fasta .= $_;
			}
		}

		close(DATA);
	}
	
	if (param('Action') =~ /FASTA/) {
		# User only asked for FASTA file, we can stop now.
		print "Content-type: text/plain\n\n";
		my $notfound = join(' ', keys(%need));
		print ">Not Found: $notfound\n" if ($notfound);
		print $fasta;
		exit;
	}

	my $blockMakerForm = get($blockMakerUrl);

	unless ($blockMakerForm) {
		fatal("Blockmaker form could not be loaded.");
	}

	print "Content-type: text/html\n\n";

	$fasta =~ s/>/&gt;/g;


	# Go through the Blockmaker form, filling in as needed.

	my $done = 0;

	foreach (split(/\n/, $blockMakerForm)) {
		if (/(HREF|ACTION|SRC)="([^"]+)"/) {
			my $arg = $1;
			my $url = $2;
			s/$arg="$url"/$arg="$blockMakerServer$url"/ if ($url =~ /^\//);
		}
		if (/<\/TEXTAREA>/) {
			s/<\/TEXTAREA>/$fasta<\/TEXTAREA>/;
		}
		print $_;
		if (/<HR>/ && (! $done)) {
			$done = 1;
			my $notfound = join(', ', sort(keys(%need)));
			if ($notfound) {
				print "<FONT COLOR=red>The following sequences were not found in SwissProt/TrEMBL: $notfound.</FONT>";
			}
		}
		print "\n";
	}

} # end of blockmkr
