#!/usr/bin/perl
#	Copyright 2003-2007 Fred Hutchinson Cancer Research Center
#	Tree viewer

#>>>>points to problems/issues
#>>>>initial $seqfile has all seqs in .pros, not just those in blocks
 
# 2003 Nick Taylor
# 2/28/07 Moved from Proweb server and modified for Blocks Server

#16112_1.tre has initial tree
#treeviewer.pl?file=../tmp/16112.tre&family=IPB000104&hash=../tmp/16112.hash&lengths=n
#Links from tree nodes:
#treeviewer.pl?file=16112&family=IPB000104&hash=16112&action=above|below
#	&contains=<seqnames>
#Clicking results in 16145.tre
# treeviewer.pl?file=16145&family=IPB000104&hash=16112&lengths=n

#Temporary files, some must be addressable from ../tmp and from /tmp:
# $$.tre	Tree file	../tmp and /tmp
# $$.png	Image file	../tmp and /tmp
# $$.hash			../tmp
# $$.blks	Blocks file	../tmp
# $$.seqs	Sequences	../tmp
# $$.cmd			../tmp
# $$.err 			../tmp

# 2/21/07 Moved from proweb. Dropped mysql and chromdb stuff.

use strict;
use CGI qw/:standard/;
use CGI::Carp;
use Storable;

#>>>uncomment for <PRE> 
#print "Content-type: text/html\n\n";
#print "<HTML><TITLE>Blocks Treeviewer</TITLE>\n";
#**********************************

#	Global variables
my($helpurl,$treeviewer,$expasyurl,$jpegquality,$tmpdir,$tmpfile,$htmlurldir,$blocksdir,$printsdir);

# URL of help page
$helpurl	=	"/treeviewer/info.html";

# Location of treeviewer binary (C program using libgd)
$treeviewer	=	"./treeviewer";

# Quality of JPEG image, leave blank for default
$jpegquality	=	30; 

$tmpdir = "../tmp";
$htmlurldir = "/tmp";		# must be soft link www/tmp -> ../tmp
$tmpfile = "$$";
$blocksdir = "../data-blocks";
$printsdir = "../data-prints";
$expasyurl = "http://us.expasy.ch/cgi-bin/get-sprot-entry?";


#--------------------------------------------------------------------------

treeviewer();

sub fatal 
{
	# Handle an error gracefully.
	my $reason = shift;
	print header();
	print "<HTML><HEAD><TITLE>Fatal Error</TITLE></HEAD><BODY BGCOLOR=white>\n";
	print "<P>Processing could not continue because: $reason</P></BODY></HTML>";
	exit(-1);
} # end of fatal


sub treeviewer 
{
	my $nameshash;	# Holds SwissProt ID numbers for SP names in tree

	#Read input from form
	my $from = param('from');

	my $family = param('family');		# blocks family name
	$family = '' unless ($family);

	my $tree = param('tree');		# user-provided tree
	$tree = '' unless ($tree);
	# If the tree has been provided, ignore family name
	if ($tree =~ /\S/)	{ $family = ''; } 
	else 			{ $tree = ''; }

	#  full name of the previous tree file, e.g. ../tmp/12065.tre
	my $file = param('file');		# previous tree file name
	$file = '' unless ($file);
	unless ($file) {
		if ($family =~ /^[A-Za-z]+\d+$/) {
			# Only do this from old-style blocks family names
			$family =~ tr/a-z/A-Z/;
		}
		$family =~ s/\s//g;
	}
	if ($file && (! -e "$file")) 
	{ fatal("stored tree file not found. $file."); }

	#  full name of the previous hash file, e.g. ../tmp/12065.hash
	#  never changes from the initial hash file name
	my $hash = param('hash');		# previous hash file
	$hash = '' unless ($hash);
	if ($hash && (! -e "$hash")) 
	{ fatal("stored data file not found. $hash."); }

	my $action = param('action');		# for $treeviewer
	$action = '' unless ($action);

	my $contains = param('contains');	# list of seq names 
	$contains = '' unless ($contains);

	my $imgfile = "$tmpfile.png";		# new image
	my $format = param('format');
	$format = 'png' unless ($format);
	if ($format =~ /^j/i) { $imgfile =~ s/\.png$/\.jpg/; }

	my $lengths = param('lengths');
	$lengths = 0 unless ($lengths);

	#	These are set later
	my ($trefile, $blkfile, $seqfile, $blkflag, $seqflag, $initname);

#print "<PRE>from=$from file=$file hash=$hash tmpfile=$tmpfile\n</PRE>";
#<<<testing, don't want to loop>>>
#if ($from ne "start") { exit(0); }

	# Has user supplied blocks and sequences?
	# If not, and $family is specified get corresponding blocks & seqs
	#	if those files have non-zero size
	#	These files are not passed on to subsequent sessions,
	#	but the hash file is - it includes $blkflag and $seqflag
	#	So each time are starting with initial blks and seqs
        $blkflag = $seqflag = 0;
	$blkfile = "$tmpdir/$tmpfile.blks";
	$seqfile = "$tmpdir/$tmpfile.seqs";
	if (param('blocks'))
	{
		$blkflag = 1;
		open(BLK, ">$blkfile");
		{  print BLK param('blocks'); }
		close(BLK);
	}

	if (param('seqs')) 
	{
		$seqflag = 1;
		open(SEQ, ">$seqfile");
		{  print SEQ param('seqs'); }
		close(SEQ);
	}

	# If this is the initial session and family is specified,
	# get those blocks and seqs
	if (! ($file && (-e "$file")))   # initial session check
        {
           if (! $blkflag)
           {
              if (($family =~ m/^IPB/) && (-s "$blocksdir/blks/$family.blks"))
              {
		$blkflag = 1;
		system("cp $blocksdir/blks/$family.blks $blkfile");
              }
              elsif (($family =~ m/^PR/) && (-s "$printsdir/blks/$family.blks"))
              {
		$blkflag = 1;
		system("cp $printsdir/blks/$family.blks $blkfile");
              }
           }
           if (! $seqflag)
           {
              if (($family =~ m/^IPB/) && (-s "$blocksdir/pros/$family.pros"))
              {
		$seqflag = 1;
		system("cp $blocksdir/pros/$family.pros $seqfile");
              }
              elsif (($family =~ m/^PR/) && (-s "$printsdir/pros/$family.pros"))
              {
		$seqflag = 1;
		system("cp $printsdir/pros/$family.pros $seqfile");
              }
#>>>at this point $seqfile has all seqs in .pros, not just those in blocks
           }
        } # end of initial session check

	#---------------------------------------------------------------
	#	Have all the parameters, start processing

	# Initial session, get the tree from the
	# db and store the hash of swissprot names -> swissprot ID's
	if (! ($file && (-e "$file"))) 
	{
		$trefile = "$tmpfile.tre";
 		$file = "$tmpdir/$trefile";
		$nameshash = getTree($family, "$file", $tree);
 		$hash = "$tmpdir/$tmpfile.hash";	# initial hash file
		store([$nameshash, $blkflag, $seqflag], "$hash");
	} 
	# Continuation of a previous session, 
#>>>>but $file and $hash but not $contains may be if select 
#>>>>"View this tree with branch lengths turned on|off"<<<
	# $file, $hash, $contains better be set
	# The name of the hash file is the initial name
	# retrieve the hash of names -> ID's from file.
	else 
	{
		($trefile) = $file =~ m/$tmpdir\/(\S+)/;
		($initname) = $hash =~ m/(\S+)\.hash/;
		($nameshash, $blkflag, $seqflag) = @{retrieve($hash)};
		# The blocks and seqs files are only created once
 		{ ($blkflag, $seqflag) = prune($initname, $contains, $blkfile, $seqfile); }
	}

	#----------------------------------------------------------------
	# Render tree using program, then display as a HTML map. 
	# $map holds HTML code for image along
	# with associated map, $list is HTML code for list of entries in tree,
	# $bootstrap is bootstrap value for root node.
        # $trefile and $imgfile are just the file names, w/o directory
	# $trefile has the complete tree from the initial session

	my ($map, $links, $list, $bootstrap) = drawTree($trefile, $imgfile, $action, $contains, $hash, $nameshash, $format, $lengths, $family);

	print header();
	print '<HTML><HEAD><TITLE>', $family ? $family : 'User-supplied', ' Tree';
	if ($bootstrap) { print " (Bootstrap Value $bootstrap)"; }

	print "</TITLE><BODY BGCOLOR=white>\n";
	print '<H2 ALIGN=CENTER>', $family ? $family : 'User-supplied Tree', "</H2>\n";
	print "<P>For help, visit the <A HREF=$helpurl TARGET=\"TreeHelp\">Tree Viewer help page</A>.</P>\n";
	print $map;
	print $links;

	print '<P>';
	if ($blkflag) 
        { print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/format_blocks.pl?$$\">Format blocks for this tree</A>] "; }
	if ($seqflag) 
#       { print " [<A TARGET=\"_blank\" HREF=\"/blocks-bin/make_blocks.sh?$seqfile\">Make blocks from sequences in this tree</A>]"; }
        { print " [<A TARGET=\"_blank\" HREF=\"/blocks-bin/make_blocks.sh?$seqfile\">Make blocks from sequences in this group</A>]"; }
	print "</P>\n";

	print $list;

	print "</BODY></HTML>\n";
} # end of treeviewer

#-------------------------------------------------------------------------
#  Draws the initial tree, or extracts a sub-tree into a new
#  treefile and draws it
#
sub drawTree 
{
	my $infile = shift;		# tree file w/o directory
	my $outfile = shift;		# image file w/o directory
	my $action = shift;		# Action to take on specified subtree
	my $contains = shift;		# Elements in specified subtree as CSV
	my $hash = shift;		# Full name of hash file ../tmp/111.hash
	my $nameshash = shift;		# Ptr to names->ids hash
	my $format = shift;		# Graphics format, begins w/ 'j' for JPEG
	my $lengths = shift;		# Use lengths? begins with 'n' for no
	my $family = shift;		# Family name
#print "<PRE>drawTree: infile=$infile outfile=$outfile action=$action\n</PRE>";
	my $tmpin = "$tmpdir/$infile";
	my $tmpout = "$tmpdir/$outfile";
 	my $errfile = "$tmpdir/$tmpfile.err";
 	my $cmdfile = "$tmpdir/$tmpfile.cmd";
 	my $subtre = "$tmpdir/$tmpfile.tre";
	my $options = '';
	my $flags = '';
	my @tree;			# All of the nodes in the tree
	my @seqs;			# All of the leaves (seqs) in the tree
	my ($map, $links);
	my ($textspacing, $boxspacing, $fontwidth, $fontheight, $clicksize, $width, $height, $bootstrap);
	my (@line,$size,$head,$url,$name,$seqs,@sprot_ids,$sprot_ids,$list);

	
	if ($format =~ /^j/i) {		# Make JPEG format
		$options = '-jpeg ';
		$options .= $jpegquality.' ' if defined($jpegquality);
		$flags .= '&format=j';
	}
	
	if ($lengths =~ /^n/i) {	# Ignore branch lengths
		$options .= ' -nolengths';
		$flags .= '&lengths=n';
	}
	
	fatal("Treeviewer binary not installed.") unless ($treeviewer);
	fatal("Treeviewer binary $treeviewer not found.") unless (-x $treeviewer);

	if ($action) {			# Make commands file if the tree
					# renderer must extract a subtree
					# or reroot the tree.
		open(COMMANDS, ">$cmdfile") || fatal("Could not create commands file, $!.");
		print COMMANDS $action, "\t", $contains, "\n";
		#  Write out the sub-tree file
 		$subtre = "$tmpdir/$tmpfile.tre";
		print COMMANDS 'write', "\t", $subtre, "\n";
		close(COMMANDS);

		open(VIEWER, "$treeviewer $tmpin $tmpout $cmdfile $options 2>$errfile |"); 
		# Change $infile so that processes
		# started from the new tree refer to the newly created tree
		$infile = "$tmpdir/$tmpfile.tre";
	} else {
					# No action to be taken, just render
					# the tree as is.
		$subtre = $tmpin;	# subtree is the whole tree
		open(VIEWER, "$treeviewer $tmpin $tmpout $options 2> $errfile |");
	}
#print "<PRE>drawTree: open VIEWER\n</PRE>";
	
	while (<VIEWER>) {		# Extract tree rendering information
					# from the output of the renderer.
		if (/^TextSpacing\s+(\d+)/) {
			$textspacing = $1;
		} elsif (/^BoxSpacing\s+(\d+)/) {
			$boxspacing = $1;
		} elsif (/^FontWidth\s+(\d+)/) {
			$fontwidth = $1;
		} elsif (/^FontHeight\s+(\d+)/) {
			$fontheight = $1;
		} elsif (/^ClickSize\s+(\d+)/) {
			$clicksize = $1;
		} elsif (/^Width\s+(\d+)/) {
			$width = $1;
		} elsif (/^Height\s+(\d+)/) {
			$height = $1;
		} elsif (/^Bootstrap\s+(\d+)/) {
			$bootstrap = $1;
		} elsif (/^Tree/) {
			@line = split(/\s+/);
			push(@tree, [ @line[1..3] ]);
		}
	} # end of treeviewer output
	close(VIEWER);

	if (-s "$errfile") {
		system("cat $errfile");
		fatal("tree rendering failed.");
	}

	$map = '<MAP NAME="tree">';	# Begin to draw the map of clickable
					# nodes.
	$size = int($clicksize / 2);
	$head = pop(@tree);		# Ignore the root node if the tree is a
					# not a leaf
	push(@tree, $head) if ($#tree < 0);
	foreach (@tree) {		# For each node in the tree, place a
					# clickable node over the boxes on the
					# graphic.
		if ($_->[2] =~ /,/) {	# Subtree node, allow zoom/reroot.
			$map .= '<AREA SHAPE=rect COORDS="';
			$map .= ($_->[0] - $size - 1).','.($_->[1] - $size).',';
			$map .= ($_->[0] - 1).','.($_->[1] + $size).'" ';
			$map .= ' TARGET="_blank" HREF="'.url().'?file='.$subtre.'&hash='.$hash.$flags.'&action=above&family='.$family.'&contains='.$_->[2]."\"\>\n";
			$map .= '<AREA SHAPE=rect COORDS="';
			$map .= ($_->[0]).','.($_->[1] - $size).',';
			$map .= ($_->[0] + $size).','.($_->[1] + $size).'" ';
			$map .= ' TARGET="_blank" HREF="'.url().'?file='.$subtre.'&hash='.$hash.$flags.'&action=below&family='.$family.'&contains='.$_->[2]."\"\>\n";
		} else {		# Leaf, allow link to SwissProt.
			if ($_->[2] =~ /^\w+\$(.*)$/) {
				$name = $1;
			} elsif ($_->[2] =~ /^[A-Z0-9_]+$/) {
				$url = $expasyurl;
				if ($nameshash->{$_->[2]}) {
					$url .= $nameshash->{$_->[2]};
				} else {
					$url .= $_->[2];
				}
				$name = $_->[2];
			} else {
				$name = $_->[2];
			}

			push(@seqs, $name);

			if ($url) {
				$map .= '<AREA SHAPE=rect COORDS="';
				$map .= ($_->[0] + $textspacing - $boxspacing).','.($_->[1] - int($fontheight / 2) - $boxspacing).',';
				$map .= ($_->[0] + $textspacing + $boxspacing + $fontwidth * length($name)).','.($_->[1] + int($fontheight / 2) + $boxspacing)."'";
				$map .= qq{ TARGET="_blank" HREF="$url">\n};
			}
		} # end of leaf
	} # end of tree node
	$map .= "</MAP>\n<IMG BORDER=0 SRC=\"".$htmlurldir.'/'.$outfile."\" USEMAP=\"#tree\" WIDTH=$width HEIGHT=$height>\n";
	$seqs = join(' ', @seqs);
					# List to pass to prune.cgi
	@sprot_ids;			# List to pass to blockmkr.cgi
	
	foreach (@seqs) {
		if ($nameshash->{$_}) {
			push(@sprot_ids, $nameshash->{$_});
		} else {
			push(@sprot_ids, $_);
		}
	}
	$sprot_ids = join(' ', @sprot_ids);

#>>>>don't get blkflag and seqflag links after select one of these<<<<
	if ($lengths =~ /^n/i) {	# Button to draw w/ branch lengths.
		$links = "<P>[<A HREF=\"".url()."?file=$subtre&family=$family&hash=$hash";
		$links .= '&format=j' if ($format =~ /^j/i);
		$links .= '">View this tree with branch lengths turned on</A>]';
	} else {			# Button to draw w/out lengths
		$links = "<P>[<A HREF=\"".url()."?file=$subtre&family=$family&hash=$hash&lengths=n";
		$links .= '&format=j' if ($format =~ /^j/i);
		$links .= '">View this tree with branch lengths turned off</A>] ';
	}
					# Other buttons to process this
					# sequences in this tree
	$links .= " [<A TARGET=\"_blank\" HREF=\"$htmlurldir/$infile\">Newick formatted tree file</A>]</P>\n";

	$list = join("<BR>\n", @seqs);	# HTML list of sequences in this tree.
	$list = "<P>There are ".($#seqs + 1)." sequences in this tree:<BR>".$list."</P>\n";

	return ($map, $links, $list, $bootstrap);
} #end of drawTree

sub getTree 
{
	# Retrieve a tree from the Blocks DB
	my $family = shift;		# Family to retrieve
	my $outfile = shift;		# File into which to put the tree
	my $newick = shift;
	my $nameshash = {};
	my ($file, $famfile, $gottree);

        $famfile = "";
        if    ($family =~ m/^IPB/) { $famfile = "$blocksdir/trees/$family.tree" }
        elsif ($family =~ m/^PR/)  { $famfile = "$printsdir/trees/$family.tree" }
#print "<PRE>getTree: family=$family outfile=$outfile famfile=$famfile\n</PRE>";

	if (! $newick ) 
        {
		if ($family =~ /^([A-Z]+)\d+$/)
		{ $file = $1; }
		elsif ($family =~ /^([A-Z][a-z]+)[A-Z]+/)
		{ $file = $1; }
			
		open(BLOCKS, $famfile) || fatal ("Could not open Blocks DB file, $!");

		$gottree = 0;
		while (<BLOCKS>) 
		{
				if (/^>/ && $gottree) {
					last;
				}
				if (/^>$family /) {
					$gottree = 1;
					next;
				}
				if ($gottree) {
					$newick .= $_;
				}
		}
		close(BLOCKS);
                if ($gottree == 0)
                {
                   print "<B>ERROR</B> Invalid $famfile. <P>\n";
                   exit(-1);
                }
		
	} # end of if not newick

	$newick =~ s/([^\s,\(\)\;:]+)/&processName($1,$nameshash)/eg;

	open(OUT, ">$outfile") || fatal("could not create tree file, $!");
	print OUT $newick;
	close(OUT);
	return($nameshash);
} # end of getTree

sub processName 
{
	# Add an entry to an anonymous hash, "A,B" adds $hash->{A} = B.
	my $name = shift;		# Input to split
	my $hash = shift;		# Hash to add to.
	return($name) unless ($name =~ /\|/);
	my @name = split(/\|/,$name);
	if ($#name == 2) {
		$hash->{$name[0]} = $name[1];
	}
	return($name[0]);
} # end of processName

#-------------------------------------------------------------------------
#	Reduces initial blkfile and/or seqfile to subset of sequences in 
#	$contains and puts them in current $blkfile and/or $seqfile
sub prune
{
	my ($initname, $contains, $blkfile, $seqfile) = @_;
	my ($bflag, $sflag, @snames, $in_rec, $out_rec, $bname, $sname, $keep);

#print "<PRE>prune: initname=$initname contains=$contains\n</PRE>";

	$bflag = $sflag = 0;
	@snames = split(/,/, $contains);
	if ($#snames >= 0)
        {
	   if (-s ("$initname.blks"))
           {
		open(IN, "<$initname.blks");
		open(OUT, ">$blkfile");
		foreach $in_rec (<IN>)
                {
		  if ($in_rec =~ m/^ID   /) { print OUT "$in_rec"; }
		  if ($in_rec =~ m/^AC   /) { print OUT "$in_rec"; }
		  if ($in_rec =~ m/^DE   /) { print OUT "$in_rec"; }
		  if ($in_rec =~ m/^BL   /) { print OUT "$in_rec"; }
		  if ($in_rec =~ m/^\/\//) { print OUT "$in_rec"; }
		  else
                  {
			($bname) = $in_rec =~ m/(\S+) /;
			#	seqname in block is generally longer
			foreach $sname (@snames)
                        {
			   if ($sname =~ m/$bname/) { print OUT "$in_rec"; }
                        }
                  }
 		}
		close(OUT);
		close(IN);
		$bflag = 1;
           }
	   if (-s ("$initname.seqs"))
           {
		open(IN, "<$initname.seqs");
		open(OUT, ">$seqfile");
		#  Assumes seqs are in fasta format
		$in_rec = <IN>;		# first record
		while (!eof(IN))
                {
		   if ($in_rec =~ m/^>/)
		   {
			$keep = 0;
		   	($bname) = $in_rec =~ m/>(\S+) /;
			foreach $sname (@snames)
                        {
			   if ($sname =~ m/$bname/) 
			   { 
				$keep = 1;
				print OUT "$in_rec"; 
				while (($in_rec = <IN>) &&
                                       !($in_rec =~ m/^>/))
				{  print OUT "$in_rec";  }
			   }
                        }
		        if ($keep == 0) { $in_rec = <IN>; }
                   }
                   else { $in_rec = <IN>;  }
		} # end of seqs file
		close(OUT);
		close(IN);
		$sflag = 1;
	   }
        } # end of non-empty list of sequence names
	return($bflag, $sflag);
} # end of prune
