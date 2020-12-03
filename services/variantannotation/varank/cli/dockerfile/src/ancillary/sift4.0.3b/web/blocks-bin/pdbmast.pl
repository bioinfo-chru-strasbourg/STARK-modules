#!/usr/bin/perl
#		pdbmast.pl <blk2pdb entry> <last_ac>
#		(replaces pdbfile.pl)
#	Insert links for pdbmast.dat entries
#---------------------------------------------------------------------------
#  7/19/02 fix 3dblocks link
#  2/18/03 More elaborate 3motif link, etc.
#  2/19/03 New format for blk2pdb.dat entries
#	Each line has: 
# <ac> <pdb>[_<chain>] <all-blocks-flag> <blocks-in-order-flag>
#  3/ 9/03 Added MMDB links
#  4/13/07 Execute getblock.pl instead of getblock.sh
#---------------------------------------------------------------------------

print "Content-type: text/html\n\n";
print "<HTML>\n";

if ( @ARGV < 2 )
{
   print "USAGE: pdbmast.pl <blk2pdb entry> <last_ac>\n";
   exit(-1);
}
$pdbmast = "$ARGV[0]";
$lastac = "$ARGV[1]";
#print "pdbmast=$pdbfile\n";
#print "lastac=$lastac\n";
open(PDB, "<$pdbmast") || die("pdbfile.pl: Could not open $pdbfile\n");
#while ($pdb_rec = <PDB>)
#{ print "$pdb_rec"; }
#close(PDB);
#exit;

$ac = $lastblk = "";
while ($pdb_rec = <PDB>)
{
   ($pdb_ac, $pdb_pdb, $pdb_all, $pdb_order) = split(/\s+/, $pdb_rec);
   if ($pdb_pdb =~ m/\_/)
   {   ($pdb, $chain) = $pdb_pdb =~ m/(\S+)\_(\S*)/; }
   else
   { $pdb = $pdb_pdb; $chain = ""; }
   if ($ac eq "")
   {
      if ($pdb_rec =~ m/^IPB/)
      {   
         $ac = substr($pdb_rec, 0, 9);  
         if ($lastac =~ m/^${ac}/ && length($lastac) > 9)
         {  $lastblk = substr($lastac, 9, 1); }
      }
      else
      {   
         $ac = substr($pdb_rec, 0, 7);  
         if ($lastac =~ m/^${ac}/ && length($lastac) > 7)
         {  $lastblk = substr($lastac, 7, 1); }
      }
#print "ac=$ac lastblk=$lastblk\n";
      print "<TITLE>Structure Links for $ac</TITLE>\n";
      print "<H1>Structure Links for $ac</H1>\n";

      print "<UL>";
      print "<LI>Structures found by searching ";
      print "<A HREF=\"/blocks-bin/getblock.pl?$ac\">$ac</A>";
      print " blocks vs <A HREF=\"http:\/\/www.rcsb.org\/pdb\/\">PDB</A> with ";
      print "<A HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">MAST</A>";
      print "<P>\n";

      print "<TABLE VSPACE=10 CELLSPACING=5 WIDTH=100%>\n";
      print "<TR>";
      print "<TH>All Blocks in Hit</TH>\n";
      print "<TH>Blocks in Order</TH>\n";
      print "<TH><A TARGET=\"_blank\" HREF=\"http://www.proweb.org/3dblocks_intro.html\">3D Blocks</A></TH>\n";
      print "<TH><A TARGET=\"_blank\" HREF=\"http://www.rcsb.org/pdb/\">Protein Data Bank</A></TH>\n";
      print "<TH><A TARGET=\"_blank\" HREF=\"http://molbio.info.nih.gov/doc/mrus/mol_r_us.html\">Molecules R US</A></TH>\n";
      print "<TH><A TARGET=\"_blank\" HREF=\"http://www.ncbi.nih.gov/Structure/MMDB/mmdb.shtml\">Entrez Structure Database</A></TH>\n";
      print "</TR>";
   }

   if ($ac ne "")
   {
      print "<TR>";

      print "<TD ALIGN=CENTER>";
      if ($pdb_all == 1) { print "Yes"; }
      else               { print "No";  }
      print "</TD>\n";

      print "<TD ALIGN=CENTER>";
      if ($pdb_order == 1) { print "Yes"; }
      elsif ($pdb_all == 1 && $pdb_order == 0) { print "No";  }
      else { print "   ";  }
      print "</TD>\n";

      print "<TD ALIGN=CENTER>";
#	&mast=yes&CUTOFF=0.0001   (searches the cobbler seq with blast now)
      print "<A TARGET=\"_blank\" HREF=\"http:\/\/www.proweb.org\/proweb-bin\/3dblocks.cgi?AC=$ac&PDB=$pdb\">$pdb_pdb</A>\n";
      print "</TD>\n";

      print "<TD ALIGN=CENTER>";
      print "<A TARGET=\"_blank\" HREF=\"http:\/\/www.rcsb.org\/pdb\/cgi\/explore.cgi?job=graphics&pdbId=$pdb\">$pdb_pdb</A>\n";
      print "</TD>\n";

      print "<TD ALIGN=CENTER>";
      print "<A TARGET=\"_blank\" HREF=\"http:\/\/molbio.info.nih.gov\/cgi-bin\/moldraw?$pdb\">$pdb_pdb</A>\n";
      print "</TD>\n";

      print "<TD ALIGN=CENTER>";
      print "<A TARGET=\"_blank\" HREF=\"http:\/\/www.ncbi.nlm.nih.gov\/Structure\/mmdb\/mmdbsrv.cgi?uid=$pdb\">$pdb_pdb</A>\n";
      print "</TD>\n";

      print "</TR>";
   }
}   # end of pdbmast
close(PDB);

if ($ac eq "")
{
   print "ERROR: Problem reading $pdbmast\n";
   close(PDB);
   exit(-1);
}


print "</TABLE>\n";
print "<P><HR>\n";

#===============3motif stuff==============================================
print "<LI><P>";
print "<A TARGET=\"_blank\" HREF=\"http://motif.stanford.edu/3motif\">Go to 3Motif at Stanford</A> ";
#print "Requires <A HREF=\"http://www.mdli.com/downloads/free.html\">Chime</A> ";
print " (Requires <A TARGET=\"_blank\" HREF=\"http://www.mdli.com/chime\">Chime</A>)<BR>\n";
print " 3Motif finds and displays structures for each block individually:<P>\n";
print "<UL>";
if ($lastblk ne "")
{
   for ($blk = "A"; $blk le $lastblk; $blk++)
   {
print "<LI><A TARGET=\"_blank\" HREF=\"http://dna.stanford.edu/cgi-bin/3matrix/nph-3matrix-browser.cgi?is_first=TRUE&passed_value_1=NULL&passed_value_2=NULL&ematrix=&pdb_id=&block=${ac}${blk}\">${ac}${blk}</A>\n";
   }
}
else
{
print "<LI><A TARGET=\"_blank\" HREF=\"http://dna.stanford.edu/cgi-bin/3matrix/nph-3matrix-browser.cgi?is_first=TRUE&passed_value_1=NULL&passed_value_2=NULL&ematrix=&pdb_id=&block=${ac}\">${ac}</A>\n";
}
print "</UL>";

print "</UL>";
print "</HTML>\n";

exit(0);

