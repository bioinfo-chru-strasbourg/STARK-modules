#!/usr/bin/perl

$initname = "../tmp/24295";
$seqfile = "../tmp/24301.seqs";
@snames = ("Q547T1_PSEAM", "ANPY_PSEAM", "ANP4_PSEAM", "ANPX_PSEAM");

           if (-s ("$initname.seqs"))
           {
                open(IN, "<$initname.seqs");
                open(OUT, ">$seqfile");
                #  Assumes seqs are in fasta format
                $in_rec = <IN>;         # first record
                while (!eof(IN))
                {
print STDERR "while: $in_rec";
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
	}
exit(0);
