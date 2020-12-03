#!/usr/bin/perl
# (C) Copyright 2001-2003, Fred Hutchinson Cancer Research Center
#	nixcheck(IP address)
#	Compare IP address against file of IPs we don't want to
#	accept traffic from

#------------------------------------------------------------------------

if (@ARGV < 1) { exit(-1); }
else { $ipaddr = @ARGV[0];  }

$nix = "../nixlist";

$done = 0;
if (-e $nix)
{
   open(NIX, "<$nix");
   #   The most current record is the last one
   while ( ($done == 0) && ($_ = <NIX>) )
   {    
      @words = split(/\s+/, $_); 
      $nixip = $words[0];
#     $nixdate = $words[1];
   #		nixlist may have 12.107.64, etc.
#      if ($nixip eq $ipaddr) { $done = 1; }
       if ($ipaddr =~ m/^$nixip/)  { $done = 1; }
   }
   close(NIX);
}
if ($done == 1)
{
   print "Blocks is currently not accepting requests from $ipaddr.\n";
}

exit($done);
