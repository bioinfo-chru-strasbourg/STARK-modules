#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

close(STDIN);
close(STDOUT);
close(STDERR);
exit if (fork());
exit if (fork());
system("/opt/www/sift/htdocs/sift-bin/start_queue.sh");


sub daemonize {
   use POSIX;
   POSIX::setsid or die "setsid: $!";
   my $pid = fork ();
   if ($pid < 0) {
      die "fork: $!";
   } elsif ($pid) {
      exit 0;
   }
   chdir "/";
   umask 0;
   foreach (0 .. (POSIX::sysconf (&POSIX::_SC_OPEN_MAX) || 1024))
      { POSIX::close $_ }
   open (STDIN, "</dev/null");
   open (STDOUT, ">/dev/null");
   open (STDERR, ">&STDOUT");
 }

