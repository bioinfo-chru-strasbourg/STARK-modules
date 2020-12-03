#!/bin/csh
#
#		logo_nongrouped_web.csh <pid> <type>
#		<type> = ps | pdf | gif
#	Makes postscript logos & converts them to individual
# pdf or gif format
#  The files colors, marks, wave and default.amino.frq must be in
#  the current directory
#	

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

set pid = $1
set tmp = "../tmp"
setenv GS_LIB /opt/sfw/share/ghostscript
setenv LD_LIBRARY_PATH /usr/openwin/lib:/usr/dt/lib:/usr/lib:/usr/ucblib:/opt/sf
w/lib:/usr/local/lib
set path = ($path /opt/sfw/bin)

# Assume all logos have been made and post-processed 
# and the postcript files exist

#-------------------------------------------------------------------------
#	See if any logo files were created
#	If so, merge them all into one big postscript file
set logos = ($tmp/$pid.*.ps)

      foreach logo ($logos)
#			      0-492 was the original logos
#	but my logos are further down on the page. 
# crop field means widthxheight+xoffset-yoffset.  my yoffset is set lower
# to get the logos. use ghostview to get the pixel #.
	convert -crop 648x300+0-200 ps:$logo pdf:$logo.pdf
	convert -crop 648x300+0-200 ps:$logo $logo.gif
      end


exit(0)
